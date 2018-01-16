""" Utilities for isgm
 Best to keep these separate from the Class modules
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import pdb
import numpy as np
import warnings

from astropy import constants as const
from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky

from linetools.analysis import absline as ltaa
from linetools.isgm.abscomponent import AbsComponent
from linetools.abund.ions import name_to_ion, ion_to_name
from linetools import utils as ltu


def chk_components(components, chk_match=False, chk_A_none=False, tol=0.2*u.arcsec):
    """ Performs checks on a list of components

    Parameters
    ----------
    components : list
      list of AbsComponent objects
    chk_match : bool, optional
      if True, require that the components match in RA/DEC, Zion, Ej, A but not velocity
    chk_A_none : bool, optional
      if True, require that A *not* be set
    tol : Quantity, optional
      Tolerance on matching SkyCoord. Default is 0.2*u.arcsec
    """
    tests = True
    # List
    if not isinstance(components, list):
        tests = False
        raise IOError('Need a list of AbsComponent objects')
    # Object
    if not all(isinstance(x, AbsComponent) for x in components):
        tests = False
        raise IOError('List needs to contain only AbsComponent objects')
    # A None
    if chk_A_none:
        if not all(x.A is None for x in components):
            tests = False
            raise IOError('Not ready for components with A set')
    # Matching?
    if chk_match:
        match = True
        comp0 = components[0]
        for comp in components[1:]:
            # RA/DEC
            match = match & bool(comp0.coord.separation(comp.coord) < tol)
            # Zion
            match = match & (comp0.Zion == comp.Zion)
            # Ej
            match = match & np.allclose(comp0.Ej.to('1/cm').value,comp.Ej.to('1/cm').value)
            # A
            match = match & (comp0.A == comp.A)
        tests = tests & match
    # Return
    return tests


def build_components_from_abslines(iabslines, clmdict=None, coord=None,
                                   **kwargs):
    """ Generate a list of AbsComponent from a list of abslines

    Groups lines with like Zion, Ej, (and A; future)

    Parameters
    ----------
    abslines : list
      List of AbsLine objects
      May be ignored if clmdict is passed in
    clmdict : dict, optional
      If present, build the abslines list from this dict
    coord : SkyCoord, optional
      Required if clmdict is used

    Returns
    -------
    components :
      list of AbsComponent objects

    """
    if clmdict is None:
        abslines = iabslines
    else:
        raise DeprecationWarning("Gone")
        abslines = []
    # Test
    if not isinstance(abslines,list):
        raise IOError('Need a list of AbsLine objects')

    # Identify unique Zion, Ej combinations in the lines
    uZiE = np.array([iline.data['Z']*1000000+iline.data['ion']*10000+
                     iline.data['Ej'].to('1/cm').value for iline in abslines])
    uniZi, auidx = np.unique(uZiE, return_index=True)

    # Loop to build components
    components = []
    for uidx in auidx:
        # Synthesize lines with like Zion, Ej
        mtZiE = np.where(uZiE == uZiE[uidx])[0]
        lines = [abslines[ii] for ii in mtZiE] # Need a list
        # Generate component
        if lines[0].data['Ej'].value > 0.:
            # Grab stars from transition name
            nstars = lines[0].name.count('*')
            if nstars == 0:
                raise ValueError("Should have at least one *")
            stars = '*'*nstars
        else:
            stars = None
        component = AbsComponent.from_abslines(lines, stars=stars, **kwargs)
        # Reset vmin, vmax
        vmin,vmax = 9999., -9999.
        for iline in lines:
            vmin = min(vmin, iline.limits.vlim[0].value)
            vmax = max(vmax, iline.limits.vlim[1].value)
        component.limits.set([vmin,vmax]*u.km/u.s)
        # Append
        components.append(component)
    # Return
    return components


def build_components_from_dict(idict, coord=None, **kwargs):
    """ Generate a list of components from an input dict

    Parameters
    ----------
    idict : dict
      Must contain either components or lines as a key
    coord : SkyCoord, optional

    Returns
    -------
    components :
      list of AbsComponent objects
      Sorted by zcomp
    """
    from linetools.spectralline import AbsLine

    components = []
    if 'components' in idict.keys():
        # Components
        for key in idict['components']:
            components.append(AbsComponent.from_dict(idict['components'][key], coord=coord, **kwargs))
    elif 'lines' in idict.keys():  # to be deprecated
        lines = []
        for key in idict['lines']:
            if isinstance(idict['lines'][key], AbsLine):
                line = idict['lines'][key]
            elif isinstance(idict['lines'][key], dict):
                line = AbsLine.from_dict(idict['lines'][key], coord=coord, **kwargs)
            else:
                raise IOError("Need those lines")
            if coord is not None:
                line.attrib['coord'] = coord
            lines.append(line)
        components = build_components_from_abslines(lines, **kwargs)
    else:
        warnings.warn("No components in this dict")
    # Sort by z -- Deals with dict keys being random
    z = [comp.zcomp for comp in components]
    isrt = np.argsort(np.array(z))
    srt_comps = []
    for idx in isrt:
        srt_comps.append(components[idx])
    # Return
    return srt_comps


def build_systems_from_components(comps, systype=None, vsys=None, **kwargs):
    """ Build a list of AbsSystems from a list of AbsComponents
    Current default implementation allows for overlapping components, i.e.
      only_overlap=True in add_component

    Parameters
    ----------
    comps : list
        List of AbsComponents
    systype : AbsSystem, optional
      Defaults to GenericAbsSystem
    vsys : Quantity, optional
      'Velocity width' of a system, used when adding components
      Passed as vtoler to add_component
      The first component will define the system redshift and all others will
      need to lie within vsys of it
    **kwargs -- Passed to add_component()

    Returns
    -------
    abs_systems : list

    """
    # Suppress warnings here as we are expecting quite a few
    import warnings
    warnings.filterwarnings("ignore")

    if systype is None:
        from linetools.isgm.abssystem import GenericAbsSystem
        systype = GenericAbsSystem
    if vsys is None:
        if 'overlap_only' not in kwargs.keys():
            kwargs['overlap_only'] = True
    else:
        kwargs['vtoler'] = vsys.to('km/s').value
    # Add
    abs_systems = []
    cpy_comps = [comp.copy() for comp in comps]
    # Loop until all components assigned
    while len(cpy_comps) > 0:
        # Use the first one
        comp = cpy_comps.pop(0)
        abssys = systype.from_components([comp])
        # Try the rest
        comps_left = []
        for icomp in cpy_comps:
            if abssys.add_component(icomp, **kwargs):
                pass
            else:
                comps_left.append(icomp)
        # Update vlim
        abssys.update_vlim()
        # Append
        abs_systems.append(abssys)
        # Save
        cpy_comps = comps_left
    # Return
    return abs_systems


def xhtbl_from_components(components, ztbl=None, NHI_obj=None):
    """ Generate a Table of XH values from a list of components
    Parameters
    ----------
    components
    ztbl
    NHI_obj

    Returns
    -------

    """
    # Get started
    tbl = Table()
    #


def complist_from_table(table):
    """
    Returns a list of AbsComponents from an input astropy.Table.

    Parameters
    ----------
    table : Table
        Table with component information (each row must correspond
        to a component).
        Each column is expecting a unit when appropriate.
           Units for vmin, vmax are assumed km/s if not given
           Units for Ej are assumed cm^-1 if not given

    Returns
    -------
    complist : list
        List of AbsComponents defined from the input table.

    Notes
    -----
    Mandatory column names: 'RA', 'DEC', 'ion_name', 'z_comp', 'vmin', 'vmax'
        These column are required.
    Special column names: 'name', 'comment', 'logN', 'sig_logN', 'flag_logN'
        These columns will fill internal attributes when corresponding.
        In order to fill in the Ntuple attribute all three 'logN', 'sig_logN', 'flag_logN'
        must be present. For convenience 'logN' and 'sig_logN' are expected to be floats
        corresponding to their values in np.log10(1/cm^2).

    Other columns: 'any_column_name'
        These will be added as attributes within the AbsComponent.attrib dictionary,
        with their respective units if given.

    """
    # Convert to QTable to handle units in individual entries more easily
    tkeys = table.keys()

    # mandatory and optional columns
    min_columns = ['RA', 'DEC', 'ion_name', 'z_comp', 'vmin', 'vmax', 'Z', 'ion', 'Ej']
    special_columns = ['name', 'comment', 'logN', 'sig_logN', 'flag_logN']
    for colname in min_columns:
        if colname not in table.keys():
            raise IOError('{} is a mandatory column. Please make sure your input table has it.'.format(colname))

    # Generate components -- Should add to linetools.isgm.utils
    complist = []
    coords = SkyCoord(ra=table['RA'], dec=table['DEC'], unit='deg')
    for kk,row in enumerate(table):
        # Setup
        if 'flag_logN' in tkeys:
            Ntuple = (row['flag_logN'], row['logN'], row['sig_logN'])
        else:
            Ntuple = None
        # Units
        vmnx = [row['vmin'], row['vmax']]
        vmnx *= u.km/u.s if table['vmin'].unit is None else table['vmin'].unit
        Ej = row['Ej'] * (1/u.cm if table['Ej'].unit is None else table['Ej'].unit)
        #
        abscomp = AbsComponent(coords[kk], (row['Z'], row['ion']),
                               row['z_comp'], vmnx, Ej=Ej, Ntup=Ntuple)
        # Extra attrib
        for key in ['comment', 'name']:
            if key in tkeys:
                setattr(abscomp, key, row[key])
        # Other
        for key in ['sig_z', 'b','sig_b','specfile']:
            try:
                abscomp.attrib[key] = row[key]
            except KeyError:
                pass
            else:
                if table[key].unit is not None:
                    abscomp.attrib[key] *= table[key].unit
                else:
                    if key in ['b', 'sig_b']:
                        warnings.warn("No units for 'b'.  Will assume km/s")
                        abscomp.attrib[key] *= u.km/u.s
        complist.append(abscomp)
    return complist


def table_from_complist(complist):
    """
    Returns a astropy.Table from an input list of AbsComponents. It only
    fills in mandatory and special attributes (see notes below).
    Other information stored in dictionary AbsComp.attrib is ignored.

    Attributes with units are stored in the Table with units

    Parameters
    ----------
    complist : list of AbsComponents
        The initial list of AbsComponents to create the Table from.

    Returns
    -------
    table : Table
        Table from the information contained in each component.

    Notes
    -----
    Mandatory columns: 'RA', 'DEC', 'ion_name', 'z_comp', 'vmin', 'vmax'
    Special columns: 'name', 'comment', 'logN', 'sig_logN', 'flag_logN', 'b',
      'sig_b', 'specfile'
      'reliability' is included if provided
    See also complist_from_table()
    """
    key_order = ['RA', 'DEC', 'comp_name', 'z_comp', 'sig_z', 'Z', 'ion', 'Ej',
                 'vmin', 'vmax','ion_name', 'flag_N', 'logN', 'sig_logN',
                 'b','sig_b', 'vel', 'sig_vel','specfile']

    tab = Table()
    # Coordinates
    coords = SkyCoord([icomp.coord for icomp in complist])
    # Basics
    if hasattr(coords, 'ra'):
        tab['RA'] = coords.ra
        tab['DEC'] = coords.dec
    elif hasattr(coords, 'b'):
        tab['b_gal'] = coords.b
        tab['l_gal'] = coords.l
        # Adjust keys
        key_order.remove('RA')
        key_order.remove('DEC')
        key_order = ['b_gal', 'l_gal'] + key_order
    else:
        raise IOError("Not ready for this coords frame: {:s}".format(coords.frame.name))

    tab['comp_name'] = [comp.name for comp in complist]
    tab['vmin'] = u.Quantity([icomp.vlim[0] for icomp in complist])
    tab['vmax'] = u.Quantity([icomp.vlim[1] for icomp in complist])
    tab['Z'] = [icomp.Zion[0] for icomp in complist]
    tab['ion'] = [icomp.Zion[1] for icomp in complist]

    # . attributes (required ones)
    for attrib in ['zcomp', 'Ej', 'flag_N', 'logN', 'sig_logN',
                   'b','sig_b','vel','sig_vel']:
        values = [getattr(icomp,attrib) for icomp in complist]
        if isinstance(values[0], u.Quantity):
            values = u.Quantity(values)
        tab[attrib] = values
    # Rename
    tab.rename_column('zcomp', 'z_comp')

    # Ion names
    ion_names = []
    for comp in complist:
        if comp.Zion == (-1,-1):
            ion_names += ["Molecule"]
        else:
            ion_names += [ion_to_name(comp.Zion)]
    tab['ion_name'] = ion_names

    # attrib dict
    for attrib in ['sig_z', 'b', 'sig_b', 'vel', 'sig_vel', 'specfile']:
        try:
            values = [icomp.attrib[attrib] for icomp in complist]
        except KeyError:
            key_order.pop(key_order.index(attrib))
            pass
        else:
            # Quantity
            if isinstance(values[0], u.Quantity):
                values = u.Quantity(values)
            tab[attrib] = values

    # Special columns
    for key in ['comment', 'reliability']:
        if hasattr(complist[0], key):
            tab[key] = [getattr(comp, key) for comp in complist]
            key_order += [key]

    assert len(key_order) == len(tab.keys())
    tab = tab[key_order]

    # May need to add an HI component here as done in the method below

    return tab


def iontable_from_components(components, ztbl=None, NHI_obj=None):
    """Generate a Table from a list of components

    Method does *not* perform logic on redshifts or vlim.
    Includes rules for adding components of like ion
    Not ready for varying atomic mass (e.g. Deuterium)

    Parameters
    ----------
    components : list
      list of AbsComponent objects
    ztbl : float, optional
      Redshift for the table
    NHI_obj : object, optional (with NHI, sig_NHI, flag_NHI attributes)
      If provided, fill HI with NHI, sig_NHI, flag_NHI

    Returns
    -------
    iontbl : Table
    """
    warnings.warn("It is likely this method will be Deprecated", DeprecationWarning)
    warnings.warn("Use table_from_complist instead!!")
    from collections import OrderedDict
    # Checks
    assert chk_components(components,chk_A_none=True)

    # Set z from mean
    if ztbl is None:
        ztbl = np.mean([comp.zcomp for comp in components])

    # Construct the Table
    cols = OrderedDict()  # Keeps columns in order
    cols['Z']=int
    cols['ion']=int
    cols['A']=int
    cols['Ej']=float
    cols['z']=float
    cols['vmin']=float
    cols['vmax']=float
    cols['flag_N']=int
    cols['logN']=float
    if isinstance(components[0].sig_logN, float):
        cols['sig_logN'] = float
    elif components[0].sig_logN.size == 2:
        cols['sig_logN'] = np.ndarray
    else:
        raise IOError("Not prepared for this type of sig_logN")
    names = cols.keys()
    dtypes = [cols[key] for key in names]
    iontbl = Table(names=names,dtype=dtypes)
    iontbl['Ej'].unit=1./u.cm
    iontbl['vmin'].unit=u.km/u.s
    iontbl['vmax'].unit=u.km/u.s

    # Identify unique Zion, Ej (not ready for A)
    uZiE = np.array([comp.Zion[0]*1000000+comp.Zion[1]*10000+
                      comp.Ej.to('1/cm').value for comp in components])
    uniZi, auidx = np.unique(uZiE, return_index=True)

    # Loop
    for uidx in auidx:
        # Synthesize components with like Zion, Ej
        mtZiE = np.where(uZiE == uZiE[uidx])[0]
        comps = [components[ii] for ii in mtZiE]  # Need a list
        synth_comp = synthesize_components(comps, zcomp=ztbl)
        # Add a row to Table
        row = dict(Z=synth_comp.Zion[0],ion=synth_comp.Zion[1],
                   z=ztbl,
                   Ej=synth_comp.Ej,vmin=synth_comp.vlim[0],
                   vmax=synth_comp.vlim[1],logN=synth_comp.logN,
                   flag_N=synth_comp.flag_N,sig_logN=synth_comp.sig_logN)
        iontbl.add_row(row)

    # NHI
    if NHI_obj is not None:
        # Existing row in Table?
        mt = np.where((iontbl['Z'] == 1) & (iontbl['ion']==1))[0]
        if len(mt) == 1:
            iontbl[mt[0]]['logN'] = NHI_obj.NHI
            iontbl[mt[0]]['sig_logN'] = np.mean(NHI_obj.sig_NHI) # Allow for two values
            iontbl[mt[0]]['flag_N'] = NHI_obj.flag_NHI
        else:
            if len(components) > 0:
                vmin=synth_comp.vlim[0]
                vmax=synth_comp.vlim[1]
            else:
                vmin = -300*u.km/u.s
                vmax = 300*u.km/u.s
            #
            row = dict(Z=1,ion=1, z=ztbl,
                       Ej=0./u.cm,vmin=vmin, vmax=vmax, logN=NHI_obj.NHI,
                       flag_N=NHI_obj.flag_NHI,sig_logN=np.mean(NHI_obj.sig_NHI))
            iontbl.add_row(row)

    # Return
    return iontbl


def synthesize_components(components, zcomp=None, vbuff=0*u.km/u.s):
    """Synthesize a list of components into one

    Requires consistent RA/DEC, Zion, Ej, (A; future)
    Is agnostic about z+vlim
    Melds column densities
    Melds velocities with an optional buffer

    Note: Could make this a way to instantiate AbsComponent

    Parameters
    ----------
    components : list
      list of AbsComponent objects
    zcomp : float, optional
      Input z to reference the synthesized component
      If not input, the mean of the input components is used
    vbuff : Quantity, optional
      Buffer for synthesizing velocities.  Deals with round off, c, etc.
    """
    # Checks
    assert chk_components(components, chk_A_none=True, chk_match=True)


    # Meld column densities
    obj = dict(flag_N=components[0].flag_N, logN=components[0].logN,
               sig_logN=components[0].sig_logN)
    for comp in components[1:]:
        if comp.flag_N != 0:
            obj['flag_N'], obj['logN'], obj['sig_logN'] = ltaa.sum_logN(obj, comp)

    # zcomp
    if zcomp is None:
        zcomp = np.mean([comp.zcomp for comp in components])

    # Set vlim by min/max  [Using non-relativistic + buffer]
    vmin = u.Quantity([(comp.zcomp-zcomp)/(1+zcomp)*const.c.to('km/s')+comp.vlim[0] for comp in components])
    vmax = u.Quantity([(comp.zcomp-zcomp)/(1+zcomp)*const.c.to('km/s')+comp.vlim[1] for comp in components])
    vlim = u.Quantity([np.min(vmin)-vbuff, np.max(vmax)+vbuff])

    # Init final component
    synth_comp = AbsComponent(components[0].coord, components[0].Zion, zcomp,
                              vlim, Ej=components[0].Ej, stars=components[0].stars,
                              Ntup=(obj['flag_N'], obj['logN'], obj['sig_logN']))

    # Return
    return synth_comp


def get_components_at_z(complist, z, dvlims):
    """In a given list of AbsComponents, it finds
    the ones that are within dvlims from a given redshift
    and returns a list of those.

    Parameters
    ----------
    complist : list
        List of AbsComponents
    z : float
        Redshift to search for components
    dvlims : Quantity array
        Rest-frame velocity limits around z
        to look for components

    Returns
    -------
    components_at_z : list
        List of AbsComponents in complist within dvlims from z
    """
    # check input
    if not isinstance(complist[0], AbsComponent):
        raise IOError('complist must be a list of AbsComponents.')
    if len(dvlims) != 2:
        raise IOError('dvlims must be a Quantity array of velocity limits (vmin, vmax).')
    else:
        try:
            dvlims_kms = dvlims.to('km/s')
        except u.UnitConversionError:
            raise IOError('dvlims must have velocity units.')

    good_complist = []
    for comp in complist:
        dv_comp = ltu.dv_from_z(comp.zcomp, z)
        if (dv_comp >= dvlims[0]) and (dv_comp <= dvlims[1]):
            good_complist += [comp]
    return good_complist


def get_wvobs_chunks(comp):
    """For a given component, it gets a list of tuples with the
    min/max observed wavelengths for each absorption line in the
    component. An error is raised if an absorption line within the
    component does not have its limits defined.

    Parameters
    ----------
    comp : AbsComponent
        The input AbsComponent object

    Returns
    -------
    wvobs_chunks : list of Quantity arrays
        A list with the wvmin, wvmax values for each absorption
        line within the component.
    """

    if not isinstance(comp, AbsComponent):
        raise ValueError('`comp` must be AbsComponent object.')

    wvobs_chunks = []
    for absline in comp._abslines:
        # Check whether the absline has already defined 'wvlim'
        if absline.limits.is_set():
            wvlim_aux = absline.limits.wvlim
            wvobs_chunks += [wvlim_aux]
        else:
            raise ValueError('{} must have its limits defined.'.format(absline))
    return wvobs_chunks


def coincident_components(comp1, comp2, tol=0.2*u.arcsec):
    """Whether two components overlap in wavelength (observed)
    space and (ra,dec) sky position. This is useful to identify
    components that may need to be fit together in a given spectrum.

    Parameters
    ----------
    comp1 : AbsComponent
        A given AbsComponent object
    comp2 : AbsComponent
        A given AbsComponent object
    tol : Quantity, optional
        Tolerance for checking whether the two components are
        in the same sky region. Default is 0.2*u.arcsec

    Returns
    -------
    answer : bool
        True if there is overlapping wavelength range and
        radec coordinates, otherwise False.
    """

    if not isinstance(comp1, AbsComponent):
        raise ValueError('comp1 must be AbsComponent object.')
    if not isinstance(comp2, AbsComponent):
        raise ValueError('comp1 must be AbsComponent object.')

    # Check whether they are in the same sky region
    if comp1.coord.separation(comp2.coord) > tol:
        return False

    # loop over abslines
    for line1 in comp1._abslines:
        for line2 in comp2._abslines:
            overlap = line1.coincident_line(line2)
            if overlap is True:
                return True
    return False


def group_coincident_components(comp_list, output_type='list'):
    """For a given input list of components, this function
    groups together components that are coincident to each other
    (including by transitivity), and returns them as a list (default)
    or dictionary of component lists.

    Parameters
    ----------
    comp_list : list of AbsComponent
        Input list of components to group
    output_type : str, optional
        Type of the output, choose either
        'list' for list or 'dict' for dictionary.

    Returns
    -------
    output : list (or dictionary) of lists of AbsComponent
        The grouped components as individual lists
        in the output.
    """
    if output_type not in ['list', 'dict', 'dictionary']:
        raise ValueError("`output_type` must be either 'list' or 'dict'.")

    ### We first want to identify and group all blended lines
    ### Sort them by observed wavelength to do this
    lst=[]
    compnos=[]
    for ii,comp in enumerate(comp_list):
        lst.extend(comp._abslines)
        compnos.extend([ii]*len(comp._abslines))
    lst=np.array(lst)
    compnos=np.array(compnos)
    wv1s=np.array([line.limits.wvlim[0].value for line in lst])
    sortidxs=np.argsort(wv1s)
    sort_lst=lst[sortidxs]
    sort_compnos=compnos[sortidxs] # This will store indices of the lines' parent comps.

    ### Identify the blended 'absline' objects
    blends = []
    ### 'blends' is a list of lists
    ### each sublist will contain the indices of consecutive blended abslines in wobs space
    for i in range(len(sortidxs)-1):
        if i == 0:
            thisblend = [i]
        if sort_lst[i].coincident_line(sort_lst[i+1]):
            thisblend.append(i+1)
        else:
            blends.append(thisblend)
            thisblend = [i+1]
        if i == (len(sortidxs) - 2):
            blends.append(thisblend)

    ### Associate the lines to their parent components
    blendnos=[]
    for blist in blends:
        blendnos.append(sort_compnos[blist])

    ### Main algorithm to group together all components with blended lines
    ## Each sublist in 'blends' will be checked to see if the parent components
    ## of the abslines have abslines in other blend sublists. These sublists
    ## that share parent components will be grouped together.   Then,
    compfound=[]  # will hold components that have been grouped
    grblends=[]  # will hold indices of blends that have been grouped
    newgroups=[]  # will hold the newly grouped abslines
    for i,bn in enumerate(blendnos):  # bn is list of indices of parent comps of blended abslines
        if i in grblends: continue  # move on if blend has already been grouped
        grblends.append(i)  # so that we don't try to group this sublist twice
        newgroups.append(bn.tolist()) # start group with parent components of this blend
        newtotry=bn  # a list of components, which have assoc. abslines, which may be in other blends
        ## check all of these components' abslines for inclusion in blend groups
        while (len(newtotry)>0):  # stop when all potential lines in grouped components have been checked
            newnewtotry=[]
            for no in newtotry:
                if no not in compfound:  # Move on if a component's lines have been checked
                    compfound.append(no)  # So that we don't group this component twice
                    blgrs=_whichgroupscontainmember(blendnos,no)  # Find which sublists contain abslines for this comp
                    for bg in blgrs:  # Go through these sublists
                        if bg not in grblends:  # If this sublist hasn't been checked, save it
                            grblends.append(bg)
                            newgroups[-1].extend(blendnos[bg].tolist())  # Add comp numbers to this group
                            # Now need to check parent comps of abslines in newly added blend sublist
                            newnewtotry.extend(np.unique(blendnos[bg]).tolist())
            newtotry=newnewtotry  # Get ready for next round of checking
            newgroups[-1]=np.unique(np.array(newgroups[-1])).tolist()  # Clean up duplicated comp indices

    out = newgroups

    # Now lets produce the final output from it
    output_list = []
    output_dict = {}
    for ii in range(len(out)):
        aux_list = []
        for jj in out[ii]:
            aux_list += [comp_list[jj]]
        output_dict['{}'.format(ii)] = aux_list
        output_list += [aux_list]

    # choose between dict of list
    if output_type == 'list':
        return output_list
    elif output_type in ['dict', 'dictionary']:
        return output_dict

def _whichgroupscontainmember(groups,member):
    ''' Once blends have been identified, find which blend a given line index belongs to.
    Parameters
    ----------
    groups : list of lists
        Input list of groups within which you want to find a member
        Here, this is used for groups of blended line indices
    member : object with same type as those in 'groups'
        Value for which to search within groups.
        Here, this corresponds to a specific line index

    Returns
    -------
    matches : list
        Indices of groups within with 'member' was found
    '''
    matches=[]
    for i,gr in enumerate(groups):
        if member in gr:
            matches.append(i)
    return matches

def group_coincident_components_old(comp_list, output_type='list'):
    """For a given input list of components, this function
    groups together components that are coincident to each other
    (including by transitivity), and returns them as a list (default)
    or dictionary of component lists.

    Parameters
    ----------
    comp_list : list of AbsComponent
        Input list of components to group
    output_type : str, optional
        Type of the output, choose either
        'list' for list or 'dict' for dictionary.

    Returns
    -------
    output : list (or dictionary) of lists of AbsComponent
        The grouped components as individual lists
        in the output.
    """
    if output_type not in ['list', 'dict', 'dictionary']:
        raise ValueError("`output_type` must be either 'list' or 'dict'.")

    # the first extreme case is that all components are independent
    # of each other, in which case we have the following output shape
    out = [[] for kk in range(len(comp_list))]

    for ii in range(len(comp_list)):
        comp_ii = comp_list[ii]
        # only append if ii does not belong to a previous round
        switch = 0
        for kk in range(len(out[:ii])):
            if ii in out[kk]:
                switch = 1
                break
        if switch == 1:
            pass
        else:
            out[ii].append(ii)

        for jj in range(ii+1, len(comp_list)):
            # print(ii,jj)
            comp_jj = comp_list[jj]
            if coincident_components(comp_ii, comp_jj):  # There is overlap between comp_ii and comp_jj
                # check in the previous ones where does jj belongs to
                switch = 0
                for kk in range(len(out[:ii])):
                    if ii in out[kk]:  # this means ii already belongs to out[kk]
                        # so jj should also go there...(if not there already)
                        if jj in out[kk]:
                            pass
                        else:
                            out[kk].append(jj)
                        switch = 1  # for not appending jj again
                        break
                if switch == 1:
                    pass  # this jj was appended already
                else:
                    # but check is not already there...
                    if jj in out[ii]:
                        pass
                    else:
                        out[ii].append(jj)
                # print(out)
            else:
                pass

    # Now we have out as a list of lists with indices, or empty lists
    # let's get rid of the empty lists
    out = [x for x in out if x != []]

    # Now lets produce the final output from it
    output_list = []
    output_dict = {}
    for ii in range(len(out)):
        aux_list = []
        for jj in out[ii]:
            aux_list += [comp_list[jj]]
        output_dict['{}'.format(ii)] = aux_list
        output_list += [aux_list]

    # choose between dict of list
    if output_type == 'list':
        return output_list
    elif output_type in ['dict', 'dictionary']:
        return output_dict


def unique_components(comps1, comps2, tol=5*u.arcsec):
    """ Identify which AbsComponent members of the comps1 list
    are *not* within the comps2 list, to given tolerances.
    Note, AbsComponent objects in the comps1 list are not examined
    against each other for uniqueness.

    Unique if any apply (test is done in this order)
      1) coord.separation > tol
      2) Z,ion,Ej set is unique
      3) redshift limits do not overlap

    Parameters
    ----------
    comps1 : list of AbsComponent objects
    comps2 : list of AbsComponent objects

    Returns
    -------
    unique : bool array
      True = members of comps1 that are not currently in comps2

    """
    c_mks = const.c.to('km/s').value
    unique = np.array([True]*len(comps1))
    # Coordinates
    ras = [icomp.coord.ra.value for icomp in comps1]
    decs = [icomp.coord.dec.value for icomp in comps1]
    coords1 = SkyCoord(ra=ras, dec=decs, unit='deg')
    ras = [icomp.coord.ra.value for icomp in comps2]
    decs = [icomp.coord.dec.value for icomp in comps2]
    coords2 = SkyCoord(ra=ras, dec=decs, unit='deg')
    # Compare
    idx, d2d, d3d = match_coordinates_sky(coords1, coords2, nthneighbor=1)
    close_enough = d2d < tol
    if np.sum(close_enough) == 0:
        return unique
    # Next step (Z, ion, Ej)
    ZiE1 = np.array([(icomp.Zion[0], icomp.Zion[1], icomp.Ej.value) for icomp in comps1])
    ZiE2 = np.array([(icomp.Zion[0], icomp.Zion[1], icomp.Ej.value) for icomp in comps2])
    indices = np.where(close_enough)[0]
    for idx in indices:  # comp1 indices
        # Match on coords
        coord_mt = np.where(coords1[idx].separation(coords2) < tol)[0]
        # Match on ZiE
        mtZiE = np.where((ZiE2[coord_mt] == ZiE1[idx]).all(axis=1))[0]
        if len(mtZiE) > 0: # Lastly redshift
            zlim_comp1 = comps1[idx].limits.zlim
            for idx2 in coord_mt[mtZiE]:
                zlim_comp2 = comps2[idx2].limits.zlim
                # Redshift overlap?
                if (zlim_comp1[0] < zlim_comp2[1]) & (zlim_comp1[1] > zlim_comp2[0]):
                    unique[idx] = False
    # Return
    return unique


