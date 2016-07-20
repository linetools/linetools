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
from astropy.units import Quantity

from linetools.analysis import absline as ltaa
from linetools.isgm.abscomponent import AbsComponent
from linetools.utils import give_dz
from linetools.spectralline import init_analy

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
      Tolerance on matching SkyCoordinates. Default is 0.2*u.arcsec
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
            vmin = min(vmin, iline.analy['vlim'][0].value)
            vmax = max(vmax, iline.analy['vlim'][1].value)
        component.vlim = [vmin,vmax]*u.km/u.s
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
                line = AbsLine.from_dict(idict['lines'][key], coord=coord)
            else:
                raise IOError("Need those lines")
            if coord is not None:
                line.attrib['coord'] = coord
            lines.append(line)
        components = build_components_from_abslines(lines, **kwargs)
    else:
        warnings.warn("No components in this dict")
    # Return
    return components


def build_systems_from_components(comps, **kwargs):
    """ Build a list of AbsSystems from a list of AbsComponents
    Current default implementation allows for overlapping components, i.e.
      only_overalp=True in add_component

    Parameters
    ----------
    comps : list

    Returns
    -------
    abs_systems : list

    """
    from linetools.isgm.abssystem import GenericAbsSystem
    if 'overlap_only' not in kwargs.keys():
        kwargs['overlap_only'] = True
    # Add
    abs_systems = []
    cpy_comps = [comp.copy() for comp in comps]
    # Loop until all components assigned
    while len(cpy_comps) > 0:
        # Use the first one
        comp = cpy_comps.pop(0)
        abssys = GenericAbsSystem.from_components([comp])
        abs_systems.append(abssys)
        # Try the rest
        comps_left = []
        for icomp in cpy_comps:
            if abssys.add_component(icomp, **kwargs):
                pass
            else:
                comps_left.append(icomp)
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
    cols['sig_logN']=float
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
        # Add a row to QTable
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
            iontbl[mt[0]]['sig_logN'] = NHI_obj.sig_NHI
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
                       flag_N=NHI_obj.flag_NHI,sig_logN=NHI_obj.sig_NHI)
            iontbl.add_row(row)

    # Add zlim to metadata
    meta = OrderedDict()
    meta['zcomp'] = ztbl

    # Return
    return iontbl


def synthesize_components(components, zcomp=None, vbuff=0*u.km/u.s):
    """Synthesize a list of components into one

    Requires consistent RA/DEC, Zion, Ej, (A; future)
    Is agnostic about z+vlim
    Melds column densities
    Melds velocities with a small buffer (10 km/s)

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

    # Init final component
    synth_comp = AbsComponent.from_component(components[0], Ntup=(components[0].flag_N, components[0].logN, components[0].sig_logN))

    # Meld column densities
    for comp in components[1:]:
        if comp.flag_N != 0:
            synth_comp.flag_N, synth_comp.logN, synth_comp.sig_logN = ltaa.sum_logN(synth_comp, comp)

    # Meld z, vlim
    # zcomp
    if zcomp is None:
        zcomp = np.mean([comp.zcomp for comp in components])
    synth_comp.zcomp = zcomp
    # Set vlim by min/max  [Using non-relativistic + buffer]
    vmin = u.Quantity([(comp.zcomp-zcomp)/(1+zcomp)*const.c.to('km/s')+comp.vlim[0] for comp in components])
    vmax = u.Quantity([(comp.zcomp-zcomp)/(1+zcomp)*const.c.to('km/s')+comp.vlim[1] for comp in components])
    synth_comp.vlim = u.Quantity([np.min(vmin)-vbuff, np.max(vmax)+vbuff])

    # Return
    return synth_comp


def overlapping_chunks(chunk1, chunk2):
    """True if there is overlap between chunks
    `chunk1` and `chunk2`. Otherwise False. Chunks are
    assumed to represent continuous coverage, so the only
    information that matters are the minimum and maximum
    values of a given chunk.

    Parameters
    ----------
    chunk1 : tuple, list, 1-d np.array, Quantity
        A given chunk, assumed to represent a contiguous region
        So only its minimum and maximum values matter
        Still, chunk must be sorted.
    chunk1 : tuple, list, 1-d np.array
        Ditto.

    Returns
    -------
    Boolean, True if there is overlap, False otherwise.

    """
    # Check units in case chunks are Quantity
    if isinstance(chunk1, Quantity):
        unit1 = chunk1.unit
        chunk1 = np.array(chunk1.value)
        if not isinstance(chunk2, Quantity):
            raise ValueError('chunk2 must be Quantity because chunk1 is!')
        try:
            chunk2 = chunk2.to(unit1)
            chunk2 = np.array(chunk2.value)  # has the same units as chunk1
        except u.core.UnitConversionError:
            raise ValueError('If chunks are given as Quantity they must have convertible units!')
    else:
        chunk1 = np.array(chunk1)  # not a quantity

    # this may be redundant but cleaner code
    if isinstance(chunk2, Quantity):
        unit2 = chunk2.unit
        chunk2 = np.array(chunk2.value)
        if not isinstance(chunk1, Quantity):
            raise ValueError('chunk1 must be Quantity because chunk2 is!')
    else:
        chunk2 = np.array(chunk2)  # not a quantity
    # here we have chunks as values (i.e no units)

    # make sure the chunks are sorted
    cond1 = np.sort(chunk1) != chunk1
    cond2 = np.sort(chunk2) != chunk2
    if (np.sum(cond1) > 0) or (np.sum(cond2) > 0):
        raise ValueError('chunks must be sorted!')

    # figure out the lowest of the chunks
    # and sort them such that chunk1 is by definition
    # the one with the lowest limit
    if np.min(chunk1) <= np.min(chunk2):
        pass
    else: # invert them if necessary
        aux = chunk1
        chunk1 = chunk2
        chunk2 = aux

    # now the chunk1 is the one with the lowest value
    # then, the only way they do not overlap is:
    if np.min(chunk2) > np.max(chunk1):
        return False
    else:  # overlap exists
        return True

def overlapping_components(comp1, comp2, tol=0.2*u.arcsec):
    """Whether two components overlap in wavelength
    space and (ra,dec) sky position. This is useful to identify
    components that may need to be fit together.

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
    Boolean : True if there is overlapping wavelength range
    and radec coordinates, otherwise False.
    """

    if not isinstance(comp1, AbsComponent):
        raise ValueError('comp1 must be AbsComponent object.')
    if not isinstance(comp2, AbsComponent):
        raise ValueError('comp1 must be AbsComponent object.')

    # Check whether they are in the same sky region
    if comp1.coord.separation(comp2.coord) > tol:
        return False

    # Define wavelength chunks, appending them from each absline
    # These will correpond to (wvmin, wvmax) tuples for each absline
    # in the respective component
    wobs_chunks_1 = []
    wobs_chunks_2 = []
    for absline in comp1._abslines:
        cond = absline.analy['vlim'] != init_analy['vlim']  # default value?
        if np.sum(cond) > 0:  # i.e. absline vlim not equal than the default
            dzlim = give_dz(absline.analy['vlim'], comp1.zcomp)
        else:  # use the vlim from component otherwise
            dzlim = give_dz(comp1.vlim, comp1.zcomp)
        wrest_aux = absline.wrest.to('AA').value  # in AA
        wrest_obs_aux = wrest_aux + (1 + comp1.zcomp + dzlim)
        wobs_chunks_1 += [tuple(wrest_obs_aux)]

    for absline in comp2._abslines:
        cond = absline.analy['vlim'] != init_analy['vlim']  # default value
        if np.sum(cond) > 0:  # i.e. absline vlim not equal than the default
            dzlim = give_dz(absline.analy['vlim'], comp2.zcomp)
        else:  # use the vlim from component otherwise
            dzlim = give_dz(comp2.vlim, comp2.zcomp)
        wrest_aux = absline.wrest.to('AA').value  # in AA
        wrest_obs_aux = wrest_aux + (1 + comp2.zcomp + dzlim)
        wobs_chunks_2 += [tuple(wrest_obs_aux)]

    # now we have two lists of wobs chunks:
    # so lets do the checking until at least 1
    # overlap is found
    for ii in range(len(wobs_chunks_1)):
        for jj in range(len(wobs_chunks_2)):
            overlap = overlapping_chunks(
                wobs_chunks_1[ii], wobs_chunks_2[jj])
            if overlap is True:
                return True
    return False
