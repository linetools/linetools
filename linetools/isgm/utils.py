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


def get_wvobs_chunks(comp):
    """For a given component, it gets a list of tuples with the
    min/max observed wavelengths for each absorption line in the
    component.

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
        cond = absline.analy['wvlim'] != init_analy['wvlim']
        if np.sum(cond) > 0: # not default, use these values then
            wvlim_aux = absline.analy['wvlim']
            wvobs_chunks += [wvlim_aux]
        else:
            # Check whether the absline has already defined 'vlim'
            cond = absline.analy['vlim'] != init_analy['vlim']
            if np.sum(cond) > 0:  # i.e. absline vlim not equal than the default, so we use these
                # define the best redshift
                if absline.attrib['z'] != 0:
                    zline = absline.attrib['z']
                else:
                    zline = comp.zcomp
                dzlim = give_dz(absline.analy['vlim'], zline)
            else:  # use the vlim from component otherwise
                zline = comp.zcomp
                dzlim = give_dz(comp.vlim, zline)
            wrest_aux = absline.wrest
            wvlim_aux = wrest_aux * (1 + zline + dzlim)
            wvobs_chunks += [wvlim_aux]
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


def group_coincident_compoments(comp_list):
    """For a given input list of components, this function
    groups together components that are coincident to each other
    (including by transitivity), and returns them as a list of
    component lists.

    Parameters
    ----------
    comp_list : list of AbsComponent
        Input list of components to group

    Returns
    -------
    output : list of lists of AbsComponent
        The grouped components as individual lists
        in the output list.
    """
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
    # return output_dict
    return output_list


