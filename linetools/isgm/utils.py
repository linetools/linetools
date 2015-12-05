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
from collections import OrderedDict

from astropy import constants as const
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.io import ascii

from linetools.spectralline import AbsLine
from linetools.analysis import absline as ltaa
from linetools.isgm.abscomponent import AbsComponent

def chk_components(components, chk_match=False, chk_A_none=False, toler=0.2*u.arcsec):
    """ Performs checks on a list of components

    Parameters
    ----------
    components : list
      list of AbsComponent objects
    chk_match : bool, optional
      if True, require that the components match in RA/DEC, Zion, Ej, A but not velocity
    chk_A_none : bool, optional
      if True, require that A *not* be set
    toler : Angle, optional
      Tolerance on matching coordinates
    """
    tests = True
    # List
    if not isinstance(components,list):
        tests = False
        raise IOError('Need a list of AbsComponent objects')
    # Object
    if not all(isinstance(x,AbsComponent) for x in components):
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
            match = match & bool(comp0.coord.separation(comp.coord) < toler)
            # Zion
            match = match & (comp0.Zion == comp.Zion)
            # Ej
            match = match & np.allclose(comp0.Ej.to('1/cm').value,comp.Ej.to('1/cm').value)
            # A
            match = match & (comp0.A == comp.A)
        tests = tests & match
    # Return
    return tests

def build_components_from_abslines(iabslines, clmdict=None, coord=None):
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
        abslines = []
        for wrest in clmdict['lines'].keys():
            clmdict['lines'][wrest].attrib['coord'] = coord
            abslines.append(clmdict['lines'][wrest])
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
        component = AbsComponent.from_abslines(lines)
        # Reset vmin, vmax
        vmin,vmax = 9999., -9999.
        for iline in lines:
            vmin = min(vmin, iline.analy['vlim'][0].value)
            vmax = max(vmax, iline.analy['vlim'][1].value)
        component.vlim = [vmin,vmax]*u.km/u.s
        # Append
        components.append(component)

    return components



def iontable_from_components(components,ztbl=None):
    """Generate a QTable from a list of components

    Method does *not* perform logic on redshifts or vlim.
    Includes rules for adding components of like ion
    Not ready for varying atomic mass (e.g. Deuterium)

    Parameters
    ----------
    components : list
      list of AbsComponent objects
    ztbl : float, optional
      Redshift for the table
    """
    from collections import OrderedDict
    # Checks
    assert chk_components(components,chk_A_none=True)

    # Set z from mean
    if ztbl is None:
        ztbl = np.mean([comp.zcomp for comp in components])

    # Construct the QTable
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
    iontbl = QTable(names=names,dtype=dtypes)
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
        comps = [components[ii] for ii in mtZiE] # Need a list
        synth_comp = synthesize_components(comps,zcomp=ztbl)
        # Add a row to QTable
        row = dict(Z=synth_comp.Zion[0],ion=synth_comp.Zion[1],
                   z=ztbl,
                   Ej=synth_comp.Ej,vmin=synth_comp.vlim[0],
                   vmax=synth_comp.vlim[1],logN=synth_comp.logN,
                   flag_N=synth_comp.flag_N,sig_logN=synth_comp.sig_logN)
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
