"""
Module related to analysis of absorption lines
  -- Methods
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from astropy import units as u
from astropy import constants as const


# Perform AODM on the line
def aodm(spec,idata):
    """  AODM calculation 
    See Savage & Sembach 1991, ApJ, 379, 245

    Parameters
    ----------
    spec -- Tuple of (wave, fx, sig)
    idata -- Tuple of (wrest, fval)

    Returns:
    N, sigN : float, float
      Column and error in linear space (cm^-2)
    flg_sat: bool
      Set to True if saturated pixels exist
    """
    flg_sat = False

    # Cut spectrum
    velo,fx,sig = spec
    wrest, fval = idata

    # dv
    delv = velo - np.roll(velo,1)
    delv[0] = delv[1]

    # Atomic data
    cst = (const.m_e.cgs*const.c.cgs / (np.pi * 
        (const.e.esu**2).cgs)).to(u.AA*u.s/(u.km*u.cm**2))
    cst = cst/(fval*wrest) #/ (u.km/u.s) / u.cm * (u.AA/u.cm)

    # Mask
    mask = (fx == fx) # True = good
    nndt = u.Quantity(np.zeros(len(fx)), unit='s/(km cm cm)')

    # Saturated?
    satp = np.where( (fx <= sig/5.) | (fx < 0.05) )[0]
    if len(satp) > 0:
        mask[satp] = False
        lim = np.where(sig[satp] > 0.)[0]
        if len(lim) > 0:
            sub = np.maximum(0.05, sig[satp[lim]]/5.)
            nndt[satp[lim]] = np.log(1./sub)*cst
            flg_sat = True
    # AODM
    nndt[mask] = np.log(1./fx[mask])*cst

    # Sum it
    ntot = np.sum( nndt*delv )
    tvar = np.sum( (delv*cst*sig/fx)**2 )


    # Return
    return ntot, np.sqrt(tvar), flg_sat
