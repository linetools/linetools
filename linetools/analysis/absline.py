"""
Module related to analysis of absorption lines
  -- Methods
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os

from astropy import units as u
from astropy import constants as const

# Atomic constant
atom_cst = (const.m_e.cgs*const.c.cgs / (np.pi * 
    (const.e.esu**2).cgs)).to(u.AA*u.s/(u.km*u.cm**2))

# Perform AODM on the line
def aodm(spec,idata):
    """  AODM calculation on an absorption line
    See Savage & Sembach 1991, ApJ, 379, 245

    Parameters
    ----------
    spec : tuple
      (vel, fx, sig)
    idata : tuple
      (wrest, fval)

    Returns
    -------
    N, sigN : float, float
      Column and error in linear space (cm^-2)
    flg_sat : bool
      Set to True if saturated pixels exist

    ToDo: 
      -- Generate a function for nndt alone
    """
    flg_sat = False

    # Cut spectrum
    velo,fx,sig = spec
    wrest, fval = idata

    # dv
    delv = velo - np.roll(velo,1)
    delv[0] = delv[1]

    # Atomic data
    cst = atom_cst/(fval*wrest) #/ (u.km/u.s) / u.cm * (u.AA/u.cm)

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

def log_clm(obj):
    """Return logN and sig_logN given linear N, sigN
    Parameters
    ---------
    obj : object
      An object with tags appropriate for the analysis
      Assumes 'logN' for column and 'sig_logN' for error for now
      Or .N and .sigN

    Returns
    -------
    logN : float
      log10 N
    sig_logN :float
      Error in log10 N
    """
    # Grab
    try:
        iN = obj['N']
    except TypeError:
        try:
            iN = obj.N
        except:
            raise IOError("Bad input to log_clm")
        else:
            isigN = obj.sigN
    else:
        isigN = obj['sigN']

    # Strip units
    try:
        N = iN.value
    except AttributeError:
        N = iN
        sigN = isigN
    else:
        sigN = isigN.value

    # Operate
    if N <= 0.:
        logN = 0.
    else:
        logN = np.log10(N)
    lgvar = ((1.0 / (np.log(10.0)*N))**2)*sigN**2
    sig_logN = np.sqrt(lgvar)
    # Return
    return logN, sig_logN

def sum_logN(obj1,obj2):
    """Add log columns and return value and errors
    Parameters
    ----------
    obj 1: object
      An object with tags appropriate for the analysis
      Assumes 'logN' for column and 'sig_logN' for error for now
    obj2 : object
      Another object with tags appropriate for the analysis

    Returns
    -------
    logN, siglogN
    """
    # Calculate
    logN = np.log10(np.sum(10.**np.array([obj1['logN'],obj2['logN']])))
    siglogN = np.sqrt(
        np.sum([(obj1['sig_logN']*(10.**obj1['logN']))**2,
        (obj2['sig_logN']*(10.**obj2['logN']))**2]))/(10.**logN)
    # Return
    return logN, siglogN