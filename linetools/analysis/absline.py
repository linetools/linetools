""" Utlities for the analysis of absorption lines
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import warnings

from astropy import units as u
from astropy import constants as const
from astropy.io import ascii

# Atomic constant
atom_cst = (const.m_e.cgs*const.c.cgs / (np.pi * 
    (const.e.esu**2).cgs)).to(u.AA*u.s/(u.km*u.cm**2))

# Perform AODM on the line
def aodm(spec, idata):
    """ AODM calculation on an absorption line

    See Savage & Sembach 1991, ApJ, 379, 245

    Parameters
    ----------
    spec : tuple
      (vel, fx, sig)
    idata : tuple
      (wrest, fval)

    Returns
    -------
    N, sig_N : float, float
      Column and error in linear space (cm^-2)
    flg_sat : bool
      Set to True if saturated pixels exist
    """
    # ToDo: 
    #   -- Generate a function for nndt alone


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
    """Return logN and sig_logN given linear N, sig_N

    Also fills the attributes

    Parameters
    ----------
    obj : object
      An object with tags appropriate for the analysis
      Assumes 'logN' for column and 'sig_logN' for error for now
      Or .N and .sig_N

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
            isig_N = obj.sig_N
    else:
        isig_N = obj['sig_N']

    # Strip units
    try:
        N = iN.value
    except AttributeError:
        N = iN
        sig_N = isig_N
    else:
        sig_N = isig_N.value

    # Operate
    if N <= 0.:
        logN = 0.
    else:
        logN = np.log10(N)
    lgvar = ((1.0 / (np.log(10.0)*N))**2)*sig_N**2
    sig_logN = np.sqrt(lgvar)

    # Fill
    try:
        obj['logN'] = logN
    except:
        obj.logN = logN
        obj.sig_logN = sig_logN
    else:
        obj['sig_logN'] = sig_logN
    # Return
    return logN, sig_logN


def linear_clm(obj):
    """Return N and sig_N given logN, sig_logN in an object

    Also fills the attributes

    Parameters
    ----------
    obj : object
      An object with tags appropriate for the analysis
      Assumes 'logN' for column and 'sig_logN' for error
      Or .logN and .sig_logN

    Returns
    -------
    N : Quantity
      column in cm^-2
    sig_N : Quantity
      error in column in cm^-2
    """
    # Grab
    try:
        logN = obj['logN']
    except TypeError:
        try:
            logN = obj.logN
        except:
            raise IOError("Bad input to linear_clm")
        else:
            sig_logN = obj.sig_logN
    else:
        sig_logN = obj['sig_logN']


    # Operate
    N = 10**logN / u.cm**2
    sig_N = sig_logN * np.log(10.) * N

    # Fill
    try:
        obj['sig_N'] = sig_N
    except:
        obj.N = N
        obj.sig_N = sig_N
    else:
        obj['N'] = N
    # Return
    return N, sig_N


def photo_cross(Z, ion, E, datfil=None, silent=False):
    """ Estimate photo-ionization cross-section using Fit parameters

    from Verner et al. 1996, ApJ, 465, 487
    JXP on 04 Nov 2014

    Parameters
    ----------
    Z : int
      Atomic number
    ion : int
      Ionization state (1=Neutral)
    E : Quantity
      Energy to calculate at [eV]

    Returns
    -------
    sigma : Quantity
      Cross-section (cm^2)
    """
    import imp
    lt_path = imp.find_module('linetools')[1]
    # Read data
    if datfil is None:
        datfil = lt_path+'/data/atomic/verner96_photoion_table1.dat'
    dat = ascii.read(datfil)

    # Deal with Units
    if not isinstance(E, u.quantity.Quantity):
        if silent is False:
            warnings.warn('photo_cross: Assuming eV for input energy')
        E = E * u.eV

    # Match
    mt = np.where((Z == dat['Z']) & (ion == dat['N']))[0]
    nmt = len(mt)
    if nmt == 0:
        raise ValueError('photo_cross: {:d},{:d} pair not in our table'.format(Z, ion))
    idx = mt[0]
    #
    x = E/(dat['E0'][idx]*u.eV) - dat['y0'][idx]
    y = np.sqrt(x**2 + dat['y1'][idx]**2)

    F = (((x-1.)**2 + dat['yw'][idx]**2) * y**(0.5*dat['P'][idx] - 5.5) *
            (1 + np.sqrt(y/dat['ya'][idx]) )**(-1.*dat['P'][idx]))

    sigma = dat['s0'][idx] * F * 1e-18 * u.cm**2

    # Energy threshold
    low = np.where(E < dat['Eth'][idx]*u.eV)[0]
    if len(low) > 0:
        sigma[low] = 0.
    # Return
    return sigma


def sum_logN(obj1, obj2):
    """Add log columns and return value and errors, with logic

    Adds two column density objects, taking into account the flags

    Parameters
    ----------
    obj1 : object
      An object with keys or attributes appropriate for the analysis
      Assumes 'logN' for column and 'sig_logN' for error for now
    obj2 : object
      Another object with keys appropriate for the analysis

    Returns
    -------
    logN, siglogN
    """
    # Check
    if not (obj1['flag_N'] in [1, 2, 3]):
        raise ValueError("flag_N in obj1 must be 1,2,3")
    if not (obj2['flag_N'] in [1, 2, 3]):
        raise ValueError("flag_N in obj2 must be 1,2,3")
    # Unpack Obj1
    flag_N, logN1, sig_logN1 = [obj1[key] for key in ['flag_N','logN','sig_logN']]
    # Sum
    logN = np.log10(np.sum(10.**np.array([obj1['logN'],obj2['logN']])))
    sig_logN = np.sqrt( np.sum([(obj1['sig_logN']*(10.**obj1['logN']))**2,
                                (obj2['sig_logN']*(10.**obj2['logN']))**2]))/(10.**logN)
    if flag_N in [1,2]: # Detection or saturated
        if obj2['flag_N'] == 2:
            flag_N = 2
        elif obj2['flag_N'] == 3:
            # No change to logN, only sig_logN
            logN = logN1
    elif flag_N == 3:
        if obj2['flag_N'] in [1,2]:
            logN = obj2['logN']
            flag_N = obj2['flag_N']
            sig_logN = obj2['sig_logN']
        else:
            pass # Take sums
    # Return
    return flag_N, logN, sig_logN
