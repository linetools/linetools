""" Utlities for the analysis of absorption lines
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import warnings

from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from astropy.utils import isiterable
from linetools.lists.linelist import LineList

# Atomic constant
atom_cst = (const.m_e.cgs*const.c.cgs / (np.pi * 
    (const.e.esu**2).cgs)).to(u.AA*u.s/(u.km*u.cm**2))

# e2/(me*c) in CGS
e2_me_c_cgs = (const.e.esu**2 / (const.c.to('cm/s') * const.m_e.to('g')))


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

    # dv (np.abs -- just in case the data weren't sorted)
    delv = np.abs(velo - np.roll(velo,1))
    delv[0] = delv[1]
    delv[-1] = delv[-2]

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
    tvar = np.sum((delv*cst*sig/fx)**2)

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


def get_tau0(wrest, fosc, N, b):
    """Get the value of the optical depth at the line center,
    tau0. Taken from Draine 2011 (see Chapter 9). It neglects stimulated
    emission which is fine for IGM or ISM except for radio-frequency
    transitions.

    Parameters
    ----------
    wrest : Quantity
        Rest-frame wavelength of the transition
    fosc : float
        Oscillator strength of the transition
    N : Quantity or Quantity array
        Column density
    b : Quantity or Quantity array of same shape as N
        Doppler parameter

    Returns
    -------
    tau0: float or array
        Optical depth at the line center. If N and b are
        arrays they must be of same shape.
    """
    # check format for N and b
    if isiterable(N):
        if np.shape(N) != np.shape(b):
            raise IOError('If N is array, b must be array of same shape.')

    # convert to CGS
    b_cgs = b.to('cm/s')
    wrest_cgs = wrest.to('cm')
    N_cgs = N.to('1/cm2')

    # tau0
    tau0 = np.sqrt(np.pi) * e2_me_c_cgs * N_cgs * fosc * wrest_cgs  / b_cgs  # eq. 9.8 Draine 2011
    # check dimensionless
    tau0 = tau0.decompose()
    if tau0.unit != u.dimensionless_unscaled:
        raise IOError('Something went wrong with the units, check input units.')
    return tau0.value


def Wr_from_N_b(N, b, wrest, fosc, gamma):
    """For a given transition with wrest, fosc and gamma, it
    returns the rest-frame equivalent width for a given
    N and b. It uses the approximation given by Draine 2011 book
    (eq. 9.27), which comes from atomic physics considerations
    See also Rodgers & Williams 1974 (NT: could not find the reference
    given by Draine)

    Parameters
    ----------
    N : Quantity or Quantity array
        Column density
    b : Quantity or Quantity array of same shape as N
        Doppler parameter
    wrest : Quantity
        Rest-frame wavelength of the transition
    fosc:  float
        Oscillator strength of the transition
    gamma: Quantity
        Gamma parameter of the transition (usually in s^-1).

    Returns
    -------
    Wr : Quantity
        Rest-frame equivalent width

    Notes
    -----
    See also Wr_from_N_b_transition()

    """

    # first calculate tau0
    tau0 = get_tau0(wrest, fosc, N, b)  # N is only used in tau0

    # convert units to CGS
    b_cgs = b.to('cm/s')
    wrest_cgs = wrest.to('cm')
    c_cgs = const.c.to('cm/s')
    gamma_cgs = gamma.to('1/s')

    # two main regimes (eq. 9.27 of Draine 2011)
    # optically thin:
    W_thin = np.sqrt(np.pi) * (b_cgs/c_cgs) * tau0 / (1 + tau0 / (2 * np.sqrt(2))) # dimensionless
    W_thin = W_thin.decompose()

    # optically thick:
    W2_thick = (2 * b_cgs/c_cgs)**2 * np.log(tau0 / np.log(2)) \
            + (b_cgs/c_cgs) * (wrest_cgs * gamma_cgs / c_cgs) * (tau0 - 1.25393) / np.sqrt(np.pi)  # dimensionless
    W_thick = np.sqrt(W2_thick)  # dimensionless
    W_thick = W_thick.decompose()

    # Check dimensionless
    if (W_thin.unit != u.dimensionless_unscaled) or (W_thick.unit != u.dimensionless_unscaled):
        raise IOError('Something went wrong with the units, check input units.')

    # combined
    W = np.where(tau0 <= 1.25393, W_thin.value, W_thick.value) # dimensionless

    Wr = W * wrest  # in wavelength units (see equation 9.4 of Draine)
    return Wr


def Wr_from_N(N, wrest, fosc):
    """For a given transition with wrest and fosc, it
    returns the rest-frame equivalent width for a given
    N. This is an approximation only valid for tau0 << 1, where
    Wr is independent on Doppler parameter and gamma (see eqs. 9.14 and 9.15 of
    Draine 2011). Please use Wr_from_N_b() if you need an approximation
    valid for a wider range in tau0.

    Parameters
    ----------
    N : Quantity or Quantity array
        Column density
    wrest : Quantity
        Rest-frame wavelength of the transition
    fosc : float
        Oscillator strength of the transition

    Returns
    -------
    Wr : Quantity
        Approximated rest-frame equivalent width valid for the tau0<<1 regime.

    Notes
    -----
    For a better approximation valid for much wider range in tau0,
    please use Wr_from_N_b(). See also Wr_from_N_transition() and Wr_from_N_b_transition().

    """

    # Check units by converting to CGS
    wrest_cgs = wrest.to('cm')
    c_cgs = const.c.to('cm/s')

    # dimensionless Wr
    W = (np.pi * e2_me_c_cgs / c_cgs) * N * fosc * wrest_cgs
    W = W.decompose()
    # Check dimensionless
    if (W.unit != u.dimensionless_unscaled):
        raise IOError('Something went wrong with the units; please check input units.')

    Wr = W * wrest  # in wavelength units (see equation 9.4 of Draine)
    return Wr

def N_from_Wr(Wr, wrest, fosc):
    """For a given transition with wrest and fosc, it
    returns the column density N, for a given rest-frame equivalent width
    Wr. This is an approximation only valid for tau0 << 1, where
    Wr is independent on Doppler parameter and gamma (see eqs. 9.14 and 9.15 of
    Draine 2011).

    Parameters
    ----------
    Wr : Quantity or Quantity array
        Rest-frame equivalent width of the AbsLine
    wrest : Quantity
        Rest-frame wavelength of the transition
    fosc : float
        Oscillator strength of the transition

    Returns
    -------
    N : Quantity
        Approximated column density N only valid for the tau0<<1 regime.

    Notes
    -----
    See also Wr_from_N()

    """

    # Check units by converting to CGS
    wrest_cgs = wrest.to('cm')
    c_cgs = const.c.to('cm/s')

    # dimensionless Wr
    W = Wr / wrest
    N = (c_cgs / (np.pi * e2_me_c_cgs)) * W / ( wrest_cgs * fosc)
    N = N.to((1/u.cm**2))
    return N

def Wr_from_N_b_transition(N, b, transition, llist='ISM'):
    """ For a given transition this function looks
    for the atomic parameters (wa0, fosc, gamma) and returns the
    rest-frame equivalent width for a given N and b. It uses the approximation given by
    Draine 2011 book (eq. 9.27), which comes from atomic physics considerations
    See also Rodgers & Williams 1974 (NT: could not find the reference
    given by Draine)

    Parameters
    ----------
    N : Quantity or Quantity array
        Column density
    b : Quantity or Quantity array of same shape as N
        Doppler parameter
    transition : str
        Name of the transition using linetools' naming
        convention, e.g. 'HI 1215'.
    llist : str
        Name of the linetools.lists.linelist.LineList
        object where to look for the transition name.
        Default is 'ISM', which means the function looks
        within `list = LineList('ISM')`.

    Returns
    -------
    Wr : Quantity
        Rest-frame equivalent width

    Notes
    -----
    See also Wr_from_N_b()

    """
    linelist = LineList(llist)
    transition_dict = linelist[transition]
    if transition_dict is None:
        raise ValueError('Transition {:s} not found within LineList {:s}'.format(transition, linelist.list))
    else:
        # get atomic parameters
        wrest = transition_dict['wrest']
        fosc = transition_dict['f']
        gamma = transition_dict['gamma']

    # return
    return Wr_from_N_b(N, b, wrest, fosc, gamma)


def Wr_from_N_transition(N, transition, llist='ISM'):
    """ For a given transition this function looks
    for the atomic parameters (wa0, fosc) and returns the
    rest-frame equivalent width for a given N. It uses the approximation given by
    Draine 2011 book (eq. 9.15), which is only valid for tau0<<1 where Wr is
    independent of Doppler parameter or gamma. Please use Wr_from_N_b_transition()
    if you need an approximation valid for a wider range in tau0.

    Parameters
    ----------
    N : Quantity or Quantity array
        Column density
    transition : str
        Name of the transition using linetools' naming
        convention, e.g. 'HI 1215', 'CIV 1550', etc.
    llist : str
        Name of the linetools.lists.linelist.LineList
        object where to look for the transition name.
        Default is 'ISM', which means the function looks
        within `list = LineList('ISM')`.

    Returns
    -------
    Wr : Quantity
        Approximated rest-frame equivalent width valid for tau0<<1.

    Notes
    -----
    See also Wr_from_N(), Wr_from_N_b_transition(), Wr_from_N_b()

    """
    linelist = LineList(llist)
    transition_dict = linelist[transition]
    if transition_dict is None:
        raise ValueError('Transition {:s} not found within LineList {:s}'.format(transition, linelist.list))
    else:
        # get atomic parameters
        wrest = transition_dict['wrest']
        fosc = transition_dict['f']

    # return
    return Wr_from_N(N, wrest, fosc)


def N_from_Wr_transition(Wr, transition, llist='ISM'):
    """ For a given transition this function looks
    for the atomic parameters (wa0, fosc) and returns the
    column density for a given rest-frame equivalent width. It uses
    the approximation given by Draine 2011 book (eq. 9.15), which is
    only valid for tau0<<1 where Wr is independent of Doppler
    parameter or gamma.

    Parameters
    ----------
    Wr : Quantity or Quantity array
        Rest-frame wavelength
    transition : str
        Name of the transition using linetools' naming
        convention, e.g. 'HI 1215', 'CIV 1550', etc.
    llist : str
        Name of the linetools.lists.linelist.LineList
        object where to look for the transition name.
        Default is 'ISM', which means the function looks
        within `list = LineList('ISM')`.

    Returns
    -------
    N : Quantity
        Approximated column density only valid for tau0<<1.

    Notes
    -----
    See also Wr_from_N(), Wr_from_N_b_transition(), Wr_from_N_b()

    """
    linelist = LineList(llist)
    transition_dict = linelist[transition]
    if transition_dict is None:
        raise ValueError('Transition {:s} not found within LineList {:s}'.format(transition, linelist.list))
    else:
        # get atomic parameters
        wrest = transition_dict['wrest']
        fosc = transition_dict['f']

    # return
    return N_from_Wr(Wr, wrest, fosc)

