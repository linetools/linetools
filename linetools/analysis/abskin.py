""" Utlities for kinematic analysis of absorption lines
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from astropy import units as u
from astropy.convolution import convolve, Box1DKernel


def generate_stau(velo, flux, sig, kbin=22.*u.km/u.s, debug=False):
    """ Generate the smoothed tau array for kinematic tests

    Parameters
    ----------
    velo : Quantity array (usually km/s)
    flux : Quantity array (flux)
    sig :  Quantity array (sig)
    kbin : Quantity (velocity), optional
      Kernel size for Gaussian smoothing of optical depth array
    debug : bool, optional

    Returns
    -------
    stau : array
       Smoothed tau array
    """
    # Velocity array
    npix = len(velo)
    pix = np.arange(npix).astype(int)

    # Calculate dv
    dv = np.abs(np.median(velo-np.roll(velo,1)))

    # Test for bad pixels
    badzero = np.where((flux == 0) | (sig <= 0))[0]
    if len(badzero) > 0:
        if np.max(badzero)-np.min(badzero) >= 5:
            raise ValueError('orig_kin: too many or too large sections of bad data')

        flux[badzero] = np.mean(np.array([flux[np.min(badzero)-1],
                                                    flux[np.max(badzero)+1]]))
        pdb.set_trace() # Should add sig too

    # Generate the tau array
    tau = np.zeros(npix)
    gd = np.where((flux > sig/2.) & (sig > 0.) )
    if len(gd) == 0:
        raise ValueError('orig_kin: Profile too saturated.')

    tau[gd] = np.log(1./flux[gd])
    sat = (pix == pix)
    sat[gd] = False
    tau[sat] = np.log(2./sig[sat])

    # Smooth
    nbin = (np.round(kbin/dv)).value
    try:
        kernel = Box1DKernel(nbin, mode='center')
    except:
        pdb.set_trace()
    stau = convolve(tau, kernel, boundary='fill', fill_value=0.)
    if debug is True:
        try:
            from xastropy.xutils import xdebug as xdb
        except ImportError:
            pdb.set_trace()
        else:
            xdb.xplot(velo, tau, stau)
            xdb.set_trace()

    # Return
    return stau


def pw97_kin(velo, stau, per=0.05, debug=False):
        """ Measure a standard suite of absorption line kinematics
        from Prochaska & Wolfe 1997

        Parameters
        ----------
        velo : Quantity array
        stau : array

        Returns
        -------
        data : dict
        """
        kin_data = {}
        # Dv (usually dv90)
        tottau = np.sum(stau)
        cumtau = np.cumsum(stau) / tottau
        lft = (np.where(cumtau > per)[0])[0]
        rgt = (np.where(cumtau > (1.-per))[0])[0] - 1
        kin_data['Dv'] = np.round(np.abs(velo[rgt]-velo[lft]))  # Nearest km/s

        if debug:
            pdb.set_trace()

        # Mean/Median
        vcen = (velo[rgt]+velo[lft])/2.
        mean = kin_data['Dv']/2.
        imn = np.argmin( np.fabs(cumtau-0.5) )
        kin_data['fmm'] = (np.abs( (velo[imn]-vcen)/mean )).value

        # fedg
        imx = np.argmax(stau)
        kin_data['fedg'] = (np.abs( (velo[imx]-vcen) / mean )).value

        # Return
        return kin_data


def cgm_kin(velo, stau, per=0.05, debug=False, cov_thresh=0.5,
            dv_zeropk=15.*u.km/u.s, do_orig_kin=False, get_stau=False):#, **kwargs):
    """ Some new tests, invented in the context of CGM studies.
    Some are thanks to John Forbes.

    Parameters
    ----------
    spec
    per
    debug
    cov_thresh : float, optional
      Parameter for the X_fcover test
    dv_zeropk
    do_orig_kin
    get_stau
    kwargs

    Returns
    -------

    """
    kdata = {}

    # voff -- Velocity centroid of profile relative to zsys
    kdata['delta_v'] = np.sum(velo*stau) / np.sum(stau)

    # X "Covering" test
    tottau = np.sum(stau)
    cumtau = np.cumsum(stau) / tottau
    lft = (np.where(cumtau > per)[0])[0]
    rgt = (np.where(cumtau > (1.-per))[0])[0] - 1

    inpix = range(lft,rgt+1)
    tau_covering = np.mean(stau[inpix])
    i_cover = np.where(stau[inpix] > cov_thresh*tau_covering)[0]

    kdata['X_fcover'] = float(len(i_cover)) / float(len(inpix))

    # Peak -- Peak optical depth velocity
    imx = np.argmax(stau)
    kdata['v_peak'] = velo[imx]

    # ###
    # Zero peak -- Ratio of peak optical depth to that within 15 km/s of zero
    tau_zero = stau[imx]
    zpix = np.where( np.abs(velo) < dv_zeropk)[0]
    if len(zpix) == 0:
        kdata['zero_pk'] = 0.
    else:
        mx_ztau = np.max(stau[zpix])
        kdata['zero_pk'] = np.max([0., np.min([mx_ztau/tau_zero,1.])])

    # ###
    # Forbes "Covering"
    dv = np.abs(velo[1]-velo[0])
    forbes_fcover = dv * np.sum(stau) / tau_zero
    kdata['JF_fcover'] = forbes_fcover

    # Return
    return kdata
