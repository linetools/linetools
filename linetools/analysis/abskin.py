""" Utlities for kinematic analysis of absorption lines
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import warnings

from astropy import units as u
from astropy import constants as const
from astropy.convolution import convolve, Box1DKernel

from linetools.spectra.xspectrum1d import XSpectrum1D


########################## ##########################
def generate_stau(velo, flux, sig, kbin=22.*u.km/u.s, debug=False):
    """ Generate the smoothed tau array for kinematic tests

    Parameters
    ----------
    velo : Quantity array (usually km/s)
    flux : Quantity array (flux)
    sig :  Quantity array (sig)

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
    badzero=np.where((flux == 0) & (sig <= 0))[0]
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
    kernel = Box1DKernel(nbin, mode='center')
    stau = convolve(tau, kernel, boundary='fill', fill_value=0.)
    if debug is True:
        from xastropy.xutils import xdebug as xdb
        xdb.xplot(velo, tau, stau)

    # Return
    return stau


def pw97_kin(velo, stau, per=0.05):
        """ Measure a standard suite of absorption line kinematics
        from Prochaska & Wolfe 1997

        Parameters
        ----------
        velo : Quantity array
        stau : array

        Returns
        -------
        kin_data : dict
        """
        kin_data = {}
        # Dv (usually dv90)
        tottau = np.sum(stau)
        cumtau = np.cumsum(stau) / tottau
        lft = (np.where(cumtau > per)[0])[0]
        rgt = (np.where(cumtau > (1.-per))[0])[0] - 1
        kin_data['Dv'] = np.round(np.abs(velo[rgt]-velo[lft]))

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
