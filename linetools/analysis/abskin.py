""" Utlities for kinematic analysis of absorption lines
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from astropy import units as u
from astropy.convolution import convolve, Box1DKernel

########################## ##########################
########################## ##########################
class AbsKin(object):
    """ Class for kinematics on an absorption line

    Attributes
    ----------
    wrest: Quantity
      Rest wavelength of line analyzed
    zabs : float
      Absorption redshift
    vmnx: Quantity tuple (vmin,vmax)
      Velocity range for analysis
    """

    # Initialize
    def __init__(self, wrest, zabs, vmnx):

        # Absorption system
        self.wrest = wrest
        self.zabs = zabs
        self.vmnx = vmnx

        # Data
        self.data = {}
        self.keys = ['flg', 'Dv', 'fedg', 'fmm', 'delta_v', 'X_fcover',
                     'v_peak', 'zero_pk', 'JF_fcover']
        self.key_dtype = ['i4', 'f4', 'f4', 'f4', 'f4', 'f4',
                          'f4', 'f4', 'f4']

        # Init
        for key in self.keys:
            self.data[key] = 0
        #xdb.set_trace()

    # Access the data and return the value
    def __getitem__(self, item):
        try:
            return self.data[item]
        except KeyError:
           raise KeyError

    ########################## ##########################
    def mk_pix_stau(self, spec, kbin=22.*u.km/u.s, debug=False, **kwargs):
        """ Generate the smoothed tau array for kinematic tests

        Parameters
        ----------
        spec : XSpectrum1D class
          Input spectrum
          velo is expected to have been filled already
        kbin : Quantity (velocity), optional
          Kernel size for Gaussian smoothing of optical depth array

        Returns
        -------
        """
        # Generate velo
        velo = spec.relative_vel((1+self.zabs)*self.wrest)
        pix = (velo >= self.vmnx[0]) & (velo <= self.vmnx[1])

        # Fill
        self.velo = velo
        self.stau = generate_stau(velo[pix], spec.flux[pix],
                                  spec.sig[pix], kbin=kbin)
        self.pix = pix


    ########################## ##########################
    def orig_kin(self, spec, kbin=22., per=0.05, get_stau=False, debug=False, **kwargs):
        """ Measure a standard suite of absorption line kinematics
        from Prochaska & Wolfe 1997

        Parameters
        ----------
        spec: Spectrum1D class
          Input spectrum
        velo is expected to have been filled already
        fill: bool (True)
          Fill the dictionary with some items that other kin programs may need

        Returns
        -------
        out_kin : dict
        Dictionary of kinematic measurements
        """
        kdata = pw97_kin(self.velo[self.pix], self.stau, per=per)
        for key,item in kdata.items():
            self.data[key] = item

        # Set flag
        if (self.data['flg'] % 2) < 1:
            self.data['flg'] = 1


    ########################## ##########################
    def cgm_kin(self, spec, per=0.05, debug=False, cov_thresh=0.5,
                dv_zeropk=15.*u.km/u.s, do_orig_kin=False, get_stau=False, **kwargs):
        """ Some new tests, invented in the context of CGM studies.
        Some are thanks to John Forbes.

        This code is usually run after orig_kin.  You should probably run them
        separately if you plan to modify the default settings of either.

        Parameters
        ----------
        spec: Spectrum1D class
          Input spectrum
        velo is expected to have been filled already
        cov_thresh: float (0.5)
          Parameter for the X_fcover test
        """
        # Generate stau and pix?
        if (self.stau is None) | (get_stau is True):
            self.mk_pix_stau(spec, **kwargs)

        # Original kin?
        if do_orig_kin is True:
            self.orig_kin(spec)

        # voff -- Velocity centroid of profile relative to zsys
        self.data['delta_v'] = np.sum(
            self.velo[self.pix] * self.stau ) / np.sum( self.stau )

        # ###
        # X "Covering" test
        tottau = np.sum( self.stau )
        cumtau = np.cumsum(self.stau) / tottau
        lft = (np.where(cumtau > per)[0])[0]
        rgt = (np.where(cumtau > (1.-per))[0])[0] - 1

        inpix = range(lft,rgt+1)
        tau_covering = np.mean( self.stau[inpix] )
        i_cover = np.where( self.stau[inpix] > cov_thresh*tau_covering)[0]

        self.data['X_fcover'] = float(len(i_cover)) / float(len(inpix))


        # ###
        # Peak -- Peak optical depth velocity
        imx = np.argmax(self.stau)
        self.data['v_peak'] = self.velo[self.pix[imx]]

        # ###
        # Zero peak -- Ratio of peak optical depth to that within 15 km/s of zero
        tau_zero = self.stau[imx]
        if (self.vmnx[0] > 0.) | (self.vmnx[1] < 0.):
            #; Not covered
            #; Assuming zero value
            self.data['zero_pk'] = 0.
        else:
            zpix = np.where( np.abs(self.velo[self.pix]) < dv_zeropk)[0]
            if len(zpix) == 0:
                raise ValueError('cgm_kin: Problem here..')
            mx_ztau = np.max(self.stau[zpix])
            self.data['zero_pk'] = np.max([0. , np.min( [mx_ztau/tau_zero,1.])])

        # ###
        # Forbes "Covering"
        dv = np.abs(self.velo[self.pix[1]]-self.velo[self.pix[0]])
        forbes_fcover = dv * np.sum( self.stau ) / tau_zero
        self.data['JF_fcover'] = forbes_fcover

        # Set flag
        if (self.data['flg'] % 4) < 2:
            self.data['flg'] += 2

    def fill_kin(self, spec, **kwargs):
        """ Perform all the measurements
        Parameters
        ----------
        spec : spectrum
        kwargs

        Returns
        -------

        """

        # Setup
        self.mk_pix_stau(spec, **kwargs)
        # Original kinematics
        self.orig_kin(spec, **kwargs)
        # New CGM kinematics
        self.cgm_kin(spec, **kwargs)

    def __repr__(self):
        return ('[{:s}: wrest={:g}, zabs={:g}, data={}]'.format(
                self.__class__.__name__, self.wrest, self.zabs, self.data))


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
    kernel = Box1DKernel(nbin, mode='center')
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
    """
    if (self.vmnx[0] > 0.) | (self.vmnx[1] < 0.):
        #; Not covered
        #; Assuming zero value
        self.data['zero_pk'] = 0.
    else:
    """
    zpix = np.where( np.abs(velo) < dv_zeropk)[0]
    if len(zpix) == 0:
        raise ValueError('cgm_kin: Problem here..')
    mx_ztau = np.max(stau[zpix])
    kdata['zero_pk'] = np.max([0., np.min([mx_ztau/tau_zero,1.])])

    # ###
    # Forbes "Covering"
    dv = np.abs(velo[1]-velo[0])
    forbes_fcover = dv * np.sum(stau) / tau_zero
    kdata['JF_fcover'] = forbes_fcover

    # Return
    return kdata
