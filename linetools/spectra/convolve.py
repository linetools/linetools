""" Functions related to convolution.
      Taken from Barak by JXP
      May replace with scipy functions
""" 

# py2.6+ compatibility
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

def convolve_psf(a, fwhm, edge='invert', replace_nan=False, debug=False):
    """ Convolve an array with a gaussian window.

    Given an array of values `a` and a gaussian full width at half
    maximum `fwhm` in pixel units, returns the convolution of the
    array with the normalised gaussian.

    Parameters
    ----------
    a : array, shape(N,)
      Array to convolve
    fwhm : float
      Gaussian full width at half maximum in pixels. This should be > 2
      to sample the gaussian PSF properly.
    replace_nan : bool, optional
      Replace NAN values [NOT YET IMPLEMENTED]

    Returns
    -------
    convolved_a : array, shape (N,)
    
    Notes
    -----
    The Gaussian kernel is calculated for as many pixels required
    until it drops to 1% of its peak value. The data will be spoiled
    at distances `n`/2 (rounding down) from the edges, where `n` is
    the width of the Gaussian in pixels.
    """
    const2   = 2.354820046             # 2*sqrt(2*ln(2))
    const100 = 3.034854259             # sqrt(2*ln(100))
    sigma = fwhm / const2
    # gaussian drops to 1/100 of maximum value at x =
    # sqrt(2*ln(100))*sigma, so number of pixels to include from
    # centre of gaussian is:
    n = np.ceil(const100 * sigma)
    #if replace_nan:
    #    a = nan2num(a, replace='interp')
    if debug:
        print("First and last {0} pixels of output will be invalid".format(n))
    x = np.linspace(-n, n, 2*n + 1)        # total no. of pixels = 2n+1
    gauss = np.exp(-0.5 * (x / sigma) ** 2 )

    return convolve_window(a, gauss, edge=edge)

def convolve_window(a, window, edge='invert'):
    """ Convolve an array with an arbitrary window.

    Parameters
    ----------
    a : array, shape (N,)
    window : array, shape (M,)
      The window array should have an odd number of elements.
    edge : {'invert', 'reflect', 'extend'} or int  (default 'invert')    
      How to mitigate edge effects. If 'invert', the edges of `a` are
      extended by inversion, similarly for reflection. 'extend' means
      the intial and final points are replicated to extend the
      array. An integer value means take the median of that many
      points at each end and extend by replicating the median value.

    Returns
    -------
    convolved_a : array, shape (N,)

    Notes
    -----
    The window is normalised before convolution.
    """
    npts = len(window)
    if not npts % 2:
        raise ValueError('`window` must have an odd number of elements!')

    n = npts // 2

    # normalise the window
    window /= window.sum()

    # Add edges to either end of the array to reduce edge effects in
    # the convolution.
    if len(a) < 2*n:
        raise ValueError('Window is too big for the array! %i %i' % (len(a), n))
    if edge == 'invert':
        temp1 = 2*a[0] - a[n:0:-1], a, 2*a[-1] - a[-2:-n-2:-1]
    elif edge == 'reflect':
        temp1 =  a[n:0:-1], a, a[-2:-n-2:-1]
    elif edge == 'extend':
        temp1 =  a[0] * np.ones(n) , a, a[-1] * np.ones(n)
    else:
        try:
            abs(int(edge))
        except TypeError:
            raise ValueError('Unknown value for edge keyword: %s' % edge)
        med1 = np.median(a[:edge])
        med2 = np.median(a[-edge:])
        temp1 =  med1 * np.ones(n) , a, med2 * np.ones(n)

    temp2 = np.convolve(np.concatenate(temp1), window, mode=str('same'))

    return temp2[n:-n]

def convolve_constant_dv(wa, fl, wa_dv=None, npix=4., vfwhm=None, dv_const=False):
    """ Convolve a wavelength array with a gaussian of constant
    velocity width.

    If `vfwhm` is specified an intermediate wavelength array with
    constant velocity pixel width is calculated. Otherwise, both
    `wa_dv` and `npix` must be given -- this is faster because no
    intermediate array needs to be calculated.

    Parameters
    ----------
    fl, wa : arrays of floats, length N
      The array to be convolved and its wavelengths.
    vfwhm : float, optional
      Full width at half maximum in velocity space (km/s) of the
      gaussian kernel with which to convolve `fl`.
    npix : float, default 4
      Number of pixels corresponding to `vfwhm` in `wa_dv` if given,
      otherwise `wa` is interpolated to an array with velocity pixel
      width = vfwhm / npix.
    wa_dv : array of floats, default `None`
      Wavelength array with a constant velocity width (this can be
      generated with make_constant_dv_wa_scale()).
    dv_const : Bool, default 'False'
      Specifies whether the input array has constant velocity

    Returns
    -------
    fl_out : array of length N
      fl convolved with the gaussian kernel with the specified FWHM.     

    Examples
    --------
    >>> from barak import sed
    >>> wa = np.arange(5000, 7000, 0.03)
    >>> fl = np.random.randn(len(wa)) + 1.
    >>> instrument_profile_fwhm = 5.5 # resolution fwhm km/s
    >>> flsmooth = convolve_constant_dv(wa, fl, vfwhm=instrument_profile)

    To avoid generating a new constant velocity scale for multiple
    convolutions of the same spectrum, you can create your own and
    pass it as the wa_dv argument:

    >>> from barak import sed
    >>> subpixkms = 1.3  # km/s    
    >>> wa_dv = sed.make_constant_dv_wa_scale(wa[0], wa[-1], subpixkms)
    >>> npix = instrument_profile / subpixkms
    >>> flsmooth = convolve_constant_dv(wa, fl, wa_dv=wa_dv, npix=npix)
    """
    # interpolate to the log-linear scale, convolve, then
    # interpolate back again.
    if vfwhm is not None:
        wa_dv = make_constant_dv_wa_scale(wa[0], wa[-1], float(vfwhm)/npix)
    if dv_const is False:
        fl_dv = np.interp(wa_dv, wa, fl)
    else: fl_dv = fl
    # Convolve
    fl_dv_smoothed = convolve_psf(fl_dv, npix)
    # Interpolate as needed
    if dv_const is False:
        fl_out = np.interp(wa, wa_dv, fl_dv_smoothed)
    else: fl_out = fl_dv_smoothed
    # Return
    return fl_out
