""" Functions related to convolution.
      Taken from Barak by JXP
      May replace with scipy functions
      NT: replace barak routines with astropy.convolution equivalents
""" 

# py2.6+ compatibility
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from astropy.convolution import convolve, Gaussian1DKernel, CustomKernel

def convolve_psf_old(a, fwhm, edge='invert', replace_nan=False, debug=False):
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

def convolve_window_old(a, window, edge='invert'):
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


#Updated functions using astropy.convolution
def convolve_psf(array, fwhm, boundary='fill', fill_value=0., normalize_kernel=False):
    """ Convolve an array with a gaussian kernel.

    Given an array of values `a` and a gaussian full width at half
    maximum `fwhm` in pixel units, returns the convolution of the
    array with the gaussian kernel.

    Parameters
    ----------
    array : array, shape(N,)
        Array to convolve
    fwhm : float
        Gaussian full width at half maximum in pixels.
    boundary : str, optional
        A flag indicating how to handle boundaries:
            * `None`
                Set the ``result`` values to zero where the kernel
                extends beyond the edge of the array (default).
            * 'fill'
                Set values outside the array boundary to ``fill_value``.
            * 'wrap'
                Periodic boundary that wrap to the other side of ``array``.
            * 'extend'
                Set values outside the array to the nearest ``array``
                value.
    fill_value : float, optional
        The value to use outside the array when using boundary='fill'
    normalize_kernel : bool, optional
        Whether to normalize the kernel prior to convolving

    Returns
    -------
    convolved_array : array, shape (N,)
    
    Notes
    -----
    This function uses astropy.convolution 
    """

    const2   = 2.354820046             # 2*sqrt(2*ln(2))
    const100 = 3.034854259             # sqrt(2*ln(100))
    sigma = fwhm / const2
    # gaussian drops to 1/100 of maximum value at x =
    # sqrt(2*ln(100))*sigma, so number of pixels to include from
    # centre of gaussian is:
    n = np.ceil(const100 * sigma)
    x_size = int(2*n) + 1 # we want this to be odd integer
    return convolve(array, Gaussian1DKernel(sigma,x_size=x_size),boundary=boundary,fill_value=fill_value,normalize_kernel=normalize_kernel)

def convolve_window(array, window, boundary='fill', fill_value=0., normalize_kernel=True):
    """ Convolve an array with an arbitrary window kernel.

    Parameters
    ----------
    array : array, shape (N,)
        Array to convolve
    window : array, shape (M,)
        The window array must have an odd number of elements.
    boundary : str, optional
        A flag indicating how to handle boundaries:
            * `None`
                Set the ``result`` values to zero where the kernel
                extends beyond the edge of the array (default).
            * 'fill'
                Set values outside the array boundary to ``fill_value``.
            * 'wrap'
                Periodic boundary that wrap to the other side of ``array``.
            * 'extend'
                Set values outside the array to the nearest ``array``
                value.
    fill_value : float, optional
        The value to use outside the array when using boundary='fill'
    normalize_kernel : bool, optional
        Whether to normalize the kernel prior to convolving

    Returns
    -------
    convolved_array : array, shape (N,)
    
    Notes
    -----
    This function uses astropy.convolution, and astropy.modeling

    If window kernel is as large or larger than the array, it raises an error.

    """
    if not len(window) % 2:
        raise ValueError('`window` must have an odd number of elements!')

    if len(window) >= len(array):
        raise ValueError('`window` is too big for the `array`!')

    return convolve(array, CustomKernel(window),boundary=boundary,fill_value=fill_value,normalize_kernel=normalize_kernel)
