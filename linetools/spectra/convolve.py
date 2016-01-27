""" Functions related to convolution.
""" 

# py2.6+ compatibility
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from astropy.convolution import convolve, Gaussian1DKernel, CustomKernel

#Updated functions using astropy.convolution
def convolve_psf(array, fwhm, boundary='fill', fill_value=0.,
                 normalize_kernel=True):
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
    return convolve(array, Gaussian1DKernel(sigma, x_size=x_size),
                    boundary=boundary, fill_value=fill_value,
                    normalize_kernel=normalize_kernel)

