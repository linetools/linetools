""" Interpolation-related tools.
""" 

# p2.6+ compatibility
from __future__ import division, print_function, unicode_literals
import numpy as np

class AkimaSpline(object):
    """ Describes an Akima Spline through a set of points.

    It must be instantiated with a set of `xvals` and `yvals` knot values,
    and then can be called with a new set of x values `x`. This is
    used by `interp_Akima`, see its documentation for more
    information.

    Parameters
    ----------
    xvals, yvals : array_like, shape (N,)
      Reference values. xvals should not contain duplicates.

    References
    ----------
    "A new method of interpolation and smooth curve fitting based
    on local procedures." Hiroshi Akima, J. ACM, October 1970, 17(4),
    589-602.

    Notes
    -----
    This is adapted from a function written by `Christoph Gohlke
    <http://www.lfd.uci.edu/~gohlke/>`_ under a BSD license:

    Copyright (c) 2007-2012, Christoph Gohlke
    Copyright (c) 2007-2012, The Regents of the University of California
    Produced at the Laboratory for Fluorescence Dynamics
    All rights reserved.
    """
    def __init__(self, xvals, yvals):
        """
        """

        x = np.asarray(xvals, dtype=np.float64)
        y = np.asarray(yvals, dtype=np.float64)
        if x.ndim != 1:
            raise ValueError("x array must be one dimensional")
     
        n = len(x)
        if n < 3:
            raise ValueError("Array too small")
        if n != len(y):
            raise ValueError("Size of x-array must match data shape")
     
        dx = np.diff(x)
        if (dx <= 0.0).any():
            isort = np.argsort(x)
            x = x[isort]
            y = y[isort]
            dx = np.diff(x)
            if (dx == 0.).any():
                raise ValueError("x array has duplicate values")

        m = np.diff(y) / dx
        mm = 2. * m[0] - m[1]
        mmm = 2. * mm - m[0]
        mp = 2. * m[n - 2] - m[n - 3]
        mpp = 2. * mp - m[n - 2]
     
        m1 = np.concatenate(([mmm], [mm], m, [mp], [mpp]))
     
        dm = np.abs(np.diff(m1))
        f1 = dm[2:n + 2]
        f2 = dm[0:n]
        f12 = f1 + f2
     
        ids = np.nonzero(f12 > 1e-9 * f12.max())[0]
        b = m1[1:n + 1]
     
        b[ids] = (f1[ids] * m1[ids + 1] + f2[ids] * m1[ids + 2]) / f12[ids]
        c = (3. * m - 2. * b[0:n - 1] - b[1:n]) / dx
        d = (b[0:n - 1] + b[1:n] - 2. * m) / dx ** 2

        self.xvals, self.yvals, self.b, self.c, self.d = x, y, b, c, d

    def __call__(self, x):
        """
        Parameters
        ----------
        x : array_like, shape (M,)
          Values at which to interpolate.

        Returns
        -------
        vals : ndarray, shape (M,)
           Interpolated values.
        """
        x = np.asarray(x, dtype=np.float64)

        if x.ndim != 1:
            raise ValueError("Array must be one dimensional")

        c0 = x < self.xvals[0]
        c2 = x > self.xvals[-1]
        c1 = ~(c0 | c2)
        if c1.sum() == 0:
            raise ValueError('All points are outside the interpolation range!')
        x1 = x[c1]
        out = np.empty_like(x)
        bins = np.digitize(x1, self.xvals)
        bins = np.minimum(bins, len(self.xvals) - 1) - 1
        b = bins[0:len(x1)]
        wj = x1 - self.xvals[b]
        out[c1] = ((wj * self.d[b] + self.c[b])*wj + self.b[b])*wj + \
                  self.yvals[b]

        # use linear extrapolation for points outside self.xvals
        if c0.any():
            y = out[c1]
            slope = (y[1] - y[0]) / (x1[1] - x1[0])
            intercept = y[0] - slope * x1[0]
            out[c0] = x[c0] *slope + intercept
        if c2.any():
            y = out[c1]
            slope = (y[-2] - y[-1]) / (x1[-2] - x1[-1])
            intercept = y[-1] - slope * x1[-1]
            out[c2] = x[c2] *slope + intercept

        return out

def interp_Akima(x_new, x, y):
    """Return interpolated data using Akima's method.

    Akima's interpolation method uses a continuously differentiable
    sub-spline built from piecewise cubic polynomials. The resultant
    curve passes through the given data points and will appear smooth
    and natural.

    Parameters
    ----------
    x_new : array_like, shape (M,)
        Values at which to interpolate.
    x, y : array_like, shape (N,)
        Reference values. x cannot contain duplicates.

    Returns
    -------
    vals : ndarray, shape (M,)
       Interpolated values.

    References
    ----------
    "A new method of interpolation and smooth curve fitting based
    on local procedures." Hiroshi Akima, J. ACM, October 1970, 17(4),
    589-602.
    """
    interpolator = AkimaSpline(x, y)
    return interpolator(x_new)
