""" Classes useful for making interactive plots.
"""
from __future__ import division, print_function, unicode_literals, absolute_import

import os
import numpy as np
from ..utils import between
from ..spectra.convolve import convolve_psf
from ..spectra.plotting import get_flux_plotrange
from .interp import AkimaSpline

from astropy.modeling import models
import astropy.units as u

import matplotlib.transforms as mtran
import matplotlib.pyplot as plt


import json

def savejson(filename, obj, overwrite=False, indent=None):
    """ Save a python object to filename using using the JSON encoder."""

    if os.path.lexists(filename) and not overwrite:
        raise IOError('%s exists' % filename)
    if filename.endswith('.gz'):
        fh = gzip.open(filename, 'wt')
    else:
        fh = open(filename, 'wt')
    try:
        json.dump(obj, fh, indent=indent)
    except:
        import pdb; pdb.set_trace()

    fh.close()

def loadjson(filename):
    """ Load a python object saved with savejson."""
    fh = open(filename, 'rt')
    try:
        obj = json.load(fh)
    except:
        import pdb; pdb.set_trace()
        
    fh.close()
    return obj


def local_median(wa, fl, er, x, npix=10, default=None):
    """ find the median flux value at x using +/- npix pixels.
    """
    i = np.searchsorted(wa, x)
    i0,i1 = i - npix, i + npix
    good = (er[i0:i1] > 0) & ~np.isnan(fl[i0:i1])
    if good.sum():
        return np.median(fl[i0:i1][good])
    else:
        return default


class PlotWrapBase(object):
    """ A base class that has all the navigation and smoothing
    keypress events.

    These attributes must be defined in a subclass:

    self.wa, self.fl    Spectrum wavelength and flux
    self.nsmooth        integer > 0 that determines the smoothing
    self.ax             Axes where spectrum is plotted
    self.fig            Figure which holds the axes.
    self.artists['fl']  The Matplotlib line artist that represents the flux.


    The keypress events need to be connected to the figure with
    something like:

    def connect(self, fig):
        cids = dict(key=[])
        # the top two are methods of PlotWrapBase
        cids['key'].append(fig.canvas.mpl_connect(
            'key_press_event', self.on_keypress_navigate))
        cids['key'].append(fig.canvas.mpl_connect(
            'key_press_event', self.on_keypress_smooth))
        self.cids.update(cids)
    """
    _help_string = """
i,o          Zoom in/out x limits
y            Zoom out y limits
Y            Guess y limits
t,b          Set y top/bottom limit
l,r          Set left/right x limit
[,]          Pan left/right
w            Plot the whole spectrum

S,U          Smooth/unsmooth spectrum
"""

    def __init__(self):
        pass

    def on_keypress_navigate(self, event):
        """ Process a keypress event. Requires attributes self.ax,
        self.fl, self.wa, self.fig
        """
        # Navigation
        if event.key == 'i' and event.inaxes:
            x0,x1 = self.ax.get_xlim()
            x = event.xdata
            dx = abs(x1 - x0)
            self.ax.set_xlim(x - 0.275*dx, x + 0.275*dx)
            self.fig.canvas.draw()
        elif event.key == 'o' and event.inaxes:
            x0,x1 = self.ax.get_xlim()
            x = event.xdata
            dx = abs(x1 - x0)
            self.ax.set_xlim(x - 0.95*dx, x + 0.95*dx)
            self.fig.canvas.draw()
        elif event.key == 'Y' and event.inaxes:
            y0,y1 = self.ax.get_ylim()
            y = event.ydata
            dy = abs(y1 - y0)
            self.ax.set_ylim(y0 - 0.05*dy, y1 + 0.4*dy)
            self.fig.canvas.draw()
        elif event.key == 'y' and event.inaxes:
            x0,x1 = self.ax.get_xlim()
            y0,y1 = get_flux_plotrange(self.fl[between(self.wa, x0, x1)])
            self.ax.set_ylim(y0, y1)
            self.fig.canvas.draw()
        elif event.key == ']':
            x0,x1 = self.ax.get_xlim()
            dx = abs(x1 - x0)
            self.ax.set_xlim(x1 - 0.1*dx, x1 + 0.9*dx)
            self.fig.canvas.draw()
        elif event.key == '[':
            x0,x1 = self.ax.get_xlim()
            dx = abs(x1 - x0)
            self.ax.set_xlim(x0 - 0.9*dx, x0 + 0.1*dx)
            self.fig.canvas.draw()
        elif event.key == 'w':
            self.ax.set_xlim(self.wa[0], self.wa[-1])
            y0,y1 = get_flux_plotrange(self.fl)
            self.ax.set_ylim(y0, y1)
            self.fig.canvas.draw()
        elif event.key == 'b' and event.inaxes:
            y0, y1 = self.ax.get_ylim()
            self.ax.set_ylim(event.ydata, y1)
            self.fig.canvas.draw()
        elif event.key == 't' and event.inaxes:
            y0, y1 = self.ax.get_ylim()
            self.ax.set_ylim(y0, event.ydata)
            self.fig.canvas.draw()
        elif event.key == 'l' and event.inaxes:
            x0, x1 = self.ax.get_xlim()
            self.ax.set_xlim(event.xdata, x1)
            self.fig.canvas.draw()
        elif event.key == 'r' and event.inaxes:
            x0, x1 = self.ax.get_xlim()
            self.ax.set_xlim(x0, event.xdata)
            self.fig.canvas.draw()

    def on_keypress_smooth(self, event):
        """ Smooth the flux with a gaussian. Requires attributes
        self.fl and self.nsmooth, self.artists['fl'] and self.fig."""

        # maybe should use boxcar smoothing?
        if event.key == 'S':
            if self.nsmooth > 0:
                self.nsmooth += 0.5
            else:
                self.nsmooth += 1
            sfl = convolve_psf(self.fl, self.nsmooth)
            self.artists['fl'].set_ydata(sfl)
            self.fig.canvas.draw()
        elif event.key == 'U':
            self.nsmooth = 0
            self.artists['fl'].set_ydata(self.fl)
            self.fig.canvas.draw()

class PlotWrapNav(PlotWrapBase):
    """ Enable simple XIDL-style navigation for plotting a spectrum.

    For example, i and o for zooming in y direction, [ and ] for
    panning, S and U for smoothing and unsmoothing.
    """
    def __init__(self, fig, ax, wa, fl, artists, printhelp=True):
        """
        Parameters
        ----------
        fig : matplotlib Figure
        ax : matplotlib axes
          The Axes where the spectrum is plotted.
        wa, fl : array
          Wavelength and flux arrays
        artists : dict
          A dictionary which must contain a key 'fl', which is the
          matplotlib artist corresponding to the flux line.
        printhelp : bool, optional
          Whether to print a help message when first called.
        """
        self.artists = artists
        self.fig = fig
        self.ax = ax
        if isinstance(wa, u.Quantity):
            wa = wa.value
        self.wa = wa
        if isinstance(fl, u.Quantity):
            fl = fl.value
        self.fl = fl
        self.nsmooth = 0
        # disable existing keypress events (like 's' for save).
        cids = list(fig.canvas.callbacks.callbacks['key_press_event'])
        for cid in cids:
            fig.canvas.callbacks.disconnect(cid)
        self.cids = {}
        self.connect()
        if printhelp:
            print(self._help_string)


    def on_keypress(self, event):
        """ Print a help message"""
        if event.key == '?':
            print(self._help_string)

    def connect(self):
        cids = dict(key=[])
        # the top two are methods of PlotWrapBase
        cids['key'].append(self.fig.canvas.mpl_connect(
            'key_press_event', self.on_keypress))
        cids['key'].append(self.fig.canvas.mpl_connect(
            'key_press_event', self.on_keypress_navigate))
        cids['key'].append(self.fig.canvas.mpl_connect(
            'key_press_event', self.on_keypress_smooth))
        self.cids.update(cids)

class InteractiveCoFit(PlotWrapBase):
    help_message = PlotWrapBase._help_string + """
a        : add a new spline knot
A        : add a new spline knot, and use a flux median to guess y position
+        : double the number of spline knots
_        : halve the number of spline knots
d        : delete the nearest knot
m        : move the nearest knot
M        : move the nearest knot, and use a flux median to guess y position

q        : quit
"""
    def __init__(self, wa, fl, er, contpoints, co=None,
                 fig=None, anchor=None):
        """ Initialise figure, plots and variables.

        Parameters
        ----------
        wa : Wavelengths
        fl : Fluxes
        er : One sigma errors
        contpoints : list of x,y tuple pairs (None)
            The points through which a cubic spline is passed,
            defining the continuum.
        co : Continuum, optional  
            The existing continuum, if one is already defined.
        anchor : bool
            Whether to prevent modification of the first and last
            spline knots. Default is None, which means anchor only if
            co is given.

        Notes
        -----
        Updates the following attributes:

         self.wa, self.fl, self.er :  wa, fl, er
         self.contpoints :  Points used to define the continuum.
         self.artists :  Dictionary of matplotlib plotting artists.
         self.connections :  Callback connections.
         self.fig :  The plotting figure instance.
        """
        #setup

        self.nsmooth = 0
        self.wa = wa
        self.fl = fl
        self.er = er
        self.anchor = anchor

        if os.path.lexists('./_knots.jsn'):
            c = raw_input('knots file exists, use this? (y) ')
            if c.lower() != 'n':
                contpoints = loadjson('./_knots.jsn')
        contpoints = sorted(tuple(cp) for cp in contpoints)

        if co is not None:
            self.continuum = np.array(co, copy=True)
            if self.anchor is None:
                self.anchor = True
        else:
            xpts, ypts = zip(*contpoints)
            if len(contpoints) >= 5:
                # need 5 points to define an Akima Spline
                spline = AkimaSpline(xpts, ypts)
                self.continuum = spline(wa)
            else:
                self.continuum = np.interp(wa, xpts, ypts)
            co = self.continuum
            if self.anchor is None:
                self.anchor = False

        if self.anchor:
            wmin = contpoints[0][0]
            wmax = contpoints[-1][0]
        else:
            wmin = wa[0]
            wmax = wa[-1]

        # add extra anchor points so the slopes match at each end of
        # the fitting region.

        i1, i2 = wa.searchsorted([wmin, wmax])
        if i1 == 0:
            i1 = 1
        if i2 == len(wa) or i2 < 0:
            i2 = len(wa) - 1
        x,y = contpoints[0]
        contpoints[0] = wa[i1], y
        x,y = contpoints[-1]
        contpoints[-1] = wa[i2], y
        self.indices = i1, i2
        if self.anchor:
            self.anchor_start = wa[i1 - 1], co[i1 - 1]
            self.anchor_end = wa[i2 + 1], co[i2 + 1]
        self.contpoints = contpoints
        self.wmin = wmin
        self.wmax = wmax

        self.artists = {}
        if fig is None:
            self.fig = plt.figure()
        else:
            self.fig = fig
        # disable any existing key press callbacks
        cids = list(self.fig.canvas.callbacks.callbacks['key_press_event'])
        for cid in cids:
            self.fig.canvas.callbacks.disconnect(cid)

        self.connections = []
        self.finished = False
        self.plotinit()
        self.update()
        self.modifypoints()
        plt.draw()

    def plotinit(self):
        """ Set up the figure and do initial plots.

        Updates the following attributes:

          self.artists
        """
        wa,fl,er = self.wa, self.fl, self.er

        # axis for spectrum & continuum
        a0 = self.fig.add_axes((0.05,0.1,0.9,0.6))
        self.ax = a0
        a0.set_autoscale_on(0)
        # axis for residuals
        a1 = self.fig.add_axes((0.05,0.75,0.9,0.2),sharex=a0)
        a1.set_autoscale_on(0)
        a1.axhline(0, color='k', alpha=0.7, zorder=99)
        a1.axhline(1, color='k', alpha=0.7, zorder=99)
        a1.axhline(-1, color='k', alpha=0.7, zorder=99)
        a1.axhline(2, color='k', linestyle='dashed', zorder=99)
        a1.axhline(-2, color='k', linestyle='dashed', zorder=99)
        m0, = a1.plot([0],[0],'.r',marker='.', mec='none', lw=0, mew=0,
                      ms=6, alpha=0.5)
        a1.set_ylim(-4, 4)
        a0.axhline(0, color='0.7')

        i0,i1 = self.indices
        art = []
        art.append(a0.axvline(wa[i0], color='r', ls='--', lw=2, zorder=10))
        art.append(a0.axvline(wa[i1], color='r', ls='--', lw=2, zorder=10))
        self.artists['indices'] = art
        a0.plot(wa, self.continuum, color='k', lw=2, ls='dashed', zorder=3)
        self.artists['fl'], = a0.plot(wa, fl, lw=1, color='0.7',
                                      linestyle='steps-mid')
        a0.plot(wa, er, lw=0.5, color='orange')
        m1, = a0.plot([0], [0], 'r', zorder=4, lw=2)
        m2, = a0.plot([0], [0], 'o', mfc='None', mew=2, ms=12, mec='r',
                      alpha=0.7)

        a0.set_xlim(wa[i0], wa[i1])
        good = (er[i0:i1] > 0) & ~np.isnan(fl[i0:i1]) & ~np.isinf(fl[i0:i1])
        ymax = 2 * np.abs(np.percentile(fl[i0:i1][good], 95))
        a0.set_ylim(-0.1*ymax, ymax)

        # for histogram
        trans = mtran.blended_transform_factory(a1.transAxes, a1.transData)
        hist, = a1.plot([], [], color='k', transform=trans)
        x = np.linspace(-3,3)

        g = models.Gaussian1D(amplitude=0.05, mean=0, stddev=1)
        a1.plot(g(x), x, color='k', transform=trans, lw=0.5)

        self.fig.canvas.draw()
        self.artists.update(contpoints=m2, cont=m1, resid=m0, hist_left=hist)

        self.finished = False


    def update(self):
        """ Calculates the new continuum, residuals and updates the plots.

        Updates the following attributes:

          self.artists
          self.continuum
        """
        wa,fl,er = self.wa, self.fl, self.er
        co = self.continuum
        if self.anchor:
            cpts = [self.anchor_start] + self.contpoints + [self.anchor_end]
        else:
            cpts = self.contpoints
        i,j = self.indices
        xpts, ypts = zip(*cpts)
        if len(cpts) >= 5:
            # need 5 points to define an Akima Spline
            spline = AkimaSpline(xpts, ypts)
            co[i:j] = spline(wa[i:j])
        else:
            co[i:j] = np.interp(wa[i:j], xpts, ypts)

        resid = (fl[i:j] - co[i:j]) / er[i:j]
        # histogram
        bins = np.arange(0, 5 + 0.1, 0.2)
        w0,w1 = self.fig.axes[1].get_xlim()
        i,j = self.indices
        x,_ = np.histogram(resid[between(wa[i:j], w0, w1)],
                           bins=bins)
        b = np.repeat(bins, 2)
        X = np.concatenate([[0], np.repeat(x,2), [0]])
        Xmax = X.max()
        X = 0.05 * X / Xmax
        self.artists['hist_left'].set_data(X, b)

        if self.anchor:
            xpts, ypts = zip(*self.contpoints[1:-1])
        else:
            xpts, ypts = zip(*self.contpoints)

        self.artists['contpoints'].set_data((xpts, ypts))
        self.artists['cont'].set_data(wa[i:j], co[i:j])
        self.artists['resid'].set_data(wa[i:j], resid)
        self.continuum = co
        savejson('_knots.jsn', self.contpoints, overwrite=True)
        self.fig.canvas.draw()

    def on_keypress(self, event):
        """ Interactive fiddling via the keyboard

        Updates:

         self.contpoints
        """
        if event.key == 'q':
            self.finished = True
            plt.close()
            return
        if event.key == '+':
            # double the number of knots
            xc, yc = zip(*self.contpoints)
            xa0, ya0 = self.contpoints[0]
            xnew = []
            xnew.extend(np.array(xc[:-1]) + 0.5*np.diff(xc))
            ynew = np.interp(xnew, xc, yc)
            ynew = [local_median(self.wa, self.fl, self.er, xnew[i],
                                 default=ynew[i])
                    for i in range(len(xnew))]
            # add to contpoints
            self.contpoints.extend(zip(xnew, ynew))
            self.contpoints.sort()
            self.update()
        if event.key == '_':
            # remove (roughly) half the number of knots
            cp = self.contpoints
            if len(cp) < 2:
                print("Too few spline knots.")
                return
            self.contpoints = [cp[0]] + cp[1:-1][1::2] + [cp[-1]]
            self.update()
        if event.inaxes != self.fig.axes[0]:
            return

        if event.key in ('a', '3'):
            if not (self.wmin < event.xdata < self.wmax):
                print('Outside fitting region')
                return
            # add a point to contpoints
            x, y = event.xdata, event.ydata
            if not self.contpoints or x not in list(zip(*self.contpoints))[0]:
                self.contpoints.append((x, y))
                self.contpoints.sort()
                self.update()
        if event.key == 'A':
            # add a point to contpoints, estimating via median
            if not (self.wmin < event.xdata < self.wmax):
                print('Outside fitting region')
                return
            x = event.xdata
            if not self.contpoints or x not in list(zip(*self.contpoints))[0]:
                y = local_median(self.wa, self.fl, self.er, x,
                                 default=event.ydata)
                self.contpoints.append((x, y))
                self.contpoints.sort()
                self.update()
        elif event.key in ('d', '4'):
            # remove a point from contpoints
            if len(self.contpoints) < 2:
                print('Need at least 1 spline knot')
                return

            contx,conty = zip(*self.ax.transData.transform(self.contpoints))
            sep = np.hypot(event.x - np.array(contx),
                           event.y - np.array(conty))
            ind = sep.argmin()
            if ind in (0, len(sep) - 1) and self.anchor:
                print('Cannot remove anchor knots')
                return
            self.contpoints.remove(self.contpoints[ind])
            self.update()
        elif event.key in ('m', 'M'):
            # Move a point
            if not between(event.xdata, self.wmin, self.wmax):
                print('Outside fitting region')
                return
            
            contx,conty = zip(*self.ax.transData.transform(self.contpoints))
            sep = np.hypot(event.x - np.array(contx),
                           event.y - np.array(conty))
            ind = np.argmin(sep)
            if self.anchor and ind == 0:
                if len(self.contpoints) > 2:
                    ind = 1
                else:
                    print('Cannot move anchor knots')
                    return
            elif self.anchor and ind == len(sep) - 1:
                if len(self.contpoints) > 2:
                    ind = len(sep) - 2
                else:
                    print('Cannot move anchor knots')
                    return

            x, y = event.xdata, event.ydata
            # if M, get y value from a local_median
            if event.key == 'M' and \
                   (not self.contpoints or
                    x not in list(zip(*self.contpoints))[0]):
                y = local_median(self.wa, self.fl, self.er, x,
                                 default=event.ydata)
            self.contpoints[ind] = x, y
            self.contpoints.sort()
            self.update()

        elif event.key == '?':
            print(self.help_message)

    def on_button_release(self, event):
        self.update()

    def modifypoints(self):
        """ Add/remove continuum points."""
        print(self.help_message)
        id1 = self.fig.canvas.mpl_connect('key_press_event',self.on_keypress)
        id2 = self.fig.canvas.mpl_connect('key_press_event',self.on_keypress_smooth)
        id3 = self.fig.canvas.mpl_connect('key_press_event',self.on_keypress_navigate)
        id4 = self.fig.canvas.mpl_connect('button_release_event',self.on_button_release)
        self.connections.extend([id1, id2, id3, id4])
