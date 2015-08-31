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
    def __init__(self, fig, ax, wa, fl, artists):
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

