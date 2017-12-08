""" Classes for making interactive plots.
"""
from __future__ import division, print_function, unicode_literals, absolute_import


import pdb

import astropy.units as u

from linetools.spectra.xspectrum1d import XSpectrum1D

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

def stack_plot(abslines, vlim=[-300,300.]*u.km/u.s, nrow=6, show=True, spec=None,
               ymnx=(-0.1,1.1), figsz=(18,11), return_fig=False,
               tight_layout=False, add_ew=False, zref=None):
    """Show a stack plot of the input lines
    Assumes the data are normalized.
    Comments about EW:
      If line is blended, upper limit is shown with error
      If EW < 3*sig_EW : upper limit is set to 2*sig_EW and shown
    Parameters
    ----------
    abslines : list
      list of AbsLine objects
    vlim : Quantities
      velocity range for the plot
    nrow : int, optional
      Maximum number of rows per column
    show : bool, optional
      Show the plot?
    ymnx : tuple, optional
      ymin, ymax
    figsz : tuple, optional
      xdim, ydim
    return_fig : bool, optional
      If True, return stackplot as plt.Figure() instance for further manipulation
    tight_layout : bool, optional
      If True, remove whitespace between panels
    add_ew : bool, optional
      If True, add EW values on the figure
    zref : float, optional
      Reference redshift for velocity scale; if None, each AbsLine.z is used

    Returns
    -------
    fig : matplotlib Figure, optional
        Figure instance containing stack plot with subplots, axes, etc.
    """
    mpl.rcParams['font.family'] = 'stixgeneral'
    mpl.rcParams['font.size'] = 15.
    # Check for spec (required)
    gdiline = []
    for iline in abslines:
        if isinstance(iline.analy['spec'], XSpectrum1D):
            gdiline.append(iline)
        else:
            if spec is not None:
                iline.analy['spec'] = spec
                gdiline.append(iline)
    nplt = len(gdiline)
    if nplt == 0:
        print("Load spectra into the absline.analy['spec']")
        return
    # Setup plot
    nrow = min(nplt, nrow)
    ncol = nplt // nrow + (nplt % nrow > 0)
    # Plot
    fig=plt.figure(figsize=figsz)
    plt.clf()
    gs = gridspec.GridSpec(nrow, ncol)
    # Loop me
    for qq, iline in enumerate(gdiline):
        #ax = plt.subplot(gs[qq % nrow, qq//nrow])
        ax = fig.add_subplot(gs[qq % nrow, qq//nrow])
        ax.clear()
        # Normalize as need be (note the spectrum can vary with iline)
        if iline.analy['spec'].co_is_set:
            co = iline.analy['spec'].co
        else:
            co = 1.
        norm_flux = iline.analy['spec'].flux/co
        norm_sig = iline.analy['spec'].sig/co
        #
        # Plot
        if zref is None:
            zref = iline.z
        velo = iline.analy['spec'].relative_vel((1+zref)*iline.wrest)
        ax.plot(velo, norm_flux,  'k-', linestyle='steps-mid')
        ax.plot(velo, norm_sig,  'r:')
        # Lines
        ax.plot([0]*2, ymnx, 'g--')
        # Axes
        ax.set_xlim(vlim.value)
        ax.set_ylim(ymnx)
        #ax.minorticks_on()
        if ((qq+1) % nrow == 0) or ((qq+1) == nplt):
            ax.set_xlabel('Relative Velocity (km/s)')
        else:
            ax.get_xaxis().set_ticks([])
        # Label
        ax.text(0.1, 0.1, iline.data['name'], transform=ax.transAxes, ha='left', va='center')#, fontsize='large')  # , bbox={'facecolor':'white'})
        if add_ew:
            if iline.attrib['flag_EW'] > 0:
                ewval = '{:05.2f}'.format(iline.attrib['EW'].value)
                ewsigval = '{:05.2f}'.format(iline.attrib['sig_EW'].value)
                compsign = '='
                if iline.analy['flg_eye'] == 1:
                    compsign = '<'
                ewtext = 'EW' + compsign + str(ewval) + r'$\pm$' + ewsigval + r'$\AA$'
                if iline.attrib['EW'].value < 3.* iline.attrib['sig_EW'].value:
                    compsign = '<'
                    ewval = '{:05.2f}'.format(iline.attrib['sig_EW'].value *2)
                    ewtext = 'EW' + compsign + str(ewval) +r'$\AA$'
                ax.text(0.95, 0.1, ewtext, transform=ax.transAxes, ha='right', va='center')
    # Handle boolean switches
    if tight_layout:
        plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
    if show:
        plt.show()
    if return_fig:
        return fig
    else:
        plt.close()

