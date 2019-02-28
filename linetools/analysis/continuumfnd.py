"""
  Module for automatically finding continuum
"""

import numpy as np
from matplotlib import pyplot as plt
import linetools.utils as ltu

from linetools.spectra.io import readspec
import linetools.spectra.xspectrum1d as xspec
import pdb


def contknots(spec, sm=10, npix=40, lchmin=20, ewsnlim=2,
              showcont=True, lam0=1500, dlam=20, yfmax=2, ymin=0,
              outf=None, overwrite=False):
    """ Finds continuum knots, and pixels that are located on continuum .
        Identifies pixels that correspond to absorption lines, and do not includes them in continuum.

        Parameters sm, npix, lchmin, ewsnlim that best define continuum might be different for different data sets.
        Recomended using lt_continuumfit script (gui) after.

    Parameters
    ----------
    spec : XSpectrum1D or str
      input spectrum, or a file that contains input XSpectrum1D spectrum
    sm : float
      number of pixels how much to smooth input spectrum when defining continuum
    npix : int
      number of pixels, such that every maximum flux value is greater than for npix pixels with lower and at
      higher wavelengths
    lchmin : int
      minimum number of pixels in the intervals that define continuum
    ewsnlim : float
      minimum value for S/N ratio of a quantity similar to equivalent width
    showcont : bool
      show continuum? if yes:
      orange line with circles - continuum and knots
      green - pixels used to define continuum knots
      black - original spectrum , gray - error in flux
      blue line - initially defined continuum (not used),
                  and initially defined continuum shifted by the error in flux in the smoothed spectrum (also not used)
      red line - smoothed spectrum
    lam0 : float
      lowest wavelength to plot
    dlam : float
      wavelenght interval width to plot; highest wavelength shown will be lam0+dlam
    yfmax : float
      maximum ylim is yfmax * np.median(flux)
    ymin : float
      minimum ylim is yfmin * np.median(flux)
    outf : str
      output file with knots (.jsn or .json type)
    overwrite : bool
      overwrite?

    Returns:
    --------
    knots : list of 2-element lists
      List of knots (wavelength, flux) that define continuum.
    knotpixs : list of int
      Indices of pixels that define continuum

    """

    # input spectrum
    if type(spec) is xspec.XSpectrum1D:
        sp = spec
    else:
        try:
            sp = readspec(spec)
        except OSError:
            print('Not valid input. Please provide XSpectrum1D object, or a file that contains XSpectrum1D object')
            return

    wav = sp.wavelength
    flx = sp.flux
    sigf = sp.sig

    # smoothed input spectrum by smsm pixels
    sp3 = sp.box_smooth(sm)
    wav3 = sp3.wavelength
    flx3 = sp3.flux
    sigf3 = sp3.sig


    # 1) Define pixels that are on the continuum.
    #    Smooth the spectrum, find pixels around maximums in the flux in the smoothed spectrum,
    #    assume that the continuum is well defined by these pixels only.
    #    This flux will be too high, and is used only in the step (1)
    ipixs = []  # continuum pixels
    for i in range(len(flx3)):
        i1 = np.max([0, i - npix])
        i2 = np.min([i + npix, (len(flx3)) - 1])
        ii = np.linspace(i1, i2, i2 - i1 + 1)
        isum = 0
        for j in ii:
            if (flx3[int(j)] <= flx3[i]):
                isum = isum + 1
            if (flx3[int(j)] >= flx3[i]):
                isum = isum - 1
        if (isum) == len(ii) - 1:
            ipixs.append(i)

    # interpolate continuum to all pixels
    cont0ipixs = flx3[ipixs]
    cont0 = np.interp(wav, wav3[ipixs], cont0ipixs)


    # 2) Identify pixels that correspond to absorption lines. Identify all other pixels.
    #    Define continuum as in (1),
    #    define absorption line as a set of adjacent pixels that are below the continuum level; flux of the boundary
    #      pixels exceeds continuum level
    #    Define EW = sum of (1 - flux/continuum) of these adjacent pixels
    #    and a condition that EW > ewsnlim * sigma(EW) ; default ewsnlim = 2

    # find pixels that correspond to absorption lines
    ind0s = [] # pixels corresponding to absorption lines
    imin = -1
    iew = 0
    isigew = 0

    for i in range(len(flx)):
        if flx[i] == 0:
            ind0s.append(i)
        else:
            if flx[i] < cont0[i]:
                iew = iew + (1 - flx[i] / cont0[i])
                isigew = (isigew ** 2 + (sigf[i] / cont0[i]) ** 2) ** 0.5
                if imin == -1:
                    imin = np.max([i - 1, 0])
            else:
                if (iew > ewsnlim * isigew):
                    for j in range(i - imin + 1):
                        ind0s.append(j + imin)
                iew = 0
                isigew = 0
                imin = -1

    # Find pixels that do not correspond to absorption lines
    ind0s = list(set(ind0s))
    allind = list(range(len(flx)))
    for i in ind0s:
        allind.remove(i) # allind - new continuum pixels


    # 3) Group continuum pixels from (2) into intervals of adjacent pixels
    chnks = []  # list of all intervals
    ich = []  # i-th interval
    for i in range(len(allind)):
        if len(ich) == 0:
            ich.append(allind[i])
        else:
            if (allind[i] - allind[i - 1]) == 1:
                ich.append(allind[i])
            else:
                chnks.append(ich)
                ich = [allind[i]]


    # 4) Define continuum knots
    #    take into account only intervals with > lchmin elements
    cwav = []
    cflx = []
    knots = []
    knotpixs = []
    for ichnk in chnks:
        if len(ichnk) >= lchmin:
            cwav.append(np.mean(wav[ichnk].value))
            cflx.append(np.mean(flx[ichnk].value))
            knots.append([np.mean(wav[ichnk].value),np.mean(flx[ichnk].value)])
            for j in ichnk:
                knotpixs.append(j)



    # plot spectrum without these lines
    if showcont:
        plt.figure(figsize=(17, 4))
        plt.plot(wav, sigf, color='gray')
        plt.plot(wav, flx, color='black')
        plt.plot(wav3, flx3, color='red')
        plt.plot(wav3[ipixs], flx3[ipixs], color='blue')

        plt.plot(wav[allind], flx[allind], color='limegreen')
        plt.plot(cwav, cflx, color='orange')

        plt.ylim(np.median(flx) * ymin, np.median(flx) * yfmax)
        plt.xlim(lam0, lam0 + dlam)
        plt.show()


    # write knots in a .json file (if defined)
    if outf is not None:
        knots = sorted(tuple(float(val) for val in cp) for
               cp in knots)
        ltu.savejson(outf, knots, overwrite=overwrite)


    return knots, knotpixs






