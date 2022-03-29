"""
  Module for automatically finding continuum
"""

import numpy as np
from matplotlib import pyplot as plt
import linetools.utils as ltu

from linetools.spectra.io import readspec
import linetools.spectra.xspectrum1d as xspec
import pdb


def contknots(spec, sm=10, npix=40, lchmin=20, ewsnlim=3, lchmax=None,
              showcont=True, lam0=None, dlam=None, yfmax=None, ymin=None,
              outf=None, overwrite=False):
    """ Finds continuum knots, and pixels that are located on continuum .
        Identifies pixels that correspond to absorption lines, and do not includes them in continuum.

        Parameters sm, npix, lchmin, ewsnlim that best define continuum might be different for different data sets.
        Recomended using lt_continuumfit script (gui) after.

    Parameters
    ----------
    spec : XSpectrum1D or str
      input spectrum, or a file that contains input XSpectrum1D spectrum
    sm : float, optional
      number of pixels how much to smooth input spectrum when defining continuum
    npix : int, optional
      number of pixels, such that every maximum flux value is greater than for npix pixels with lower and at
      higher wavelengths
    lchmin : int, optional
      minimum number of pixels in the intervals that define continuum
    ewsnlim : float, optional
      minimum value for S/N ratio of a quantity similar to equivalent width
    lchmax : int, optional
      ~maximum number of pixels in the intervals that define continuum
    showcont : bool, optional
      show continuum? if yes:
      orange line with circles - continuum and knots
      green - pixels used to define continuum knots
      black - original spectrum , gray - error in flux
      blue line - initially defined continuum (not used),
                  and initially defined continuum shifted by the error in flux in the smoothed spectrum (also not used)
      red line - smoothed spectrum
    lam0 : float, optional
      lowest wavelength to plot
    dlam : float, optional
      wavelenght interval width to plot; highest wavelength shown will be lam0+dlam
    yfmax : float, optional
      maximum ylim is yfmax * np.median(flux)
    ymin : float, optional
      minimum ylim is yfmin * np.median(flux)
    outf : str, optional
      output file with knots (.jsn or .json type)
    overwrite : bool, optional
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

        iflx3max = np.max(flx3[i1:i2+1])
        if iflx3max == flx3[i]:
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

    # find pixels that correspond to the continuum (= all pixels except the ones that correspond to absorption lines)

    # use only pixels where cont > 0
    nz = np.where(cont0 > 0)[0]
    flxcp = np.asarray(flx[nz]) # we will modify this
    flxcp1 = np.asarray(flx[nz])
    cont0cp = np.asarray(cont0[nz])
    sigfcp = np.asarray(sigf[nz])
    wavcp = wav[nz]

    # calculate iews and isigews
    k = np.where(flxcp > cont0cp)[0] # absorption lines will be between these pixels
    flxcp[k]=cont0cp[k]
    iews = 1 - np.asarray(flxcp) / np.asarray(cont0cp)
    isigew2s = (np.asarray(sigfcp) / np.asarray(cont0cp))**2

    # identify absorption lines, and find index (in flxcp) of the bluest pixel of each line
    cumiews = np.cumsum(iews)
    cumisigew2s = np.cumsum(isigew2s)
    kcumiews = cumiews[k]
    kcumisigew2s = cumisigew2s[k]
    dkcumiews = np.diff(kcumiews)
    dkcumisigew2s = np.diff(kcumisigew2s)
    jj = np.where(dkcumiews > ewsnlim * (dkcumisigew2s**0.5))[0]
    knew = k[jj] # index (in flxcp) of the bluest pixel of absorption line

    # identify pixels on continuum
    arr = np.repeat(0,len(flxcp))
    arr[knew] = 1
    arr[k[np.asarray(jj) + 1]] = arr[k[np.asarray(jj) + 1]] -1
    cumarr = np.cumsum(arr)
    # ind0s = np.where(cumarr == 1)[0] # absorption lines pixels
    allind = np.where((cumarr != 1) & (flxcp != 0))[0] # continuum pixels


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
    chnks.append(ich)

    if lchmax is None:
        lchmax = 5 * lchmin


    # 4) Define continuum knots
    #    take into account only intervals with > lchmin elements
    cwav = []
    cflx = []
    knots = []
    knotpixs = []
    for ichnk in chnks:
        if len(ichnk) >= lchmin:
            for j in ichnk:
                knotpixs.append(j)
            #
            if len(ichnk) < lchmax:
                cwav.append(np.mean(wavcp[ichnk].value))
                cflx.append(np.mean(flxcp1[ichnk])) #.value))
                knots.append([np.mean(wavcp[ichnk].value),np.mean(flxcp1[ichnk])]) #  .value)])
                #for j in ichnk:
                #    knotpixs.append(j)
            else:
                # break into smaller intervals first
                nch = round(len(ichnk)/lchmax) # + 1
                for j in range(nch):
                    cwav.append(np.mean(wavcp[ichnk[j*lchmax:(j+1)*lchmax]].value))
                    cflx.append(np.mean(flxcp1[ichnk[j*lchmax:(j+1)*lchmax]]))  # .value))
                    knots.append([np.mean(wavcp[ichnk[j*lchmax:(j+1)*lchmax]].value), np.mean(flxcp1[ichnk[j*lchmax:(j+1)*lchmax]])])  # .value)])
                #
                cwav.append(np.mean(wavcp[ichnk[nch * lchmax:]].value))
                cflx.append(np.mean(flxcp1[ichnk[nch * lchmax:]]))  # .value))
                knots.append([np.mean(wavcp[ichnk[nch * lchmax:]].value),
                              np.mean(flxcp1[ichnk[nch * lchmax:]])])  # .value)])



    # input parameters lam0, dlam, yfmax and ymin
    if lam0 is None:
        lam0 = np.min(wav.value)
    if dlam is None:
        dlam = np.max(wav.value) - lam0
    #
    ii = np.where((wav.value > lam0) & (wav.value < lam0 + dlam))[0]
    if yfmax is None:
        yfmax = np.max(flx[ii])/np.median(flx)
    if ymin is None:
        ymin = np.min(flx[ii])/np.median(flx)



    # plot spectrum without these lines
    if showcont:
        plt.figure(figsize=(17, 4))
        plt.plot(wav, sigf, color='gray')
        plt.plot(wav, flx, color='black')
        plt.plot(wav3, flx3, color='red')
        plt.plot(wav3[ipixs], flx3[ipixs], color='blue')

        plt.plot(wavcp[allind], flxcp1[allind], color='limegreen')
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






