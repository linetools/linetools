#!/usr/bin/env python

"""
Plot an absorption line with given parameters.
Also print the line data
  Examples:
  lt_absline 1334.3    13.5 7
  lt_absline 1215.6701 14.0 30
"""
import pdb

def plot_absline(wrest,logN,b, show=True):
    """Plot an absorption line with N,b properties

    Parameters
    ----------
    wrest : float
      Rest wavelength (Ang)
    logN : float
      Log10 column
    b : float
      Doppler parameter (km/s)
    show : bool
      Whether to display the plot (set False for running
      tests). Default True.
    """
    import numpy as np
    from linetools.spectralline import AbsLine
    from astropy import units as u

    # Search for the closest absline
    aline = AbsLine(wrest*u.AA, closest=True)

    # Generate a fake wavelength array near the line
    wvoff = 50. # Ang
    dwv = wrest/100000. # Ang (echelle)
    wave = np.arange(wrest-wvoff, wrest+wvoff, dwv)

    # Generate spectrum with voigt
    aline.attrib['N'] = 10**logN * u.cm**-2
    aline.attrib['b'] = b * u.km/u.s
    xspec = aline.generate_voigt(wave=wave*u.AA)
    # get the plotting limits
    ind = np.flatnonzero(xspec.flux.value < 0.9)
    wmin = xspec.wavelength[max(0, ind[1] - 10)]
    wmax = xspec.wavelength[min(len(xspec.flux) - 1,  ind[-2] + 10)]
    #import pdb; pdb.set_trace()
    xspec.constant_sig(0.1) # S/N = 10 per pix

    # Calculate EW
    aline.analy['spec'] = xspec
    aline.analy['wvlim'] = np.array([wrest-15., wrest+15])*u.AA
    aline.measure_ew()
    print(aline)
    print('EW = {:g}'.format(aline.attrib['EW']))

    # Plot
    xspec.plot(xlim=(wmin.to(u.AA).value, wmax.to(u.AA).value), show=show)

def main(args=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Plot an absorption line with the given parameters.')
    parser.add_argument("wrest", type=float, help="Rest wavelength in Angstroms")
    parser.add_argument("logN", type=float, help="log10 column density")
    parser.add_argument("b", type=float, help="b-value in km/s")
     
    pargs = parser.parse_args()
    plot_absline(pargs.wrest, pargs.logN, pargs.b)

