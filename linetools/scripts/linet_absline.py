#!/usr/bin/env python

"""
Plot an absorption line with given parameters.
Also print the line data
  Examples:
  linet_absline 1334.3    13.5 7
  linet_absline 1215.6701 14.0 30
"""

import pdb


def plot_absline(wrest,logN,b):
    """Plot an absorption line with N,b properties
    Parmaeters
    ----------
    wrest : float
      Rest wavelength (Ang)
    logN : float
      Log10 column
    b : float
      Doppler parameter (km/s)
    """
    import numpy as np
    from linetools.spectra.xspectrum1d import XSpectrum1D
    from linetools.lists.linelist import LineList
    from linetools.spectralline import AbsLine
    from linetools.analysis import voigt as lav
    from astropy import units as u

    # Search for the closest absline
    aline = AbsLine(wrest*u.AA, closest=True)

    # Generate a fake wavelength array near the line
    wvoff = 50. # Ang
    dwv = wrest/100000. # Ang (echelle)
    wave = np.arange(wrest-wvoff, wrest+wvoff, dwv)

    # Generate spectrum with voigt
    aline.attrib['N'] = logN
    aline.attrib['b'] = b * u.km/u.s
    xspec = aline.generate_voigt(wave=wave*u.AA)
    xspec.constant_sig(0.1) # S/N = 10 per pix


    # Calculate EW
    aline.analy['spec'] = xspec
    aline.analy['wvlim'] = np.array([wrest-15., wrest+15])*u.AA
    aline.measure_ew()
    print(aline)
    print('EW = {:g}'.format(aline.attrib['EW']))

    # Plot
    xspec.plot()

def main(args=None):
    from astropy.utils.compat import argparse
    # Parse
    parser = argparse.ArgumentParser(description='Parser for linet_absline')
    parser.add_argument("wrest", type=float, help="Rest wavelength in Angstroms")
    parser.add_argument("logN", type=float, help="log10 column density")
    parser.add_argument("b", type=float, help="b-value in km/s")
     
    pargs = parser.parse_args()
    plot_absline(pargs.wrest, pargs.logN, pargs.b)

