""" Utlities for the analysis of emission lines
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import warnings

from astropy import units as u
from astropy import constants as const
from astropy.io import ascii
from astropy.utils import isiterable
from linetools.lists.linelist import LineList

# Atomic constant
atom_cst = (const.m_e.cgs*const.c.cgs / (np.pi * 
    (const.e.esu**2).cgs)).to(u.AA*u.s/(u.km*u.cm**2))

# e2/(me*c) in CGS
e2_me_c_cgs = (const.e.esu**2 / (const.c.to('cm/s') * const.m_e.to('g')))


# Estimate metallicity from a tried-and-true method
def metallicity(method, emsystem):
    """ Calculate a nebular abundance from an input set of Emission lines

    Parameters
    ----------
    method : str
      Name of the method to be applied
      PG16 -- http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1601.08217
          Requires Hbeta, [OII], [OIII], [NII], [SII]
          Return r_val and s_val
    emsystem : EmSystem

    Returns
    -------
    """
    if method == 'PG16':
        # Requires Hbeta, [OII], [OIII], [NII], [SII]
        R2 = (emsystem.get_emline('[OII] 3726').attrib['flux'] +
              emsystem.get_emline('[OII] 3729').attrib['flux']) / emsystem.get_emline('Hbeta').attrib['flux']
        R3 = (emsystem.get_emline('[OIII] 4959').attrib['flux'] +
              emsystem.get_emline('[OIII] 5007').attrib['flux']) / emsystem.get_emline('Hbeta').attrib['flux']
        N2 = (emsystem.get_emline('[NII] 6548').attrib['flux'] +
              emsystem.get_emline('[NII] 6584').attrib['flux']) / emsystem.get_emline('Hbeta').attrib['flux']
        S2 = (emsystem.get_emline('[SII] 6716').attrib['flux'] +
              emsystem.get_emline('[SII] 6731').attrib['flux']) / emsystem.get_emline('Hbeta').attrib['flux']
        # Proceed
        if np.log10(N2) < -0.6:
            r_val = 7.932 + 0.944*np.log10(R3/R2) + 0.695*np.log10(N2) + \
                ((0.97 - 0.291*np.log10(R3/R2)) - 0.019*np.log10(N2))*np.log10(R2)

            s_val = 8.072 + 0.789*np.log10(R3/S2) + 0.726*np.log10(N2) + \
                (1.069 - 0.170*np.log10(R3/S2) +0.022*np.log10(N2))*np.log10(S2)
        else:
            r_val = 8.589 + 0.022*np.log10(R3/R2) + 0.399*np.log10(N2) + \
                (-0.137 + 0.164*np.log10(R3/R2) + 0.589*np.log10(N2))*np.log10(R2)

            s_val = 8.424 + 0.030*np.log10(R3/S2) + 0.751*np.log10(N2) + \
                (-0.349 + 0.182*np.log10(R3/S2) +0.508*np.log10(N2))*np.log10(S2)
        return r_val.decompose().value, s_val.decompose().value

