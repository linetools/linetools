"""
#;+ 
#; NAME:
#; spectralline
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for SpectralLine class
#;   23-Jun-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import numpy as np
import pdb
from abc import ABCMeta, abstractmethod
import copy, imp

from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord

from linetools.analysis import utils as lau
from linetools.analysis import absline as laa
from linetools.spectra import io as lsio
from linetools.lists.linelist import LineList

#import xastropy.atomic as xatom
#from xastropy.stats import basic as xsb
#from xastropy.xutils import xdebug as xdb

# class SpectralLine(object):
# class AbsLine(SpectralLine):
# class AbsComponens(AbsLine):

# Class for Spectral line
class SpectralLine(object):
    """Class for a spectral line.  Emission or absorption 

    Attributes:
        ltype: str
          type of line, e.g.  Abs, Emiss
        wrest : Quantity
          Rest wavelength
        z: float
          Redshift
    """
    __metaclass__ = ABCMeta

    # Initialize with wavelength
    def __init__(self, ltype, trans, linelist=None, closest=False,
        z=0.):
        """  Initiator

        Parameters
        ----------
        ltype : string
          Type of Spectral line, 'Abs'
        trans: Quantity or str
          Quantity: Rest wavelength (e.g. 1215.6700*u.AA)
          str: Name of transition (e.g. 'CIV 1548')
        linelist : LineList, optional
          Class of linelist or str setting LineList
        closest : bool, optional
          Take the closest line to input wavelength? [False]
        """

        # Required
        self.ltype = ltype
        if ltype not in ['Abs']:
            raise ValueError('spec/lines: Not ready for type {:s}'.format(ltype))

        # Init
        if not isinstance(trans,(Quantity,basestring)):
            raise ValueError('Rest wavelength must be a Quantity or str')

        # Other
        self.data = {} # Atomic/Moleculare Data (e.g. f-value, A coefficient, Elow)
        self.analy = {'spec': None, # Analysis inputs (e.g. spectrum; from .clm file or AbsID)
            'wvlim': [0., 0.]*u.AA, # Wavelength interval about the line (observed)
            'vlim': [0., 0.]*u.km/u.s, # Velocity limit of line, relative to self.attrib['z']
            'do_analysis': 1 # Analyze
            }
        self.attrib = {   # Properties (e.g. column, EW, centroid)
                       'coord': SkyCoord(ra=0.*u.deg, dec=0.*u.deg),  #  Coords
                       'z': z, 'zsig': 0.,  #  Redshift
                       'v': 0.*u.km/u.s, 'vsig': 0.*u.km/u.s,  #  Velocity relative to z
                       'EW': 0.*u.AA, 'EWsig': 0.*u.AA, 'flgEW': 0 # EW
                       }

        # Fill data
        self.fill_data(trans, linelist=linelist, closest=closest)

    def ismatch(self,inp,Zion=None,RADec=None):
        '''Query whether input line matches on:  z, Z, ion, RA, Dec
        Parameters:
        ----------
        inp: SpectralLine or tuple
          SpectralLine -- Other spectral line for comparison
          tuple -- (z,wrest) float,Quantity
             e.g. (1.3123, 1215.670*u.AA)
        Zion: tuple of ints, optional
          Generally used with tuple input, e.g. (6,2)
        RADec: tuple of Quantities, optional
          Generally used with tuple input e.g. (124.132*u.deg, 29.231*u.deg)

        Returns:
        -------
        answer: bool
          True if a match, else False
        '''
        coord = None
        if isinstance(inp,SpectralLine):
            wrest = inp.wrest
            z = inp.attrib['z']
            if Zion is None:
                Zion = (inp.data['Z'], inp.data['ion'])
            if RADec is None:
                coord = inp.attrib['coord']
        elif isinstance(inp,tuple):
            z = inp[0]
            wrest = inp[1]
        else:
            raise ValueError('ismatch: Bad input')

        # Queries
        answer = ( np.allclose(self.wrest.to(u.AA).value,
                               wrest.to(u.AA).value) &
            np.allclose(self.attrib['z'], z, rtol=1e-6))
        if Zion is not None:
            answer = answer & (self.data['Z'] == Zion[0]) & (self.data['ion'] == Zion[1])
        if (coord is not None) or (RADec is not None):
            if coord is None:
                coord = SkyCoord(ra=RADec[0], dec=RADec[1])
            answer = (answer & (coord.separation(self.attrib['coord']) < 0.1*u.arcsec))

        # Return
        return answer

    def cut_spec(self, normalize=False, relvel=False):
        '''Setup spectrum for analysis.  Splice.  Normalize too (as desired)

        Parameters:
        ----------
        normalize: bool, optional
          Normalize if true (and continuum exists)
        relvel: bool, optional
          Calculate and return relative velocity [False]

        Returns:
        ----------
        fx, sig, dict(wave,velo) -- 
          Arrays (numpy or Quantity) of flux, error, and wavelength/velocity
        '''
        # Checks
        if self.analy['spec'] is None:
            raise ValueError('spectralline.cut_spec: Need to set spectrum!')
        if self.analy['spec'].wcs.unit == 1.:
            raise ValueError('Expecting a unit!')
                    # Velocity

        # Pixels for evaluation
        if np.sum(self.analy['wvlim'].value > 0.):
            pix = self.analy['spec'].pix_minmax(self.analy['wvlim'])[0]
        elif np.sum(np.abs(self.analy['vlim'].value) > 0.):
            pix = self.analy['spec'].pix_minmax(
                self.attrib['z'], self.wrest, self.analy['vlim'])[0]
        else:
            raise ValueError('spectralline.cut_spec: Need to set wvlim or vlim!')
        self.analy['pix'] = pix

        # Cut for analysis
        fx = self.analy['spec'].flux[pix]
        sig = self.analy['spec'].sig[pix]
        wave = self.analy['spec'].dispersion[pix]

        # Velocity array
        self.analy['spec'].velo = self.analy['spec'].relative_vel(
            self.wrest*(1+self.attrib['z']))
        velo = self.analy['spec'].velo[pix] 

        # Normalize?
        if normalize:
            try:
                fx = fx / self.analy['spec'].conti[pix]
            except AttributeError:
                pass
            else:
                sig = sig / self.analy['spec'].conti[pix]

       # Return
        return fx, sig, dict(wave=wave, velo=velo)


    # EW 
    def measure_ew(self, flg=1, initial_guesses=None):
        """EW calculation
        Default is simple boxcar integration
        Observer frame, not rest-frame (use measure_restew below)
          wvlim must be set!
          spec must be set!

        Parameters
        ----------
        flg : int, optional
          1-- Boxcar integration
          2-- Gaussian fit
        
        initial_guesses, optional: tuple of floats
          if a model is chosen (e.g. flg=2, Gaussian) a tuple of (amplitude, mean, stddev)
          can be specified. 

        Fills
        -----
        self.attrib[ 'EW', 'sigEW' ]
          EW and error in observer frame
        """
        imp.reload(lau)
        # Cut spectrum
        fx, sig, xdict = self.cut_spec(normalize=True)
        wv = xdict['wave']

        # Calculate
        if flg == 1: # Boxcar
            EW, sigEW = lau.box_ew( (wv, fx, sig) )
        elif flg == 2: #Gaussian
            EW, sigEW = lau.gaussian_ew( (wv, fx, sig), self.ltype, initial_guesses=initial_guesses)
        else:
            raise ValueError('measure_ew: Not ready for this flag {:d}'.format(flg))

        # Fill
        self.attrib['EW'] = EW 
        self.attrib['sigEW'] = sigEW 

    # EW 
    def measure_restew(self,**kwargs):
        """  Rest EW calculation
        Return rest-frame.  See "measure_ew" above for details
        """
        # Standard call
        self.measure_ew(**kwargs)

        # Push to rest-frame
        self.attrib['EW'] = self.attrib['EW'] / (self.attrib['z']+1)
        self.attrib['sigEW'] = self.attrib['sigEW'] / (self.attrib['z']+1)

    # Output
    def __repr__(self):
        txt = '[{:s}:'.format(self.__class__.__name__)
        try:
            txt = txt+' {:s},'.format(self.data['name'])
        except KeyError:
            pass
        txt = txt + ' wrest={:g}'.format(self.wrest)
        txt = txt + ']'
        return (txt)

# ###########################################
# Class for Generic Absorption Line System
class AbsLine(SpectralLine):
    """Spectral absorption line
    trans: Quantity or str
      Quantity: Rest wavelength (e.g. 1215.6700*u.AA)
      str: Name of transition (e.g. 'CIV 1548')
        [Note: for an unknown transition use string 'unknown']
    """
    # Initialize with a .dat file
    def __init__(self, trans, **kwargs):
        # Generate with type
        SpectralLine.__init__(self,'Abs', trans, **kwargs)

    def print_specline_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'AbsLine'

    def fill_data(self,trans, linelist=None, closest=False):
        ''' Fill atomic data and setup analy
        Parameters:
        -----------
        trans: Quantity or str
          Quantity: Rest wavelength (e.g. 1215.6700*u.AA)
          str: Name of transition (e.g. 'CIV 1548')
            [Note: for an unknown transition use string 'unknown']
        linelist : LineList, optional
          Class of linelist or str setting LineList
        closest : bool, optional
          Take the closest line to input wavelength? [False]
        '''

        # Deal with LineList
        if linelist is None:
            llist = LineList('ISM')
        elif isinstance(linelist,basestring):
            llist = LineList(linelist)
        elif isinstance(linelist,LineList):
            llist = linelist
        else:
            raise ValueError('Bad input for linelist')

        # Closest?
        llist.closest = closest

        # Data
        newline = llist[trans]
        self.data.update(newline)

        # Update
        self.wrest = self.data['wrest']
        self.trans = self.data['name']

        #
        self.analy.update( {
            'flg_eye': 0,
            'flg_limit': 0, # No limit
            'datafile': '', 
            'name': self.data['name']
            })

        # Additional attributes for Absorption Line
        self.attrib.update({'N': 0., 'Nsig': 0., 'flagN': 0, # Column
                       'b': 0.*u.km/u.s, 'bsig': 0.*u.km/u.s  # Doppler
                       } )
    # Voigt
    def generate_voigt(self, wave=None, **kwargs):
        """  Generate a Voigt profile model for the absorption line
        in a given spectrum.

        Parameters:
        ----------
        wave: Quantity array
          Wavelength array on which to calculate the line
          Must be set if self.analy['spec'] is not filled

        Returns:
        ----------
        spec: XSpectrum1D
          Spectrum with the input wavelength and the absorbed flux
        """
        from linetools.analysis import voigt as lav
        reload(lav)
        # Checks
        if self.attrib['N'] < 1.:
            raise ValueError("Need to initialize log column density in attrib['N']")
        if self.attrib['b'] < 1.*u.km/u.s:
            raise ValueError("Need to initialize Doppler parameter in attrib['b']")
        if wave is None:
            # Assume a spectrum has been loaded already
            try:
                wave = self.analy['spec'].dispersion
            except:
                raise ('You must provide a wavelength array in generate_voigt')

        # Main call
        spec = lav.voigt_from_abslines(wave, self, **kwargs)
        return spec

    # AODM
    def measure_aodm(self, nsig=3.):
        """  AODM calculation
        Parameters
        ----------
        nsig: float, optional
          Number of sigma significance required for a "detection"

        Fills:
        -------
        self.attrib[ 'N', 'sigN', 'logN', 'sig_logN' ]  
          Column densities and errors, linear and log
        """
        # Cut spectrum
        fx, sig, xdict = self.cut_spec(normalize=True)
        velo = xdict['velo']

        # Calculate
        N,sigN,flg_sat = laa.aodm( (velo, fx, sig), (self.wrest,self.data['f']) )

        # Flag
        if flg_sat:
            self.attrib['flagN'] = 2
        else:
            if N > nsig*sigN:
                self.attrib['flagN'] = 1
            else:
                self.attrib['flagN'] = 3

        # Values
        self.attrib['N'] = N
        self.attrib['sigN'] = sigN

        # Log
        logN, sig_logN = laa.log_clm(self.attrib)
        self.attrib['logN'] = logN # Dimensionless
        self.attrib['sig_logN'] = sig_logN # Dimensionless


    # Output
    def __repr__(self):
        txt = '[{:s}:'.format(self.__class__.__name__)
        # Name
        try:
            txt = txt+' {:s},'.format(self.data['name'])
        except KeyError:
            pass
        # wrest
        txt = txt + ' wrest={:.4f}'.format(self.wrest)
        # fval
        try:
            txt = txt+', f={:g}'.format(self.data['fval'])
        except KeyError:
            pass
        txt = txt + ']'
        return (txt)

def many_abslines(all_wrest, llist):
    '''Generate a list of AbsLine objects
    Useful for when you have many lines (>1000) to generate
    that have similar wrest.  Uses deepcopy

    Parameters:
    -----------
    all_wrest: list of lines
    llist: LineList

    Returns:
    ----------
    abs_lines : list 
      of AbsLine Objects
    '''
    # Find unique lines
    wrestv =  np.array([iwrest.value for iwrest in all_wrest]) 
    uniq_wrest = np.unique( wrestv )

    # Generate a simple dict
    adict = {}
    unit = all_wrest[0].unit
    for iuni in uniq_wrest:
        adict[iuni] = AbsLine(iuni*unit,linelist=llist)

    # Copy em up
    abs_lines = []
    for iwrestv in wrestv:
        abs_lines.append(copy.deepcopy(adict[iwrestv]))

    # Return
    return abs_lines
