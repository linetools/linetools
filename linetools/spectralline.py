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

import numpy as np
import pdb
from abc import ABCMeta, abstractmethod
import copy

from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity

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
                       'RA': 0.*u.deg, 'Dec': 0.*u.deg,  #  Coords
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
        if isinstance(inp,SpectralLine):
            wrest = inp.wrest
            z = inp.attrib['z']
            if Zion is None:
                Zion = (inp.data['Z'], inp.data['ion'])
            if RADec is None:
                RADec = (inp.attrib['RA'], inp.attrib['Dec'])
        elif isinstance(inp,tuple):
            z = inp[0]
            wrest = inp[1]
        else:
            raise ValueError('ismatch: Bad input')

        # Queries
        answer = ( np.allclose(self.wrest, wrest) &
            np.allclose(self.attrib['z'], z, rtol=1e-6))
        if Zion is not None:
            answer = answer & (self.data['Z'] == Zion[0]) & (self.data['ion'] == Zion[1])
        if RADec is not None:
            answer = (answer & np.allclose(self.attrib['RA'], RADec[0]) &
                np.allclose(self.attrib['Dec'], RADec[1]) )

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
        """  EW calculation
        Default is simple boxcar integration
        Observer frame, not rest-frame (use measure_restew below)
          wvlim must be set!
          spec must be set!

        Parameters
        ----------
        flg: int, optional
          1: Boxcar integration
          2: Gaussian fit
        
        initial_guesses, optional: tuple of floats
          if a model is chosen (e.g. flg=2, Gaussian) a tuple of (amplitude, mean, stddev)
          can be specified. 

        Fills:
        -------
        self.attrib[ 'EW', 'sigEW' ] : 
          EW and error in observer frame
        """
        reload(lau)
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
            self.llist = LineList('ISM')
        elif isinstance(linelist,basestring):
            self.llist = LineList(linelist)
        elif isinstance(linelist,LineList):
            self.llist = linelist
        else:
            raise ValueError('Bad input for linelist')

        # Closest?
        self.llist.closest = closest

        # Data
        newline = self.llist[trans]
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
        reload(laa)

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
        logN = np.log10( self.attrib['N'].value ) 
        lgvar = ((1. / (np.log(10.)*self.attrib['N'].value))**2) * self.attrib['sigN'].value**2
        sig_logN = np.sqrt(lgvar)
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


 ###########################################
# Class for Generic Component of Absorption Line Systems
class AbsComponent(object):
    
    """A component corresponds to list of AbsLines from a common ion species 
    (e.g. HI Lyman series, CIV doublet, etc), which is supposed to characterize a 
    single observational entity; therefore it has a single unique position 
    (RA, DEC, redshift), column density, velocity range, etc. 

    [under construction]

    """

    def __init__(self, abslines):
        """  
        Parameters
        ----------
        abslines : a list of AbsLine objects or a single AbsLine object
    
        """
        #Unify format and always return a list of AbsLine objects
        abslines = self.unify_format(abslines)
        
        #check whether individual AbsLines are from the same position
        #check whether individual AbsLines are coming from the same ion species
        for absline in abslines[1:]:
            self.check_same_position(abslines[0],absline)
            self.check_same_ion_species(abslines[0],absline)
            
        #check for duplicates, if found raise an error
        for i in range(len(abslines)-1):
            for j in range(i+1,len(abslines)):
                self.check_no_same_transition(abslines[i],abslines[j])

        #all checks passed, then we are good to initialize self.abslines
        self.abslines = abslines

        # Additional attributes for Absorption Component
        # [to be implemented]

    def unify_format(self,abslines):
        """If this is a single AbsLine it returns it as a list, 
        otherwise keep it as a list"""
        if isinstance(abslines,AbsLine):
            abslines = [abslines]
        elif isinstance(abslines,list):
            abslines = abslines
        else:
            raise ValueError('AbsComponent: abslines is not a list of AbsLine objects nor a single AbsLine object!')
        return abslines

    def check_same_ion_species(self,absline_1,absline_2):
        """For 2 given AbsLine objects, this function checks whether these 
        are coming same ion species, otherwise it raises an error.

        Input
        -----
        absline_1 : AbsLine object
        absline_2 : AbsLine object

        """
        if isinstance(absline_1,AbsLine):
            pass
        else:
            raise ValueError('absline_1 is not an AbsLine object')
        if isinstance(absline_2,AbsLine):
            pass
        else:
            raise ValueError('absline_2 is not an AbsLine object')
        
        data_1 = absline_1.data
        data_2 = absline_2.data
        
        if (data_1['Z'] != data_2['Z']):
            raise ValueError('AbsLine {} and {} are not from the same atomic number!'.format(data_1['name'],data_2['name']))
        if (data_1['ion'] != data_2['ion']):
            raise ValueError('Absline {} and {} are not from the same ionization state!'.format(data_1['name'],data_2['name']))
        if (data_1['Ej'] != data_2['Ej']):
            raise ValueError('Absline {} and {} are not from the same lower energy level!'.format(data_1['name'],data_2['name']))
        
    def check_same_position(self,absline_1,absline_2,ra_tol=None,dec_tol=None,z_tol=None):
        """Checks whether 2 given absorption lines are coming from the
        same (RA,DEC,z) position, otherwise it raises an error
        
        Input
        -----
        absline_1 : AbsLine object
        absline_2 : AbsLine object

        [Note: we need to decide what to do about tolerances] 

        """
        #check format
        if isinstance(absline_1,AbsLine):
            pass
        else:
            raise ValueError('absline_1 is not an AbsLine object')
        if isinstance(absline_2,AbsLine):
            pass
        else:
            raise ValueError('absline_2 is not an AbsLine object')
        
        #check RA
        ra1 = absline_1.attrib['RA']
        ra2 = absline_2.attrib['RA']
        if np.allclose(ra1,ra2):
            pass
        else:
            raise ValueError('AbsLine {} and {} do not have the same RA coordinates!'.format(data_1['name'],data_2['name']))
        
        #check Dec.
        dec1 = absline_1.attrib['Dec']
        dec2 = absline_2.attrib['Dec']
        if np.allclose(dec1,dec2):
            pass
        else:
            raise ValueError('AbsLine {} and {} do not have the same Dec. coordinates!'.format(data_1['name'],data_2['name']))
        
        #Check redshift
        z1 = absline_1.attrib['z']
        z2 = absline_2.attrib['z']
        if np.allclose(z1,z2):
            pass
        else:
            raise ValueError('AbsLine {} and {} do not have the same Dec. coordinates!'.format(data_1['name'],data_2['name']))

    def check_no_same_transition(self,absline_1,absline_2):
        """For 2 given AbsLine objects, this function checks whether these 
        are not the same transition, otherwise it raises an error.

        Input
        -----
        absline_1 : AbsLine object
        absline_2 : AbsLine object

        """
        if isinstance(absline_1,AbsLine):
            pass
        else:
            raise ValueError('absline_1 is not an AbsLine object')
        if isinstance(absline_2,AbsLine):
            pass
        else:
            raise ValueError('absline_2 is not an AbsLine object')
        
        data_1 = absline_1.data
        data_2 = absline_2.data
        
        if (data_1['name'] == data_2['name']):
            raise ValueError('Two AbsLine are the exact same transition: {}!'.format(data_1['name']))    

    def add_abslines(self,abslines):
        """Adds an AbsLine object to the current AbsComponent list.
        The new absline has to be of the same ion species as the current 
        ones, otherwise an error is raised

        Parameters
        ----------
        absline : a single AbsLine object or a list of AbsLine objects
    
        """
        #Unify format and always return a list of AbsLine objects
        abslines = self.unify_format(abslines)
        
        # check all elements are of the same position and species as the original one
        for absline in abslines:
            self.check_same_position(self.abslines[0],absline)
            self.check_same_ion_species(self.abslines[0],absline)

        #check for duplicates within the new abslines, if found raise an error
        for i in range(len(abslines)-1):
            for j in range(i+1,len(abslines)):
                self.check_no_same_transition(abslines[i],abslines[j])
        
        #check for duplicates between the new and original abslines
        for i in range(len(abslines)):
            for j in range(len(self.abslines)):
                self.check_no_same_transition(abslines[i],self.abslines[j])
        
        # all test passed, append the new abslines to self.abslines
        for absline in abslines:
            self.abslines += [absline]

    def remove_absline(self,absline_name):
        """Removes absline by name from the current self.abslines list

        Parameters
        ----------
        absline_name : a single string with the name of transition to remove

        """

        #check format
        if isinstance(absline_name,str):
            pass
        else:
            raise ValueError('Input is not a string!')
        
        for i in range(len(self.abslines)):
            name = self.abslines[i].data['name']
            if name == absline_name:
                #remove such instance
                del self.abslines[i]
                break

    def remove_abslines(self,abslines_names):
        """Removes multiple AbsLines in the Component given a list of 
        names

        Parameters
        ----------
        abslines_names : a list of strings or a single string with the name of transition to remove

        """
        #check format
        if isinstance(abslines_names,str):
            abslines_names = [abslines_names] #if single string make it a list
        elif isinstance(abslines_names,list):
            if all(isinstance(absline_name,str) for absline_name in abslines_names):
                pass
            else:    
                raise ValueError('An element of abslines_names is not a string!')

        for absline_name in abslines_names:
            self.remove_absline(absline_name)


