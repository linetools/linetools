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
from abc import ABCMeta, abstractmethod

from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity

from linetools.spectra import io as lsio
from linetools.lists.linelist import LineList

#import xastropy.atomic as xatom
#from xastropy.stats import basic as xsb
#from xastropy.xutils import xdebug as xdb

# class SpectralLine(object):
# class AbsLine(SpectralLine):

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
            'wvlim': [0., 0.], # Wavelength interval about the line (observed)
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
        if np.sum(self.analy['wvlim']) == 0.:
            raise ValueError('spectralline.cut_spec: Need to set wvlim!') # Could use VMNX
        if self.analy['spec'] is None:
            raise ValueError('spectralline.cut_spec: Need to set spectrum!')
        if self.analy['spec'].wcs.unit == 1.:
            raise ValueError('Expecting a unit!')

        # Pixels for evaluation
        pix = self.analy['spec'].pix_minmax(self.analy['wvlim'])[0]

        # Cut for analysis
        fx = self.analy['spec'].flux[pix]
        sig = self.analy['spec'].sig[pix]
        wave = self.analy['spec'].dispersion[pix]


        # Normalize
        if normalize:
            try:
                fx = fx / self.analy['spec'].conti[pix]
            except AttributeError:
                pass
            else:
                sig = sig / self.analy['spec'].conti[pix]

        # Velocity
        self.analy['spec'].velo = self.analy['spec'].relative_vel(
            self.wrest*(1+self.attrib['z']))
        # Cut
        velo = self.analy['spec'].velo[pix] 
        # Return
        return fx, sig, dict(wave=wave, velo=velo)


    # EW 
    def box_ew(self):
        """  EW calculation
        Default is simple boxcar integration
        Observer frame, not rest-frame
          wvlim must be set!
          spec must be set!

        Parameters
        ----------

        Returns:
          EW, sigEW : EW and error in observer frame
        """
        # Cut spectrum
        fx, sig, xdict = self.cut_spec(normalize=True)
        wv = xdict['wave']

        # dwv
        dwv = wv - np.roll(wv,1)
        dwv[0] = dwv[1]


        # Simple boxcar
        EW = np.sum( dwv * (1. - fx) ) 
        varEW = np.sum( dwv**2 * sig**2 )
        sigEW = np.sqrt(varEW) 


        # Fill
        self.attrib['EW'] = EW 
        self.attrib['sigEW'] = sigEW 

        # Return
        return EW, sigEW
            
    # EW 
    def restew(self):
        """  Rest EW calculation
        Return rest-frame.  See "box_ew" above for details
        """
        # Standard call
        EW,sigEW = self.box_ew()
        # Push to rest-frame
        self.attrib['EW'] = EW / (self.attrib['z']+1)
        self.attrib['sigEW'] = sigEW / (self.attrib['z']+1)

        # Return
        return self.attrib['EW'], self.attrib['sigEW'] 

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
                       'b': 0., 'bsig': 0.  # Doppler
                       } )

    '''
    # Perform AODM on the line
    def aodm(self, flg_sat=None):
        """  AODM calculation

        Parameters
        ----------

        Returns:
          N, sigN : Column and error in linear space
          flg_sat:  Set to True if saturated pixels exist
        """
        # Cut spectrum
        fx, sig, velo = self.cut_spec(normalize=True, relvel=True)

        # dv
        delv = velo - np.roll(velo,1)
        delv[0] = delv[1]

        # Atomic data
        cst = (const.m_e.cgs*const.c.cgs / (np.pi * 
            (const.e.esu**2).cgs)).to(u.AA*u.s/(u.km*u.cm**2))
        cst = cst/(self.data['f']*self.wrest) #/ (u.km/u.s) / u.cm * (u.AA/u.cm)

        # Mask
        mask = (fx == fx) # True = good
        nndt = Quantity(np.zeros(len(fx)), unit='s/(km cm cm)')

        # Saturated?
        satp = np.where( (fx <= sig/5.) | (fx < 0.05) )[0]
        if len(satp) > 0:
            mask[satp] = False
            lim = np.where(sig[satp] > 0.)[0]
            if len(lim) > 0:
                sub = np.maximum(0.05, sig[satp[lim]]/5.)
                nndt[satp[lim]] = np.log(1./sub)*cst
                flg_sat = True
        # AODM
        nndt[mask] = np.log(1./fx[mask])*cst

        # Sum it
        ntot = np.sum( nndt*delv )
        tvar = np.sum( (delv*cst*sig/fx)**2 )

        # Fill
        self.attrib['N'] = ntot
        self.attrib['sigN'] = np.sqrt(tvar)

        # Log
        logN = np.log10( self.attrib['N'].value ) 
        lgvar = ((1. / (np.log(10.)*self.attrib['N'].value))**2) * self.attrib['sigN'].value**2
        sig_logN = np.sqrt(lgvar)
        self.attrib['logN'] = logN # Dimensionless
        self.attrib['sig_logN'] = sig_logN # Dimensionless

        # Return
        return self.attrib['N'], self.attrib['sigN']
    '''


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

