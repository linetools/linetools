""" Classes for an emission or absorption spectral line
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import numpy as np
import copy
import pdb
import warnings

from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord

from linetools.analysis import utils as lau
from linetools.analysis import absline as laa
from linetools.lists.linelist import LineList
from linetools import utils as ltu

# Class for Spectral line
class SpectralLine(object):
    """
    Class for a spectral line. Emission or absorption.

    Parameters
    ----------
    ltype : str
        Type of Spectral line. Just `Abs` implemented for now but
        we expect to include `Emiss` in future.
    trans : Quantity, str
        Either the rest wavelength (e.g. 1215.6700*u.AA) or the
        transition name (e.g. 'CIV 1548')
    linelist : LineList, str; optional
        A LineList instance or string input to LineList (e.g. 'HI')
    closest : bool; optional
        Take the closest line to input wavelength (Default is False)

    Attributes
    ----------
    ltype : str
        Type of line, either 'Abs' or 'Emiss' (not currently implemented though)
    wrest : Quantity
        Rest wavelength
    z : float
        Redshift
    attrib : dict
        Line properties (e.g. column, EW, centroid, RA, Dec)
    analy : dict
        Analysis inputs (e.g. a spectrum, wavelength limits)
    data : dict
        Line atomic/molecular data (e.g. f-value, A coefficient, Elow)
    """

    @classmethod
    def from_dict(cls, idict, warn_only=False):
        """ Initialize from a dict (usually read from disk)

        Parameters
        ----------
        idict : dict
          dict with the Line parameters
        warn_only : bool, optional

        Returns
        -------
        sline : SpectralLine
         SpectralLine of the proper type

        """
        # Init
        if idict['ltype'] == 'Abs':
            # TODO: remove this try/except eventually
            try:
                sline = AbsLine(idict['name'])
            except: #  This is to be compatible JSON files already written with old notation
                sline = AbsLine(idict['trans'])
        else:
            raise ValueError("Not prepared for type {:s}.".format(idict['ltype']))
        # Check data
        for key in idict['data']:
            if isinstance(idict['data'][key], dict):  # Assume Quantity
                val = idict['data'][key]['value']
            else:
                val = idict['data'][key]
            try:
                assert sline.data[key] == val
            except AssertionError:
                if warn_only:
                    warnings.warn("Different data value for {:s}: {}, {}".format(key,sline.data[key],val))
        # Set analy
        for key in idict['analy']:
            if isinstance(idict['analy'][key], dict):  # Assume Quantity
                #sline.analy[key] = Quantity(idict['analy'][key]['value'],
                #                             unit=idict['analy'][key]['unit'])
                sline.analy[key] = ltu.convert_quantity_in_dict(
                        idict['analy'][key])
            elif key == 'spec_file':
                warnings.warn("You will need to load {:s} into attrib['spec'] yourself".format(
                        idict['analy'][key]))
            else:
                sline.analy[key] = idict['analy'][key]
        # Set attrib
        for key in idict['attrib']:
            if isinstance(idict['attrib'][key], dict):
                sline.attrib[key] = ltu.convert_quantity_in_dict(
                        idict['attrib'][key])
            elif key in ['RA','DEC']:
                sline.attrib['coord'] = SkyCoord(ra=idict['attrib']['RA']*u.deg,
                                              dec=idict['attrib']['DEC']*u.deg)
            else:
                sline.attrib[key] = idict['attrib'][key]
        return sline

    # Initialize with wavelength
    def __init__(self, ltype, trans, linelist=None, closest=False, z=0.,
                 verbose=True):

        # Required
        self.ltype = ltype
        if ltype not in ['Abs']:
            raise ValueError('spec/lines: Not ready for type {:s}'.format(ltype))

        # Init
        if not isinstance(trans,(Quantity,basestring)):
            raise ValueError('Rest wavelength must be a Quantity or str')

        # Other
        self.data = {} # Atomic/Molecular Data (e.g. f-value, A coefficient, Elow)
        self.analy = {
            'spec': None,              # Analysis inputs (e.g. spectrum; from .clm file or AbsID)
            'wvlim': [0., 0.]*u.AA,    # Wavelength interval about the line (observed)
            'vlim': [0., 0.]*u.km/u.s, # Rest-frame velocity limit of line, relative to self.attrib['z']
            'flag_kin': 0,             # Use for kinematic analysis?
            'do_analysis': 1           # Analyze?
            }
        self.attrib = {
            'coord': SkyCoord(ra=0.*u.deg, dec=0.*u.deg),  # Coords
            'z': z, 'sig_z': 0.,                           # Redshift
            'v': 0.*u.km/u.s, 'sig_v': 0.*u.km/u.s,        # rest-frame velocity relative to z
            'EW': 0.*u.AA, 'sig_EW': 0.*u.AA, 'flag_EW': 0 # EW
            }

        # Fill data
        self.fill_data(trans, linelist=linelist, closest=closest, verbose=verbose)

    def ismatch(self, inp, Zion=None, RADec=None):
        """Query whether input line matches on:  z, Z, ion, RA, Dec

        Parameters
        ----------
        inp : SpectralLine or tuple
          * SpectralLine -- Other spectral line for comparison
          * tuple -- (z, wrest) float, Quantity
            e.g. (1.3123, 1215.670*u.AA)
        Zion : tuple of ints, optional
          Generally used with tuple input, e.g. (6,2)
        RADec : tuple of Quantities, optional
          Generally used with tuple input e.g. (124.132*u.deg, 29.231*u.deg)

        Returns
        -------
        answer : bool
          True if a match, else False
        """
        coord = None
        if isinstance(inp, SpectralLine):
            wrest = inp.wrest
            z = inp.attrib['z']
            if Zion is None:
                Zion = (inp.data['Z'], inp.data['ion'])
            if RADec is None:
                coord = inp.attrib['coord']
        elif isinstance(inp, tuple):
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

    def cut_spec(self, normalize=False):
        """ Cut out a chunk of the spectrum around this line.

        Parameters
        ----------
        normalize : bool
            Whether to normalize the spectrum (if continuum exists)
            Default is False

        Returns
        -------
        fx, sig, dict(wave, velo, pix)
            Arrays (numpy or Quantity) of flux, error, and wavelength/velocity
            The velocity is calculated relative to self.attrib['z']
        """

        # Checks
        if self.analy['spec'] is None:
            raise ValueError('spectralline.cut_spec: Need to set spectrum!')
        if self.analy['spec'].wavelength.unit == 1.:
            raise ValueError('Expecting a unit!')

        # Pixels for evaluation
        if np.sum(self.analy['wvlim'].value > 0.):
            pix = self.analy['spec'].pix_minmax(self.analy['wvlim'])[0]
        elif np.sum(np.abs(self.analy['vlim'].value) > 0.):
            pix = self.analy['spec'].pix_minmax(
                self.attrib['z'], self.wrest, self.analy['vlim'])[0]
        else:
            raise ValueError('spectralline.cut_spec: Need to set wvlim or vlim!')
        self.analy['pix'] = pix

        # Normalize?
        if normalize:
            sv_normed = self.analy['spec'].normed
            self.analy['spec'].normed = True

        # Cut for analysis
        fx = self.analy['spec'].flux[pix]
        sig = self.analy['spec'].sig[pix]
        wave = self.analy['spec'].wavelength[pix]

        # Velocity array created within the XSpectrum1D class and cut afterwards
        self.analy['spec'].velo = self.analy['spec'].relative_vel(
            self.wrest*(1 + self.attrib['z']))
        velo = self.analy['spec'].velo[pix]

        # Set it back
        if normalize:
            self.analy['spec'].normed = sv_normed

        # Return
        return fx, sig, dict(wave=wave, velo=velo, pix=pix)

    def measure_ew(self, flg=1, initial_guesses=None):
        """ Measures the observer frame equivalent width

        Note this requires the keys `wvlim` and `spec` in analy to
        be set! Default is simple boxcar integration.
        Observer frame, not rest-frame (use measure_restew()
        for rest-frame).

        It sets these attributes:
           * self.attrib['EW', 'sig_EW']:
             The EW and error in observer frame

        Parameters
        ----------
        flg : int
            * 1 -- Boxcar integration (default)
            * 2 -- Gaussian fit
        initial_guesses : tuple of floats [None]
          If a model is chosen (e.g. flg=2, Gaussian) a tuple of
          (amplitude, mean, stddev) can be specified.

        """
        # Cut spectrum
        fx, sig, xdict = self.cut_spec(normalize=True)
        wv = xdict['wave']

        # Calculate
        if flg == 1: # Boxcar
            EW, sig_EW = lau.box_ew( (wv, fx, sig) )
        elif flg == 2: #Gaussian
            EW, sig_EW = lau.gaussian_ew( (wv, fx, sig), self.ltype, initial_guesses=initial_guesses)
        else:
            raise ValueError('measure_ew: Not ready for this flag {:d}'.format(flg))

        # Fill
        self.attrib['EW'] = EW 
        self.attrib['sig_EW'] = sig_EW

    def measure_restew(self, **kwargs):
        """  Measure the rest-frame equivalent width

        See `~measure_ew` for details.
        """
        # Standard call
        self.measure_ew(**kwargs)

        # Push to rest-frame
        self.attrib['EW'] = self.attrib['EW'] / (self.attrib['z']+1)
        self.attrib['sig_EW'] = self.attrib['sig_EW'] / (self.attrib['z']+1)

    def measure_kin(self, **kwargs):
        """  Measure Kinematics
        """
        from linetools.analysis import abskin
        fx, sig, cdict = self.cut_spec()
        stau = abskin.generate_stau(cdict['velo'], fx, sig, **kwargs)
        # Measure
        kin_data = abskin.pw97_kin(cdict['velo'], stau, **kwargs)
        # Save
        self.attrib['kin'] = kin_data

    # EW
    def to_dict(self):
        """ Convert class to dict

        Returns
        -------
        adict : dict
         dict representation of the SpectralLine
        """
        from numpy.ma.core import MaskedConstant
        from astropy.units import Quantity
        # Starting
        adict = dict(ltype=self.ltype, analy=dict(), attrib=dict(), data=dict(),
                     name=self.name, wrest=dict(value=self.wrest.value,
                                                  unit=self.wrest.unit.to_string()))
        # Data
        for key in self.data:
            # Skip masked values
            if isinstance(self.data[key], MaskedConstant):
                continue
            # Quantity
            elif isinstance(self.data[key], Quantity):
                adict['data'][key] = dict(value=self.data[key].value,
                                          unit=self.data[key].unit.to_string())
            else:
                adict['data'][key] = self.data[key]
        # Attrib
        for key in self.attrib:
            if key == 'coord':
                adict['attrib']['RA'] = self.attrib['coord'].ra.value
                adict['attrib']['DEC'] = self.attrib['coord'].dec.value
            elif isinstance(self.attrib[key], Quantity):
                adict['attrib'][key] = dict(value=self.attrib[key].value,
                                            unit=self.attrib[key].unit.to_string())
            else:
                adict['attrib'][key] = self.attrib[key]
        # Analysis
        for key in self.analy:
            if key == 'spec':
                if self.analy['spec'] is not None:
                    adict['analy']['spec_file'] = self.analy['spec'].filename
            elif isinstance(self.analy[key], Quantity):
                adict['analy'][key] = dict(value=self.analy[key].value,
                                            unit=self.analy[key].unit.to_string())
            else:
                adict['analy'][key] = self.analy[key]

        # Polish for JSON
        adict = ltu.jsonify(adict)
        # Return
        return adict

    # Output
    def __repr__(self):
        txt = '<{:s}:'.format(self.__class__.__name__)
        try:
            txt = txt+' {:s},'.format(self.data['name'])
        except KeyError:
            pass
        txt = txt + ' wrest={:g}'.format(self.wrest)
        txt = txt + '>'
        return (txt)

# ###########################################
# Class for Generic Absorption Line System
class AbsLine(SpectralLine):
    """Class representing a spectral absorption line.

    Parameters
    ----------
    trans : Quantity or str
        Quantity -- Rest wavelength (e.g. 1215.6700*u.AA)
        str -- Name of transition (e.g. 'CIV 1548'). For an
        unknown transition use string 'unknown'.
    """
    # Initialize with a .dat file
    def __init__(self, trans, **kwargs):
        # Generate with type

        # need to use super here. (See
        # http://docs.astropy.org/en/stable/development/codeguide.html#super-vs-direct-example)
        super(AbsLine, self).__init__('Abs', trans, **kwargs)

    def print_specline_type(self):
        """ Return a string representing the type of vehicle this is."""
        return 'AbsLine'

    def fill_data(self, trans, linelist=None, closest=False, verbose=True):
        """ Fill atomic data and setup analy.

        Parameters
        ----------
        trans : Quantity or str
          Either a rest wavelength (e.g. 1215.6700*u.AA) or the name
          of a transition (e.g. 'CIV 1548'). For an unknown transition
          use string 'unknown'.
        linelist : LineList, optional
          Class of linelist or str setting LineList
        closest : bool, optional
          Take the closest line to input wavelength? [False]
        """

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
        self.name = self.data['name']

        #
        self.analy.update( {
            'flg_eye': 0,
            'flg_limit': 0, # No limit
            'datafile': '', 
            'name': self.data['name']
            })

        # Additional fundamental attributes for Absorption Line
        self.attrib.update({'N': 0./u.cm**2, 'sig_N': 0./u.cm**2, 'flag_N': 0, # Column
                       'b': 0.*u.km/u.s, 'sig_b': 0.*u.km/u.s  # Doppler
                       } )

    # Voigt
    def generate_voigt(self, wave=None, **kwargs):
        """ Generate a Voigt profile model for the absorption line in
        a given spectrum.

        Parameters
        ----------
        wave : Quantity array
          Wavelength array on which to calculate the line
          Must be set if self.analy['spec'] is not filled

        Returns
        -------
        spec : XSpectrum1D
          Spectrum with the input wavelength and the absorbed flux
        """
        from linetools.analysis import voigt as lav
        # Checks
        if self.attrib['N'] < 1./u.cm**2:
            raise ValueError("Need to initialize log column density in attrib['N']")
        if self.attrib['b'] < 1.*u.km/u.s:
            raise ValueError("Need to initialize Doppler parameter in attrib['b']")
        if wave is None:
            # Assume a spectrum has been loaded already
            try:
                wave = self.analy['spec'].wavelength
            except:
                raise ('You must provide a wavelength array in generate_voigt')

        # Main call
        spec = lav.voigt_from_abslines(wave, self, **kwargs)
        return spec

    # AODM
    def measure_aodm(self, nsig=3., normalize=True):
        """ AODM calculation

        It sets these attributes:
          * self.attrib[ 'N', 'sig_N', 'logN', 'sig_logN' ]:
            Column densities and errors, linear and log

        Parameters
        ----------
        nsig : float, optional
          Number of sigma significance required for a "detection"
        normalize : bool, optional
          Normalize first?
        """
        # Cut spectrum
        fx, sig, xdict = self.cut_spec(normalize=normalize)
        velo = xdict['velo']

        # Calculate
        N, sig_N, flg_sat = laa.aodm((velo, fx, sig), (self.wrest,self.data['f']))

        # Flag
        if flg_sat:
            self.attrib['flag_N'] = 2
        else:
            if N > nsig*sig_N:
                self.attrib['flag_N'] = 1
            else:
                self.attrib['flag_N'] = 3

        # Values
        self.attrib['N'] = N
        self.attrib['sig_N'] = sig_N

        # Log
        laa.log_clm(self.attrib)

    # Output
    def __repr__(self):
        txt = '<{:s}:'.format(self.__class__.__name__)
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
        txt = txt + '>'
        return (txt)

def many_abslines(all_wrest, llist):
    """Generate a list of AbsLine objects.

    Useful for when you have many lines (>1000) to generate that have
    similar wrest.  Uses deepcopy.

    Parameters
    ----------
    all_wrest : list of lines
    llist : LineList

    Returns
    -------
    abs_lines : list of AbsLine Objects
    """
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
