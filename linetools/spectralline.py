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
from linetools.analysis.linelimits import LineLimits
from linetools.lists.linelist import LineList
from linetools import utils as ltu
from linetools.spectra.xspectrum1d import XSpectrum1D

# A few globals to speed up performance (astropy.Quantity issue)
zero_coord = SkyCoord(ra=0.*u.deg, dec=0.*u.deg)  # Coords
init_analy = {
            'spec': None,              # Analysis inputs (e.g. spectrum; from .clm file or AbsID)
            'flag_kin': 0,             # Use for kinematic analysis?
            'do_analysis': 1           # Analyze?
            }
init_attrib = {
            'coord': zero_coord,                           # Coords
            'EW': 0.*u.AA, 'sig_EW': 0.*u.AA, 'flag_EW': 0 # EW
            }

abs_attrib = {'N': 0./u.cm**2, 'sig_N': 0./u.cm**2, 'flag_N': 0, # Column    ## NOT ENOUGH SPEED-UP
              'logN': 0., 'sig_logN': 0.,
                    'b': 0.*u.km/u.s, 'sig_b': 0.*u.km/u.s  # Doppler
                    }

emiss_attrib = {'flux': 0.*u.erg/u.s, 'sig_flux': 0.*u.erg/u.s, 'flag_flux': 0,
                }


class SpectralLine(object):
    """ Class for a spectral line. Emission or absorption.

    Parameters
    ----------
    ltype : str
        Type of Spectral line. `Abs` and `Em` implemented
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
        Type of line, either 'Abs' or 'Em'
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
    limits : LineLimits
        Limits including zlim, vlim, wvlim.
    """

    @classmethod
    def from_dict(cls, idict, coord=None, warn_only=False, chk_data=True, **kwargs):
        """ Initialize from a dict (usually read from disk)

        Parameters
        ----------
        idict : dict
          dict with the Line parameters
        chk_data : bool, optional
          Check atomic data in dict against current values in LineList
        warn_only : bool, optional
          If the chk_data is performed and the values do not match, only
          throw a Warning as opposed to crashing

        Returns
        -------
        sline : SpectralLine
         SpectralLine of the proper type

        """
        # Init
        if idict['ltype'] == 'Abs':
            # TODO: remove this try/except eventually
            try:
                sline = AbsLine(idict['name'], **kwargs)
            except KeyError: #  This is to be compatible JSON files already written with old notation (e.g. DLA H100)
                sline = AbsLine(idict['trans'], **kwargs)
        elif idict['ltype'] == 'Em':
            sline = EmLine(idict['name'], **kwargs)
        else:
            raise ValueError("Not prepared for type {:s}.".format(idict['ltype']))
        # Check data
        if chk_data:
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
        for key in idict['analy'].keys():
            if isinstance(idict['analy'][key], dict):  # Assume Quantity
                #sline.analy[key] = Quantity(idict['analy'][key]['value'],
                #                             unit=idict['analy'][key]['unit'])
                #pdb.set_trace()
                sline.analy[key] = ltu.convert_quantity_in_dict(idict['analy'][key])
            elif key == 'spec_file':
                # spec_file is intended to be the name of the spectrum file
                # spec is intendended to hold an XSpectrum1D object
                warnings.warn("You will need to load {:s} into analy['spec'] yourself".format(
                        idict['analy'][key]))
                sline.analy[key] = idict['analy'][key]
            else:
                sline.analy[key] = idict['analy'][key]

        # Set attrib
        for key in idict['attrib'].keys():
            if isinstance(idict['attrib'][key], dict):
                sline.attrib[key] = ltu.convert_quantity_in_dict(idict['attrib'][key])
            elif key in ['RA','DEC']:
                if coord is None:
                    sline.attrib['coord'] = SkyCoord(ra=idict['attrib']['RA']*u.deg,
                                                  dec=idict['attrib']['DEC']*u.deg)
                else:
                    sline.attrib['coord'] = coord
            else:
                sline.attrib[key] = idict['attrib'][key]

        # Set z and limits
        if 'z' in sline.attrib.keys():  # Backwards compatability
            z = sline.attrib.pop('z')
        else:
            z = 0.
        if 'limits' in idict.keys():
            if 'wrest' not in idict['limits'].keys(): # compatibility with IGMGuesses
                # import pdb; pdb.set_trace()
                idict['limits']['wrest'] = ltu.jsonify(sline.wrest)
                idict['limits']['z'] = z
            sline.limits = LineLimits.from_dict(idict['limits'])
        else:
            sline.limits = LineLimits(sline.wrest, z, [z,z])
            if 'vlim' in sline.analy.keys():  # Backwards compatability
                if sline.analy['vlim'][1] > sline.analy['vlim'][0]:
                    sline.limits.set(sline.analy['vlim'])
            elif 'wvlim' in sline.analy.keys():  # Backwards compatability
                if sline.analy['wvlim'][1] > sline.analy['wvlim'][0]:
                    sline.limits.set(sline.analy['wvlim'])
        return sline

    # Initialize with wavelength
    def __init__(self, ltype, trans, linelist=None, closest=False, z=0.,
                 verbose=True, **kwargs):

        # Required
        self.ltype = ltype
        if ltype not in ['Abs', 'Em']:
            raise ValueError('spec/lines: Not ready for type {:s}'.format(ltype))

        # Init
        if not isinstance(trans,(Quantity,basestring)):
            raise ValueError('Rest wavelength must be a Quantity or str')

        # Other
        self.data = {} # Atomic/Molecular Data (e.g. f-value, A coefficient, Elow)
        self.analy = init_analy.copy()
        self.attrib = init_attrib.copy()

        # Fill data
        self.fill_data(trans, linelist=linelist, closest=closest, verbose=verbose)
        # Limits
        try:
            zlim = kwargs['zlim']
        except KeyError:
            zlim = [z,z]
        if ltype in ['Abs', 'Em']:
            self.limits = LineLimits.from_specline(self, z, zlim)
        else:
            raise ValueError('Not ready to set limits for this type')

    @property
    def z(self):
        """ Return z
        """
        return self.limits.z

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
            if self.ltype == 'Abs':
                llist = LineList('ISM')
            elif self.ltype == 'Em':
                llist = LineList('Galaxy')
            else:
                raise ValueError("Not ready for ltype = {:s}".format(self.ltype))
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
        if newline is None:
            raise ValueError("Transition {} not found in LineList {:s}".format(trans, llist.list))
        try:
            self.data.update(newline)  # Expected to be a LineList dict object
        except TypeError:
            raise TypeError("Probably should not be here")


        # Update
        self.wrest = self.data['wrest']
        self.name = self.data['name']

        #
        self.update()  # is this used ?

    def setz(self, z):
        """ Set redshift wherever it is needed/expected
        Parameters
        ----------
        z : float

        Returns
        -------

        """
        if not isinstance(z,float):
            raise IOError("Input redshift needs to be a float")
        # Set
        self.limits._z = z
        # Warning?
        if self.limits.is_set():
            warnings.warn("Consider whether to update the limits of this line")


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
            z = inp.z
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
            np.allclose(self.z, z, rtol=1e-6))
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
            The velocity is calculated relative to self.z
        """

        # Checks
        if self.analy['spec'] is None:
            raise ValueError('spectralline.cut_spec: Need to set spectrum!')
        if self.analy['spec'].wavelength.unit == 1.:
            raise ValueError('Expecting a unit!')

        # Pixels for evaluation
        if self.limits.is_set():
            pix = self.analy['spec'].pix_minmax(self.limits.wvlim)[0]
        else:
            raise ValueError('spectralline.cut_spec: Need to set limits!')
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
            self.wrest*(1 + self.z))
        velo = self.analy['spec'].velo[pix]

        # Set it back
        if normalize:
            self.analy['spec'].normed = sv_normed

        # Return
        return fx, sig, dict(wave=wave, velo=velo, pix=pix)

    def measure_ew(self, flg=1, initial_guesses=None):
        """ Measures the observer frame equivalent width

        Note this requires self.limits to be initialized
        Default is simple boxcar integration.
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
        # Check that there is sufficient data
        if len(fx) <= 1:
            warnings.warn("Spectrum does not cover {:g}".format(self.wrest))
            self.attrib['EW'] = 0.
            self.attrib['sig_EW'] = -1
            return

        # Calculate
        if flg == 1: # Boxcar
            EW, sig_EW = lau.box_ew( (wv, fx, sig) )
        elif flg == 2: #Gaussian
            EW, sig_EW = lau.gaussian_ew( (wv, fx, sig), self.ltype, initial_guesses=initial_guesses)
        else:
            raise ValueError('measure_ew: Not ready for this flag {:d}'.format(flg))

        # Fill
        self.attrib['flag_EW'] = 1
        self.attrib['EW'] = EW
        self.attrib['sig_EW'] = sig_EW

    def measure_restew(self, **kwargs):
        """  Measure the rest-frame equivalent width

        See `~measure_ew` for details.
        """
        # Standard call
        self.measure_ew(**kwargs)

        # Push to rest-frame
        self.attrib['EW'] = self.attrib['EW'] / (self.z+1)
        self.attrib['sig_EW'] = self.attrib['sig_EW'] / (self.z+1)

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
        adict = dict(ltype=self.ltype, analy=dict(), attrib=dict(), data=dict(),\
                     limits=dict(), name=self.name, wrest=dict(value=self.wrest.value,\
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
                if isinstance(self.analy['spec'], basestring):
                    adict['analy']['spec_file'] = self.analy['spec']
                elif isinstance(self.analy['spec'], XSpectrum1D):
                    adict['analy']['spec_file'] = self.analy['spec'].filename
                else:
                    pass
            elif isinstance(self.analy[key], Quantity):
                adict['analy'][key] = dict(value=self.analy[key].value,
                                            unit=self.analy[key].unit.to_string())
            else:
                adict['analy'][key] = self.analy[key]

        # Limits
        adict['limits'] = self.limits.to_dict()

        # Polish for JSON
        adict = ltu.jsonify(adict)
        # Return
        return adict

    def coincident_line(self, specline):
        """Whether the current SpectralLine overlaps in
        observed wavelength space with the given input SpectralLine

        Parameters
        ----------
        specline : SpectralLine
            A SpectralLine object

        Returns
        -------
        answer : bool
          True if there is overlap in wvobs space, False otherwise.

        """
        if not self.limits.is_set():
            raise ValueError("{} has not set its limits!".format(self.__repr__()))
        if not specline.limits.is_set():
            raise ValueError("{} has not set its limits!".format(specline.__repr__()))

        return ltu.overlapping_chunks(self.limits.wvlim, specline.limits.wvlim)

    def copy(self):
        """ Generate a copy

        Returns
        -------
        newline : SpectralLine
          copy of the object
        """
        return copy.deepcopy(self)  # Cheat for now

    def __repr__(self):
        txt = '<{:s}:'.format(self.__class__.__name__)
        try:
            txt = txt+' {:s},'.format(self.data['name'])
        except KeyError:
            pass
        # z
        txt = txt + ' z={:.4f}'.format(self.z)
        #
        txt = txt + ' wrest={:g}'.format(self.wrest)
        txt = txt + '>'
        return (txt)


class AbsLine(SpectralLine):
    """Class representing a spectral absorption line.

    Parameters
    ----------
    trans : Quantity or str
        Quantity -- Rest wavelength (e.g. 1215.6700*u.AA)
        str -- Name of transition (e.g. 'CIV 1548'). For an
        unknown transition use string 'unknown'.
    """
    def __init__(self, trans, **kwargs):
        # Generate with type

        # need to use super here. (See
        # http://docs.astropy.org/en/stable/development/codeguide.html#super-vs-direct-example)
        super(AbsLine, self).__init__('Abs', trans, **kwargs)

    def print_specline_type(self):
        """ Return a string representing the type of vehicle this is."""
        return 'AbsLine'

    def update(self):
        """ Fill in a few additional dicts
        """
        self.analy.update( {
            'flg_eye': 0,
            'flg_limit': 0, # No limit
            'datafile': '',
            'name': self.data['name']
        })

        # Additional fundamental attributes for Absorption Line
        self.attrib.update(abs_attrib.copy())

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
        # Check that there is sufficient data
        if len(fx) <= 1:
            warnings.warn("Spectrum does not cover {:g}".format(self.wrest))
            self.attrib['flag_N'] = 0
            return

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

    def get_tau0(self, N, b):
        """It returns the optical depth at the line center, tau0,
        for a given column density and Doppler parameter. It uses
        approximation given in Draine 2011 (see Chapter 9).
        It neglects stimulated emission which is fine for IGM or ISM
        except for radio-frequency transitions.

        Parameters
        ----------
        N : Quantity or Quantity array
            Column density
        b : Quantity or Quantity array of same shape as N
            Doppler parameter

        Returns
        -------
        tau0: float or array
            Optical depth at the line center. If N and b are
            arrays they must be of same shape.

        Notes
        -----
        This is a wrapper to linetools.analysis.absline.get_tau0()
        """
        try:
            fosc = self.data['f']
        except KeyError:
            raise NotImplementedError('AbsLine {} has not set its oscillator strength.'.format(self.__repr__))
        return laa.get_tau0(self.wrest, fosc, N, b)

    def get_Wr_from_N_b(self, N, b):
        """It returns the rest-frame equivalent width for a given
        N and b. It uses the approximation given by Draine 2011 book
        (eq. 9.27), which comes from atomic physics considerations
        See also Rodgers & Williams 1974 (NT: could not find the reference
        given by Draine)

        Parameters
        ----------
        N : Quantity or Quantity array
            Column density
        b : Quantity or Quantity array of same shape as N
            Doppler parameter

        Returns
        -------
        Wr : Quantity
            Rest-frame equivalent width

        Notes
        -----
        This is a wrapper to linetools.analysis.absline.Wr_from_N_b().
        See also linetools.analysis.absline.Wr_from_N_b_transition().
        """
        try:
            fosc = self.data['f']
        except KeyError:
            raise NotImplementedError('AbsLine {} has not set its oscillator strength.'.format(self.__repr__))
        try:
            gamma = self.data['gamma']
        except KeyError:
            raise NotImplementedError('AbsLine {} has not set its gamma value.'.format(self.__repr__))
        return laa.Wr_from_N_b(N, b, self.wrest, fosc, gamma)

    def get_Wr_from_N(self, N):
        """It returns the approximated rest-frame equivalent width for
        a given N. It uses the approximation given by Draine 2011 book
        (eq. 9.15), which is valid for tau0<<1 where Wr is independent
        of Doppler parameter or gamma. Use get_Wr_from_N_b() for a
        better approximation valid for wider range in tau0.

        Parameters
        ----------
        N : Quantity or Quantity array
            Column density

        Returns
        -------
        Wr : Quantity
            Approximated rest-frame equivalent width valid for tau0<<1.

        Notes
        -----
        This is a wrapper to linetools.analysis.absline.Wr_from_N().
        See also get_Wr_from_N_b() and linetools.analysis.absline.Wr_from_N_transition().
        """
        try:
            fosc = self.data['f']
        except KeyError:
            raise NotImplementedError('AbsLine {} has not set its oscillator strength.'.format(self.__repr__))
        return laa.Wr_from_N(N, self.wrest, fosc)

    def get_N_from_Wr(self, Wr):
        """It returns the approximated column density N, for a given rest-frame equivalent width
        Wr. This is an approximation only valid for tau0 << 1, where
        Wr is independent on Doppler parameter and gamma (see eqs. 9.14 and 9.15 of
        Draine 2011). This may be useful to put upper limits on non-detections.

        Parameters
        ----------
        Wr : Quantity or Quantity array
            Rest-frame equivalent width of the AbsLine

        Returns
        -------
        N : Quantity
            Approximated column density N valid for the tau0<<1 regime.

        Notes
        -----
        This is a wrapper to linetools.analysis.absline.Wr_from_N(). See also self.get_Wr_from_N()
        """
        try:
            fosc = self.data['f']
        except KeyError:
            raise NotImplementedError('AbsLine {} has not set its oscillator strength.'.format(self.__repr__))
        return laa.N_from_Wr(Wr, self.wrest, fosc)

    def __repr__(self):
        txt = '<{:s}:'.format(self.__class__.__name__)
        # Name
        try:
            txt = txt+' {:s},'.format(self.data['name'])
        except KeyError:
            pass
        # z
        txt = txt + ' z={:.4f}'.format(self.z)
        # wrest
        txt = txt + ' wrest={:.4f}'.format(self.wrest)
        # fval
        try:
            txt = txt+', f={:g}'.format(self.data['f'])
        except (KeyError, ValueError):
            pass
        txt = txt + '>'
        return (txt)


class EmLine(SpectralLine):
    """Class representing a spectral emission line.

    Parameters
    ----------
    trans : Quantity or str
        Quantity -- Rest wavelength (e.g. 1215.6700*u.AA)
        str -- Name of transition (e.g. 'CIV 1548'). For an
        unknown transition use string 'unknown'.
    """
    def __init__(self, trans, **kwargs):
        # Generate with type
        super(EmLine, self).__init__('Em', trans, **kwargs)

    def update(self):
        """ Fill in a few additional dicts
        """
        #
        self.analy.update( {
            'datafile': '',
            'name': self.data['name']
            })

        # Additional fundamental attributes for Emission Line
        self.attrib.update(emiss_attrib.copy())


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
