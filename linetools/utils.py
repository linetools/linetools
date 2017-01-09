""" Module for general utilities which don't belong in another sub-package.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import json
import gzip, os
import warnings
import pdb

import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity, Unit

try:
    basestring
except NameError:  # For Python 3
    basestring = str


def between(a, vmin, vmax):
    """ Return a boolean array True where vmin <= a < vmax.

    Parameters
    ----------
    a : array of shape N
      Input array.
    vmin, vmix: float
      Minimum and maximum values to test between.

    Returns
    -------
    c : array of shape N
      Boolean array true where vmin < a < vmax.

    Notes
    -----
    This is a convenience function equivalent to (vmin <= a) & (a < vmax).
    Be careful of floating point issues when dealing with equalities.
    """
    a = np.asarray(a)
    c = a < vmax
    c &= a >= vmin
    return c


def convert_quantity_in_dict(idict):
    """ Return a dict where Quantities (usually from a JSON file)
    have been converted from unit/value

    Parameters
    ----------
    idict : dict
      Input dict

    Returns
    -------
    obj : dict or Quantity
    """
    if 'unit' in idict.keys():  # Simple dict of Quantity (e.g. from jsonify)
        obj = Quantity(idict['value'], unit=idict['unit'])
        return obj
    else:  # Nested dict (possibly)
        for key in idict.keys():
            if isinstance(idict[key], dict):
                idict[key] = convert_quantity_in_dict(idict[key])
        return idict


def name_from_coord(coord, precision=(2,1)):
    """ Generate a standard JXXXXXX.XX+XXXXXX.X name from a SkyCoord object

    Parameters
    ----------
    coord : SkyCoord
    precision : tuple, optional
      Number of decimal places to include in name

    Returns
    -------
    name : str
      In JXX format
    """
    name = 'J{:s}{:s}'.format(coord.ra.to_string(unit=u.hour,sep='',pad=True,precision=precision[0]),
            coord.dec.to_string(sep='',pad=True,alwayssign=True,precision=precision[1]))
    # Return
    return name


def radec_to_coord(radec):
    """ Converts one of many of Celestial Coordinates
    `radec` formats to an astropy SkyCoord object. Assumes
    J2000 equinox.

    Parameters
    ----------
    radec : str or tuple or SkyCoord
        Examples:
        'J124511+144523',
        '124511+144523',
        'J12:45:11+14:45:23',
        ('12:45:11','+14:45:23')
        ('12 45 11', +14 45 23)
        ('12:45:11','14:45:23')  -- Assumes positive DEC
        (123.123, 12.1224) -- Assumed deg

    Returns
    -------
    coord : SkyCoord
      Converts to astropy.coordinate.SkyCoord (as needed)

    """
    from astropy.coordinates import SkyCoord

    # RA/DEC
    if isinstance(radec, (tuple)):
        if isinstance(radec[0], basestring):
            if radec[1][0] not in ['+', '-']:  #
                DEC = '+'+radec[1]
                warnings.warn("Assuming your DEC is +")
            else:
                DEC = radec[1]
            #
            coord = SkyCoord(radec[0]+DEC, frame='fk5',
                                  unit=(u.hourangle, u.deg))
        else:
            coord = SkyCoord(ra=radec[0], dec=radec[1], unit='deg')
    elif isinstance(radec,SkyCoord):
        coord = radec
    elif isinstance(radec,basestring):
        # Find first instance of a number (i.e. strip J, SDSS, etc.)
        for ii in range(len(radec)):
            if radec[ii].isdigit():
                break
        radec = radec[ii:]
        #
        if ':' in radec:
            coord = SkyCoord(radec, frame='fk5', unit=(u.hourangle, u.deg))
        else:  # Add in :
            if ('+' in radec) or ('-' in radec):
                sign = max(radec.find('+'), radec.find('-'))
            else:
                raise ValueError("radec must include + or - for DEC")
            newradec = (radec[0:2]+':'+radec[2:4]+':'+radec[4:sign+3] +':'+radec[sign+3:sign+5]+':'+radec[sign+5:])
            coord = SkyCoord(newradec, frame='fk5', unit=(u.hourangle, u.deg))
    # Return
    return coord


def scipy_rebin(a, *args):
    """ Simple script to rebin an input array to a new shape.

    Akin to IDL's routine Taken from scipy documentation:
    http://wiki.scipy.org/Cookbook/Rebinning As in IDL, the new shape
    must be a factor of the old one.  The ugly 'evList trick' builds
    and executes a python command.

    """
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
    #print ''.join(evList)
    return eval(''.join(evList))


def jsonify(obj, debug=False):
    """ Recursively process an object so it can be serialised in json
    format.

    WARNING - the input object may be modified if it's a dictionary or
    list!

    Parameters
    ----------
    obj : any object
    debug : bool, optional

    Returns
    -------
    obj - the same obj is json_friendly format (arrays turned to
    lists, np.int64 converted to int, np.float64 to float, and so on).

    """
    if isinstance(obj, np.float64):
        obj = float(obj)
    elif isinstance(obj, np.float32):
        obj = float(obj)
    elif isinstance(obj, np.int32):
        obj = int(obj)
    elif isinstance(obj, np.int64):
        obj = int(obj)
    elif isinstance(obj, np.int16):
        obj = int(obj)
    elif isinstance(obj, np.bool_):
        obj = bool(obj)
    elif isinstance(obj, np.string_):
        obj = str(obj)
    elif isinstance(obj, Quantity):
        obj = dict(value=obj.value, unit=obj.unit.to_string())
    elif isinstance(obj, np.ndarray):  # Must come after Quantity
        obj = obj.tolist()
    elif isinstance(obj, dict):
        for key, value in obj.items():
            obj[key] = jsonify(value, debug=debug)
    elif isinstance(obj, list):
        for i,item in enumerate(obj):
            obj[i] = jsonify(item, debug=debug)
    elif isinstance(obj, tuple):
        obj = list(obj)
        for i,item in enumerate(obj):
            obj[i] = jsonify(item, debug=debug)
        obj = tuple(obj)
    elif isinstance(obj, Unit):
        obj = obj.name
    elif obj is u.dimensionless_unscaled:
        obj = 'dimensionless_unit'

    if debug:
        print(type(obj))
    return obj


def savejson(filename, obj, overwrite=False, indent=None, easy_to_read=False,
             **kwargs):
    """ Save a python object to filename using the JSON encoder.

    Parameters
    ----------
    filename : str
    obj : object
      Frequently a dict
    overwrite : bool, optional
    indent : int, optional
      Input to json.dump
    easy_to_read : bool, optional
      Another approach and obj must be a dict
    kwargs : optional
      Passed to json.dump

    Returns
    -------

    """
    import io

    if os.path.lexists(filename) and not overwrite:
        raise IOError('%s exists' % filename)
    if easy_to_read:
        if not isinstance(obj, dict):
            raise IOError("This approach requires obj to be a dict")
        with io.open(filename, 'w', encoding='utf-8') as f:
            f.write(json.dumps(obj, sort_keys=True, indent=4,
                               separators=(',', ': '), **kwargs))
    else:
        if filename.endswith('.gz'):
            with gzip.open(filename, 'wt') as fh:
                json.dump(obj, fh, indent=indent, **kwargs)
        else:
            with open(filename, 'wt') as fh:
                json.dump(obj, fh, indent=indent, **kwargs)


def loadjson(filename):
    """ Load a python object saved with savejson."""
    if filename.endswith('.gz'):
        with gzip.open(filename, "rb") as f:
            obj = json.loads(f.read().decode("ascii"))
    else:
        with open(filename, 'rt') as fh:
            obj = json.load(fh)

    return obj


def rel_vel(wavelength, wv_obs):
    """ Simple relative velocity method

    Parameters
    ----------
    wavelength : Quantity array
    wv_obs : Quantity

    Returns
    -------

    """
    if not isinstance(wavelength, Quantity):
        raise ValueError('Input wavelength array needs to be a Quantity array')
    if not isinstance(wv_obs, Quantity):
        raise ValueError('Input wv_obs needs to be a Quantity')
    return ((wavelength - wv_obs) * const.c / wv_obs).to('km/s')


def v_from_z(z1, z2):
    """ Find the relativistic velocity between 2 redshifts.

    Parameters
    ----------
    z1 : float
       One redshift.
    z2 : float or array
       Other redshift(s)

    Returns
    -------
    v : Quantity (km/s)
      Velocity

    Notes
    -----
    """

    """
    R = (1+z1) / (1+z2)
    v = const.c * (R**2 - 1)/(1+R**2)

    return v.to('km/s')"""
    raise DeprecationWarning("This function is deprecated, please use dv_from_z() instead.")


def z_from_v(z, v):
    """ Find the redshift given z and v

    Parameters
    ----------
    z : float or array
       Redshift
    v : Quantity or float or array or array of Quantity
       Velocities. If not Quantity it assumes km/s units.

    Returns
    -------
    z : float or array
      New redshifts

    Notes
    -----
    """

    """
    # Check for unit
    if not isinstance(v, u.quantity.Quantity):
        # Assume km/s
        v = v * u.Unit('km/s')

    # b
    bval = (v / const.c.to('km/s'))

    # R
    R = np.sqrt((1-bval)/(1+bval))
    # Finally
    znew = (1+z)/R - 1

    return znew.value
    """
    raise DeprecationWarning("This function is deprecated, instead please use either z_from_dv() or dz_from_dv()"
                             " depending on your needs.")


# Slightly different functions for passing from dv to dz, and viceversa (that NT prefers).
# May need to agree on one kind of conversion in the future
def give_dz(dv, zref, rel=True):
    """Same as dz_from_dv. This function will be deprecated."""
    DeprecationWarning("This function will be deprecated. Please use instead dz_from_dv().")
    raise DeprecationWarning('Gone :P')
    return dz_from_dv(dv, zref, rel=rel)


def give_dv(z, zref, rel=True):
    """Same as dv_from_dz. This function will be deprecated."""
    warnings.warn("This function will be deprecated. Please use instead dv_from_z().")
    raise DeprecationWarning('Gone :P')
    return dv_from_z(z, zref, rel=rel)


def dv_from_z(z, zref, rel=True):
    """Gives the rest-frame velocity difference dv
    between z and zref. dv=0 at zref by definition.

    Parameters
    ----------
    z : float or np.ndarray or list or tuple
        Redshifts to calculate dv on
    zref : float or np.ndarray or list or tuple
        Reference redshift where dv=0 by definition.
        If the shape of zref is equal to the shape of z,
        each dv is calculated at each zref, otherwise zref
        is expected to be float.
    rel : bool, optional
        Whether to apply relativistic correction for
        a locally flat space-time. Default is True.

    Returns
    -------
    dv : Quantity or Quantity array
        Rest-frame velocity difference between z and zref (dv=0 at zref by definition).
        It has the same shape as z.
        """
    # check format
    if not isinstance(z, (float, np.ndarray, list, tuple)):
        raise IOError('z must be float or np.ndarray or list or tuple.')
    if not isinstance(zref, (float, np.ndarray, list, tuple)):
        raise IOError('zref must be float or np.ndarray or list or tuple.')
    if (not isinstance(zref, float)) and (np.shape(zref) != np.shape(z)):
        raise IOError('If zref is not float, it must be of same shape as z.')

    z = np.array(z)
    zref = np.array(zref)

    if rel:
        dv = ((1 + z)**2 - (1 + zref)**2) / ((1 + z)**2 + (1 + zref)**2)
    else:
        dv = (z - zref) / (1. + zref)

    return dv * const.c.to('km/s')


def dz_from_dv(dv, zref, rel=True):
    """Gives redshift difference for a given
    velocity difference(s) with respect to zref.

    Parameters
    ----------
    dv : Quantity or Quantity array
        Rest-frame velocity difference with respect to zref
    zref : float or np.ndarray or list or tuple
        Reference redshift where dv=0.
        If shape of zref is equal than shape of dv,
        each dz is calculated at each zref, otherwise zref
        is expected to be float
    rel : bool, optional
        Whether to apply relativistic correction for
        a locally flat space-time. Default is True.

    Returns
    -------
    dz : np.array
        Redshift difference for a given dv with respect to zref.
        Same shape as dv.

    Notes
    -----
    See also linetools.utils.z_from_dv()
    """
    if not isinstance(dv, u.quantity.Quantity):
        raise IOError('dv must be Quantity or Quantity array.')
    if not isinstance(zref, (float, np.ndarray, list, tuple)):
        raise IOError('zref must be float or np.ndarray or list or tuple.')
    if (not isinstance(zref, float)) and (np.shape(zref) != np.shape(dv)):
        raise IOError('If zref is not float, it must be of same shape as dv.')

    zref = np.array(zref)

    beta = dv / const.c
    beta = beta.decompose()
    # check dimensionless
    if beta.unit != u.dimensionless_unscaled:
        raise IOError('dv must have velocity units.')
    beta = beta.value  # beta is dimensionless

    if rel:
        aux = np.sqrt((1. + beta) / (1. - beta))
        dz = (1. + zref) * (aux - 1.)
    else:
        dz = beta * (1. + zref)
    return dz


def z_from_dv(dv, zref, rel=True):
    """Gives the redshift for a given
    velocity difference(s) with respect to zref.

    Parameters
    ----------
    dv : Quantity or Quantity array
        Rest-frame velocity difference with respect to zref
    zref : float or np.ndarray or list or tuple
        Reference redshift where dv=0.
        If shape of zref is equal than shape of dv,
        each z is calculated at each zref, otherwise zref
        is expected to be float
    rel : bool, optional
        Whether to apply relativistic correction for
        a locally flat space-time. Default is True.

    Returns
    z : np.array
        Absolute redshift for a given dv with respect to zref.
        Same shape as dv.

    Notes
    -----
    See also linetools.utils.dz_from_dv()
    """
    return dz_from_dv(dv, zref, rel=rel) + zref


def overlapping_chunks(chunk1, chunk2):
    """True if there is overlap between chunks
    `chunk1` and `chunk2`. Otherwise False. Chunks are
    assumed to represent continuous coverage, so the only
    information that matters are the minimum and maximum
    values of a given chunk. Chunks must be sorted though.

    Parameters
    ----------
    chunk1 : tuple, list, 1-d np.array, Quantity, Quantity array
        A given chunk, assumed to represent a contiguous region
        so only its minimum and maximum values matter. Still,
        chunk must be sorted.
    chunk2 : tuple, list, 1-d np.array, Quantity, Quantity array
        Ditto.

    Returns
    -------
    answer : bool
        True if there is overlap, False otherwise.

    """
    # Check units in case chunks are Quantity
    if isinstance(chunk1, Quantity):
        unit1 = chunk1.unit
        chunk1 = np.array(chunk1.value)
        if not isinstance(chunk2, Quantity):
            raise ValueError('chunk2 must be Quantity because chunk1 is!')
        try:
            chunk2 = chunk2.to(unit1)
            chunk2 = np.array(chunk2.value)  # has the same units as chunk1
        except u.core.UnitConversionError:
            raise ValueError('If chunks are given as Quantity they must have convertible units!')
    else:
        chunk1 = np.array(chunk1)  # not a quantity

    # this may be redundant but cleaner code
    if isinstance(chunk2, Quantity):
        unit2 = chunk2.unit
        chunk2 = np.array(chunk2.value)
        if not isinstance(chunk1, Quantity):
            raise ValueError('chunk1 must be Quantity because chunk2 is!')
    else:
        chunk2 = np.array(chunk2)  # not a quantity
    # here we have chunks as values (i.e no units)

    # make sure the chunks are sorted
    cond1 = np.sort(chunk1) != chunk1
    cond2 = np.sort(chunk2) != chunk2
    if (np.sum(cond1) > 0) or (np.sum(cond2) > 0):
        raise ValueError('chunks must be sorted!')

    # figure out the lowest of the chunks
    # and sort them such that chunk1 is by definition
    # the one with the lowest limit
    if np.min(chunk1) <= np.min(chunk2):
        pass
    else: # invert them if necessary
        aux = chunk1
        chunk1 = chunk2
        chunk2 = aux

    # now the chunk1 is the one with the lowest value
    # then, the only way they do not overlap is:
    if np.min(chunk2) > np.max(chunk1):
        return False
    else:  # overlap exists
        return True