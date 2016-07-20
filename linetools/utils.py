""" Module for general utilities which don't belong in another sub-package.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import json
import gzip, os
import warnings
import pdb

import numpy as np
from astropy.constants import c
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


def radec_to_coord(radec):
    """ Converts one of many of Celestial Coordinates
    `radec` formats to an astropy SkyCoord object. Assumes
    J2000 equinox.

    Parameters
    ----------
    radec : str or tuple
        Examples:
        'J124511+144523',
        '124511+144523',
        'J12:45:11+14:45:23',
        ('12:45:11','+14:45:23')
        ('12 45 11', +14 45 23)
        ('12:45:11','14:45:23')  -- Assumes positive DEC

    Returns
    -------
    coord : SkyCoord

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


def savejson(filename, obj, overwrite=False, indent=None):
    """ Save a python object to filename using using the JSON encoder."""

    if os.path.lexists(filename) and not overwrite:
        raise IOError('%s exists' % filename)
    if filename.endswith('.gz'):
        with gzip.open(filename, 'wt') as fh:
            json.dump(obj, fh, indent=indent)
    else:
        with open(filename, 'wt') as fh:
            json.dump(obj, fh, indent=indent)


def loadjson(filename):
    """ Load a python object saved with savejson."""
    if filename.endswith('.gz'):
        with gzip.open(filename, "rb") as f:
            obj = json.loads(f.read().decode("ascii"))
    else:
        with open(filename, 'rt') as fh:
            obj = json.load(fh)

    return obj


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
    R = (1+z1) / (1+z2)
    v = c * (R**2 - 1)/(1+R**2)

    return v.to('km/s')

# #####
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
    # Check for unit
    if not isinstance(v, u.quantity.Quantity):
        # Assume km/s
        v = v * u.Unit('km/s')

    # b
    bval = (v/c.to('km/s'))

    # R
    R = np.sqrt((1-bval)/(1+bval))
    # Finally
    znew = (1+z)/R - 1

    return znew.value

# Slightly different functions for passing from dv to dz, and viceversa (that NT prefers).
# May need to agree on one kind of conversion in the future
def give_dv(z, zmean, rel=True):
    """Gives velocity difference between z and zmean.

    Parameters
    ---------
    z : float or np.array
        Redshifts to calculate dv on
    zmean : float or np.array
        Rest-frame redshift to perform the calculation.
        If shape of zmean is equal than shape of z,
        each dv is calculated at each zmean, otherwise zmean
        is expected to be float
    rel : bool, optional
        Whether to apply relativistic correction for
        a locally flat space-time. Default is True.

    Returns
    -------
    dv : Quantity or Quantity array
        Rest frame velocity difference between z and zmean, at
        zmean. It has same shape as z.
        """
    z = np.array(z)
    zmean = np.array(zmean)

    if rel:
        dv = ((1 + z)**2 - (1 + zmean)**2) / ((1 + z)**2 + (1 + zmean)**2)
    else:
        dv = (z - zmean) / (1. + zmean)

    return dv * c.to('km/s')

def give_dz(dv, zmean, rel=True):
    """Gives redshift difference for a given
    velocity difference(s) at zmean.

    Parameters
    ---------
    dv : Quantity or Quantity array
        Rest-frame velocity at zmean to calculate
        the corresponding redshift difference, dz
    zmean : float or np.array
        Rest-frame redshift to perform the calculation.
        If shape of zmean is equal than shape of dv,
        each dv is calculated at each zmean, otherwise zmean
        is expected to be float
    rel : bool, optional
        Whether to apply relativistic correction for
        a locally flat space-time. Default is True.

    Returns
    dz : np.array
        Redshift difference between dv and zmean. Same shape as
        dv.
        """
    if not isinstance(dv, u.quantity.Quantity):
        raise ValueError('dv must be Quantity or Quantity array!')
    try:
        dv = dv.to('km/s') # dv in km/s
        dv = np.array(dv.value)
    except UnitConversionError:
        raise ValueError('dv must have velocity units!')

    zmean = np.array(zmean)

    if rel:
        beta = dv / c.to('km/s').value
        aux = np.sqrt((1.+ beta)/(1.- beta))
        dz = (1. + zmean) * (aux - 1.)
    else:
        dz = dv * (1. + zmean) / c.to('km/s').value
    return dz
