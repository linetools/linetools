""" Module for integrating specdb within linetools.  Mainly widgets and scripts"""

try:
    import specdb
except ImportError:
    flg_specdb = False
else:
    flg_specdb = True
    from specdb.specdb import SpecDB

from astropy.coordinates import SkyCoord


def load_specb(specdb_file):
    """
    Load a specdb_file

    Parameters
    ----------
    specdb_file: str
      Full path to the specdb file

    Returns
    -------
    spdb: SpecDB or None
      Returns None if specdb cannot be imported

    """
    # Required import?
    if not flg_specdb:
        return None
    #
    spdb = SpecDB(db_file=specdb_file)
    return spdb


def load_xspec(spdb, ra, dec, group=None, **kwargs):
    """
    Load an XSpectrum1D class from a SpecDB object

    Parameters
    ----------
    spdb: SpecDB
    ra: float
     RA coordinate (ICRS)
    dec: float
     DEC coordinate (ICRS)
    group: str, optional
     Data group to grab spectrum from
    **kwargs
     Passed to XSpectrum1D, e.g. masking


    Returns
    -------
    spec: XSpectrum1D or None
      Returns None if specdb cannot be imported

    """
    # Required import?
    if not flg_specdb:
        return None

    # Need a list of group(s)
    if group is None:
        group = []
    else:
        group = [group]

    # Coord
    scoord = SkyCoord(ra=ra, dec=dec, unit='deg')

    spec, _ = spdb.spectra_from_coord(scoord, groups=group, **kwargs)

    # Return
    return spec
