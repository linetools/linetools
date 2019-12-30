""" Used to make set lists. Only intended for developer use.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import imp, glob
import copy

from astropy.io import ascii

from linetools.lists import parse as llp


lt_path = imp.find_module('linetools')[1]


def mk_hi(infil=None, outfil=None, stop=True):
    """ Make the HI list ISM + HI

    Parameters
    ----------
    infil : str, optional
      Starting file.  Should use latest llist_vX.X.ascii
    outfil : str, optional
      Output file.  Default is to use infil
    """
    if infil is None:
        fils = glob.glob(lt_path+'/lists/sets/llist_v*')
        fils.sort()
        infil = fils[-1] # Should grab the lateset

    # Read
    data = ascii.read(infil, format='fixed_width')

    # Parse on HI
    idx = []
    for kk,row in enumerate(data):
        if row['name'][0:3] == 'HI ':  # Necessary to avoid DI
            data[kk]['fHI'] = 1

    # Write
    if outfil is None:
        outfil = infil

    if stop:
        import pdb
        pdb.set_trace()
    data.write(outfil, format='ascii.fixed_width')


def add_galaxy_lines(outfil, infil=None, stop=True):
    """ Append galaxy lines (as necessary)

    Parameters
    ----------
    outfil : str
      Output file.
    infil : str, optional
      Starting file.  Should use latest llist_vX.X.ascii
    """
    import pdb
    if infil is None:
        fils = glob.glob(lt_path+'/lists/sets/llist_v*')
        fils.sort()
        infil = fils[-1]  # Should grab the lateset

    # Read set file
    data = ascii.read(infil, format='fixed_width')
    # reformat name column to avoid truncation of names
    data['name'] = data['name'].astype("|S20")

    # Read galaxy lines (emission)
    forbidden = llp.read_forbidden()
    recomb = llp.read_recomb()

    tmp_row = copy.deepcopy(data[0])
    for key in ['fISM', 'fSI', 'fHI', 'fEUV', 'fAGN']:
        tmp_row[key] = 0
    tmp_row['fgE'] = 1
    # Add if new
    for row in forbidden:
        if np.sum(np.abs(row['wrest']-data['wrest']) < 0.0001) == 0:
            tmp_row['wrest'] = row['wrest']
            tmp_row['name'] = row['name']
            data.add_row(tmp_row)
    # Add if new
    for row in recomb:
        if np.sum(np.abs(row['wrest']-data['wrest']) < 0.0001) == 0:
            tmp_row['wrest'] = row['wrest']
            tmp_row['name'] = row['name'].replace('_',' ')
            data.add_row(tmp_row)

    # Write
    print('Make sure you want to do this!')
    if stop:
        import pdb
        pdb.set_trace()
    data.write(outfil, format='ascii.fixed_width')


def add_xray_lines(outfil, infil=None, stop=True):
    """ Pre-pend X-ray lines (as necessary)

    Parameters
    ----------
    outfil : str
      Output file.
    infil : str, optional
      Starting file.  Should use latest llist_vX.X.ascii
    """
    if infil is None:
        fils = glob.glob(lt_path+'/lists/sets/llist_v*')
        fils.sort()
        infil = fils[-1]  # Should grab the latest

    # Read set file
    data = ascii.read(infil, format='fixed_width')

    # Read galaxy lines (emission)
    v96 = llp.parse_verner96()

    tmp_row = copy.deepcopy(data[0])
    for key in ['fISM', 'fSI', 'fHI', 'fAGN']:
        tmp_row[key] = 0
    tmp_row['fEUV'] = 1
    # Add if new
    for row in v96:
        if row['wrest'] > 100.:
            continue
        if np.min(np.abs(row['wrest']-data['wrest'])) > 0.0001:
            tmp_row['wrest'] = row['wrest']
            import pdb; pdb.set_trace()
            tmp_row['name'] = row['name']
            data.add_row(tmp_row)

    # Sort
    data.sort('wrest')

    # Write
    print('Make sure you want to do this!')
    if stop:
        import pdb; pdb.set_trace()
    data.write(outfil, format='ascii.fixed_width', overwrite=True)


# Test
if __name__ == '__main__':
    flg = 0
    #flg += 2**0   # X-ray lines
    flg += 2**1   # v1.3

    if flg & (2**0):
        add_galaxy_lines('sets/llist_v1.2.ascii')
        add_xray_lines('sets/llist_v1.2.ascii')

    # JXP -- 2019 May 02  (only added HeII 4687, but whatever)
    if flg & (2**1):
        add_galaxy_lines('sets/llist_v1.3.ascii')
        add_xray_lines('sets/llist_v1.3.ascii')
