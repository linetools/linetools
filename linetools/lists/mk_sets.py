""" Used to make set lists. Only intended for developer use.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp, glob

from astropy import units as u
from astropy.units.quantity import Quantity
from astropy import constants as const
from astropy.io import fits, ascii
from astropy.table import QTable, Column, Table


lt_path = imp.find_module('linetools')[1]
xa_path = imp.find_module('xastropy')[1]

#
def mk_ism(outfil=None, overwrite=False):
    ''' Make the ISM list from grb.lst in xastropy

    SHOULD BE DEPRECATED AFTER v0.0

    Parameters
    ----------
    outfil: str, optional
      Outfil
    '''
    raise ValueError('BAD IDEA TO RUN THIS AGAIN')
    raise ValueError('REALLY')
    # Read
    fil = xa_path+'/data/spec_lines/grb.lst'
    data = ascii.read(fil, format='fixed_width_no_header',data_start=1,
                names=('wrest', 'name', 'fval'),
                col_starts=(0,9,22), col_ends=(8,20,32))
    # Write?
    if outfil is None:
        outfil = lt_path+'/lists/sets/llist_v0.0.ascii'
    chkf = glob.glob(outfil)
    if (len(chkf) > 0) & (not overwrite):
        print('Not over-writing {:s}'.format(outfil))
        return

    # Add columns
    nrow = len(data)
    cism = Column([1]*nrow, name='fISM')
    oflgs = ['fSI', 'fHI', 'fH2', 'fCO', 'fEUV',
        'fgE', 'fgA', 'fAGN']
    allc = [cism]
    for oflg in oflgs:
       allc.append(Column([0]*nrow, name=oflg))
    data.add_columns(allc)

    # Write
    data[['wrest','name','fISM']+oflgs].write(outfil, format='ascii.fixed_width')
    print('mk_ism: Wrote {:s}'.format(outfil))

def mk_strong(infil=None, outfil=None):
    ''' Make the Strong ISM list from lls.lst in xastropy

    SHOULD BE DEPRECATED AFTER v0.0

    Parameters
    ----------
    infil : str, optional
      Starting file.  Should use latest llist_vX.X.ascii
    outfil : str, optional
      Outfil
    '''
    if infil is None:
        fils = glob.glob(lt_path+'/lists/sets/llist_v*')
        infil = fils[-1] # Should grab the lateset

    # Read
    data = ascii.read(infil, format='fixed_width')

    # Use lls.lst from xastropy    
    fil = xa_path+'/data/spec_lines/lls.lst'
    lls = ascii.read(fil, format='fixed_width_no_header',data_start=1,
                names=('wrest', 'name', 'fval'),
                col_starts=(0,9,22), col_ends=(8,20,32))

    # Set flag
    for row in lls:
        mt = np.where(np.abs(row['wrest']-data['wrest']) < 1e-3)[0]
        if len(mt) == 1:
            data[mt[0]]['fSI'] = 1
        elif len(mt) == 0:
            print('No line in Table with wrest={:g}'.format(row['wrest']))
        else:
            raise ValueError('Multiple lines in Table {:g}'.format(row['wrest']))

    # Write
    if outfil is None:
        outfil = infil
    data.write(outfil, format='ascii.fixed_width')
    print('mk_strong: Wrote {:s}'.format(outfil))

def mk_hi(infil=None, outfil=None):
    ''' Make the HI list ISM + HI

    Parameters
    ----------
    infil : str, optional
      Starting file.  Should use latest llist_vX.X.ascii
    outfil : str, optional
      Output file.  Default is to use infil
    '''
    if infil is None:
        fils = glob.glob(lt_path+'/lists/sets/llist_v*')
        fils.sort()
        infil = fils[srt][-1] # Should grab the lateset

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

    import pdb
    pdb.set_trace()
    data.write(outfil, format='ascii.fixed_width')

