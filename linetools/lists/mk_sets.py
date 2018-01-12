""" Used to make set lists. Only intended for developer use.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import imp, glob
import pdb
import copy

from astropy.io import ascii

from linetools.lists import parse as llp


lt_path = imp.find_module('linetools')[1]


'''
def mk_ism(outfil=None, overwrite=False):
    """ Make the ISM list from grb.lst in xastropy

    SHOULD BE DEPRECATED AFTER v0.0

    Parameters
    ----------
    outfil : str, optional
      Outfil
    """
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
'''

'''
def mk_strong(infil=None, outfil=None):
    """ Make the Strong ISM list from lls.lst in xastropy

    SHOULD BE DEPRECATED AFTER v0.0

    Parameters
    ----------
    infil : str, optional
      Starting file.  Should use latest llist_vX.X.ascii
    outfil : str, optional
      Outfil
    """
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
'''

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


def add_cashman17(outfil, infil=None, stop=True):
    """ Pre-pend Cashman+17 lines (as necessary)

    Parameters
    ----------
    outfil : str
      Output file.
    infil : str, optional
      Starting file.  Should use latest llist_vX.X.ascii
    """
    import pdb

    mt_toler = llp.llist_toler  # THIS IS IMPORTANT

    if infil is None:
        fils = glob.glob(lt_path+'/lists/sets/llist_v*')
        fils.sort()
        infil = fils[-1]  # Should grab the latest

    # Read set file
    data = ascii.read(infil, format='fixed_width')

    # Read galaxy lines (emission)
    cashman = llp.read_cashman17()

    # reference row
    tmp_row = copy.deepcopy(data[0])

    tchr = [row['name'][0:2] for row in data]
    atchr = np.array(tchr)

    nupdate = 0
    for row in cashman:
        # Add if new
        try:
            mtch = (np.abs(row['wrest'] - data['wrest']) < mt_toler) & (atchr == row['name'][0:2])
        except:
            pdb.set_trace()

        if np.sum(mtch) == 0:
            tmp_row['wrest'] = row['wrest']
            tmp_row['name'] = row['name'].replace('_', ' ')
            tmp_row['fISM'] = 1
            #print("Adding: {:s}".format(row['name']))
            tchr += [row['name'][0:2]]
            atchr = np.array(tchr)
            data.add_row(tmp_row)

        # update if already there
        elif np.sum(mtch) == 1:
            if row['name'][0:2] != 'CI':
                # get indice and replace values
                ind = np.where(mtch)[0]
                print("Updating wavelength: {:s} with {:s}, {:f} {:f}".format(
                    data['name'][ind][0], row['name'], data['wrest'][ind][0], row['wrest']))
                data['wrest'][ind] = row['wrest']
                #data['name'][ind] = row['name']
                nupdate += 1
        # this only happens if >1 line exist within 0.0001 Angstroms
        else:
            import pdb;pdb.set_trace()

    print("We updated {:d} lines".format(nupdate))

    # Add in SIV 1072 for NT
    tmp_row['wrest'] = 1072.973
    tmp_row['name'] = 'SIV 1072'
    tmp_row['fISM'] = 1
    data.add_row(tmp_row)
    tmp_row['wrest'] = 1073.518
    tmp_row['name'] = 'SIV 1073'
    tmp_row['fISM'] = 1
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
    flg += 2**0   # X-ray lines

    if flg & (2**0):
        add_cashman17('sets/llist_v1.3.ascii', infil='sets/llist_v1.2.ascii')
        #add_galaxy_lines('sets/llist_v1.3.ascii')
        #add_xray_lines('sets/llist_v1.3.ascii')
