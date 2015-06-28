"""
Module for parsing Line List data
  Includes the Dict Definition for the Data
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp, glob, pdb
import subprocess

from astropy import units as u
from astropy.units.quantity import Quantity
from astropy import constants as const
from astropy.io import fits, ascii
from astropy.table import QTable, Column, Table
from astropy.io.votable import parse as vo_parse

lt_path = imp.find_module('linetools')[1]

# def line_data
# def read_sets
# def read_H2
# def read_CO
# def read_verner94
# def parse_morton00
# def parse_morton03
# def mktab_morton03
# def update_fval
# def roman_to_number

# TODO
# Ingest CO data
# Ingest Galaxy lines
# Ingest AGN lines


#
def line_data(nrows=1):
    ''' Defines the dict (and/or Table) for spectral line Data
    Parameters:
    ----------
    nrows: int, optional
      Number of rows in Table [default = 1]

    Group definition:
    -----------------
        0: None
        1: "All" ISM (intended to be all atomic lines ever observed)
        2: Strong ISM
        4: HI Lyman series
        8: H2
       16: CO
       32: EUV
       64: Galaxy Emission
      128: Galaxy Absorption
      256: AGN
      512: ??
     1024: User1 (Reserved)
     2048: User2 (Reserved)
    '''
    ldict = {
        'name': ' '*20,           # Name
        'wrest': 0.*u.AA,     # Rest Wavelength (Quantity)
        'f':  0.,             # Oscillator strength
        'gk': 0.,             # Degeneracy of the upper level
        'Ej': 0./u.cm,        # Energy of lower level (relative to ground state)
        'Ek': 0./u.cm,        # Energy of upper level (relative to ground state)
        'Ex': 0./u.cm,        # Excitation energy (cm^-1)
        'A': 0./u.s,          # Einstein coefficient
        'gj': 0,              # Lower statistical weight (2J+1)
        'gk': 0,              # Upper statistical weight (2J+1)
        'gamma': 0./u.s,      # Sum of A
        'nj': 0,              # Orbital level of lower state (or vibrational level)
        'nk': 0,              # Orbital level of upper state (or vibrational level)
        'Jj': 0.,             # Tot ang mom (z projection) of lower state (or rotation level)
        'Jk': 0.,             # Tot ang mom (z projection) of upper state (or rotation level)
        'el': 0,              # Electronic transition (2=Lyman (B-X), 3=Werner (C-X)) 
        'Z': 0,               # Atomic number (for atoms)
        'Am': 0,              # Mass number (often written as "A"; only used for D)
        'ion': 0,             # Ionic state (1=Neutral)
        'mol': ' '*10,        # Molecular name (H2, HD, CO, C13O)
        'Ref': ' '*50,        # References
        'group': 0            # Flag for grouping
        }

    # Table
    clms = []
    for key in ldict.keys():
        if type(ldict[key]) is Quantity:
            clm = Column( ([ldict[key].value]*nrows)*ldict[key].unit, name=key)
        else:
            clm = Column( [ldict[key]]*nrows, name=key)
        # Append
        clms.append(clm)
    tbl = Table(clms)

    return ldict, tbl

#   
def read_sets(infil=None):
    ''' Read sets file
    Parameters:
    ---------
    infil: str, optional
      Set file
    '''
    if infil is None:
        fils = glob.glob(lt_path+'/lists/sets/llist_v*')
        infil = fils[-1] # Should grab the lateset
    # Read
    set_data = ascii.read(infil, format='fixed_width')

    # Return
    return set_data

#
def read_H2():
    ''' Simple def to read H2 data

    Returns:
    --------
    QTable of H2 lines
    '''
    H2_fil = lt_path + '/data/lines/H2_resonance.ascii'
    print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(H2_fil))
    data = Table.read(H2_fil, format='ascii', guess=False, comment=';')

    # Units
    data['wrest'].unit = u.AA

    # Rename some columns
    data.rename_column('Jp', 'Jj')
    data.rename_column('Jpp', 'Jk')
    data.rename_column('np', 'nj')
    data.rename_column('npp', 'nk')

    # Molecule column
    cmol = Column(['H2']*len(data), name='mol')
    data.add_column(cmol)

    # Group
    cgroup = Column(np.ones(len(data),dtype='int')*(2**3), name='group')
    data.add_column(cgroup)

    # Return
    return data

#
def read_CO():
    ''' Simple def to read CO UV data

    Returns:
    --------
    Table of CO lines
    '''
    CO_fil = lt_path + '/data/lines/CO_UV.ascii'
    print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(CO_fil))
    data = ascii.read(CO_fil)

    # Units

    # Rename some columns
    data.rename_column('Jp', 'Jj')
    data.rename_column('Jpp', 'Jk')
    data.rename_column('np', 'nj')
    data.rename_column('npp', 'nk')
    data.rename_column('iso', 'Am') # Isotope
    data.rename_column('wave', 'wrest') 

    data['wrest'].unit = u.AA

    # Fvalues
    data['fv'] = 10.**data['fv']
    data.rename_column('fv', 'f')

    # Molecule column
    cmol = Column(['CO']*len(data), name='mol')
    data.add_column(cmol)

    # Group
    cgroup = Column(np.ones(len(data),dtype='int')*(2**4), name='group')
    data.add_column(cgroup)

    # Return
    return data


#
def read_verner94():
    '''
    Read Verner1994 Table
    '''
    # Read
    verner94 = lt_path + '/data/lines/verner94_tab6.fits'
    print(
        'linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
            verner94))
    tbl_6 = QTable.read(verner94)


    # My table
    ldict, data = line_data(nrows=len(tbl_6))

    # Fill
    data['wrest'] = tbl_6['lambda']
    data['f'] = tbl_6['Fik']
    data['gj'] = tbl_6['Gi']
    data['gk'] = tbl_6['Gk']
    data['Z'] = tbl_6['Z']
    data['ion'] = tbl_6['Z'] - tbl_6['N'] + 1
    for ii,row in enumerate(tbl_6):
        data[ii]['name'] = (
            row['Species'][0:2].strip() + row['Species'][2:].strip() + 
            ' {:d}'.format(int(row['lambda'].value)))
        #xdb.set_trace()

    #  Finish
    data['group'] = 1
    data['Ref'] = 'Verner1994'
    data['mol'] = ''

    # Return
    return data

def parse_morton00(orig=False):
    '''Parse tables from Morton 2000, ApJS, 130, 403

    Parameters:
    -----------
    orig:
      Use original code to parse the ASCII file

    Returns:
    -----------
    data:  Table
      Atomic data
    '''
    # Look for FITS
    fitsf = lt_path + '/data/lines/morton00_table2.fits.gz'
    morton00_tab2 = glob.glob(fitsf)

    if (len(morton00_tab2) > 0) & (not orig):
        print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
            morton00_tab2[0]))
        data = Table.read(morton00_tab2[0])
    else:
        # File
        morton00_tab2 = lt_path + '/data/lines/morton00_table2.dat'
        # Call
        data = parse_morton03(orig=True, tab_fil=morton00_tab2)
        # Update
        data['Ref'] = 'Morton2000'

    # Return
    return data

#
def parse_morton03(orig=False, tab_fil=None, HIcombine=True):
    '''Parse tables from Morton 2003, ApJS, 149, 205 

    Parameters:
    -----------
    orig:
      Use original code to parse the ASCII file
    tab_fil: str, optional
      Filename to use.  Default = /data/lines/morton03_table2.dat 
    HIcombine: bool, optional
      Combine doublet for HI [True]

    Returns:
    -----------
    data:  Table
      Atomic data
    '''
    # Look for FITS
    fitsf = lt_path + '/data/lines/morton03_table2.fits.gz'
    morton03_tab2 = glob.glob(fitsf)

    if (len(morton03_tab2) > 0) & (not orig):
        print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
            morton03_tab2[0]))
        data = Table.read(morton03_tab2[0])
    else:

        ## Read Table 2
        if tab_fil is None:
            morton03_tab2 = lt_path + '/data/lines/morton03_table2.dat'
        else:
            morton03_tab2 = tab_fil
        print(
            'linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
                morton03_tab2))
        f = open(morton03_tab2, 'r')
        lines = f.readlines()
        f.close()

        ## Find Elements and Ions
        elmi = []
        elmZ = []
        elmc = []
        ioni = []
        isoi = []
        ionv = []
        for kk,line in enumerate(lines):
            #print('kk = {:d}'.format(kk))
            try:
                tmp = ('Z = ' in line) & ('A = ' in line)  # Deals with bad Byte in Morton00
            except UnicodeDecodeError:
                tmp = False
            if tmp:
                # Grab Z
                ipos = line.find('Z = ')
                elmZ.append(int(line[ipos+4:ipos+7]))
                ipos2 = line.find('= ')
                elmc.append(line[ipos2+2:ipos].strip())
                #xdb.set_trace()
                # Line index
                elmi.append(kk)

            # ISOTOPE and ION
            try:
                tmp2 = ( (('I ' in line[0:13]) | ('V ' in line[0:13]))
                    & (line[0:3] not in ['IOD','VAN']) & (line[0:2] != 'I ') )
            except UnicodeDecodeError:
                tmp2 = False
            if tmp2:
                # Grab ion
                ipos = line[0:10].find(' ')
                if ipos > 4:
                    ipos3 = line[0:10].find('I')
                    iionv = line[ipos3:ipos]
                else:
                    iionv = line[ipos:6].strip()
                if (len(iionv) == 0) | (iionv == '5s') | (iionv == 'B I'):
                    pdb.set_trace()
                ionv.append(iionv)
                if iionv == 'Z =':
                    pdb.set_trace()

                # Line index
                ioni.append(kk)

                # Deal with Isotope
                if line[0] in ['0','1','2','3','4','5','6','7','8','9']:
                    isoi.append(kk)
                # Deuterium
                if line[0] == 'D':
                    Dline = kk

        #pdb.set_trace()
        ## Initialize table
        ldict, tbl = line_data(nrows=len(lines))

        ## Parse lines with UV rest wavelength
        count = 0
        for kk,line in enumerate(lines):
            try:
                tmp = line[23] == '.'
            except IndexError:
                pass
            else:
                if tmp: # UV wavelength?
                    # Parse

                    # Ion/Isotope
                    if kk > np.max(ioni):
                        gdi = len(ioni)-1
                    else:
                        gdi = np.where( (kk > np.array(ioni)) & (kk < np.roll(np.array(ioni),-1)))[0]
                        if len(gdi) != 1:
                            pdb.set_trace()
                            raise ValueError('Uh oh ion')
                    if ioni[gdi] in isoi: # Isotope
                        continue
                    # Isotope (Atomic number)
                    if ioni[gdi] == Dline:
                        tbl[count]['Am'] = 2
                    # Ion
                    tbl[count]['ion'] = roman_to_number(ionv[gdi])

                    # Wavelength
                    tbl[count]['wrest'] = float(line[19:28]) #* u.AA
                    # Z
                    gdZ = np.where( (kk > np.array(elmi)) & (kk < np.roll(np.array(elmi),-1)))[0]
                    if len(gdZ) != 1:
                        if kk > np.max(elmi):
                            gdZ = len(elmi)-1
                        else:
                            #xdb.set_trace()
                            raise ValueError('Uh oh elm')
                    tbl[count]['Z'] = elmZ[gdZ]
                    # Name
                    tbl[count]['name'] = elmc[gdZ]+ionv[gdi]+' {:d}'.format(
                        int(tbl[count]['wrest']))
                    # f
                    try:
                        tbl[count]['f'] = float(line[79:89])
                    except ValueError:
                        continue # Skip ones without f-value
                    # Ej, Ek
                    tbl[count]['Ej'] = float(line[29:38]) #/ u.cm
                    tbl[count]['Ek'] = float(line[40:50]) #/ u.cm
                    # A
                    try:
                        tbl[count]['A'] = float(line[59:68]) #/ u.s
                    except ValueError:
                        pass
                    # gamma
                    try:
                        tbl[count]['gamma'] = float(line[69:79]) #/ u.s
                    except ValueError:
                        pass
                    # gl, gu
                    tbl[count]['gj'] = int(line[52:54])
                    tbl[count]['gk'] = int(line[56:58])           

                    # Only use combined HI lines
                    if HIcombine:
                        if ((tbl[count]['Z'] == 1) & (tbl[count]['ion']==1) 
                            & (tbl[count]['gk'] != 6)):
                            #print('Skipping HI line {:g}'.format(tbl[count]['wrest']))
                            continue 
                    # Ex
                    #all_dict[count]['Ex'] = 0.  # Zero out units (for Table)

                    # Increment
                    count += 1
        # Trim
        data = tbl[0:count]

        # Last
        data['group'] = 1
        data['Ref'] = 'Morton2003'
        data['mol'] = ''

    # Return
    return data

def mktab_morton03(do_this=False, outfil=None, fits=True):
    '''Used to generate a VO or FITS Table for the Morton2003 paper
    Only intended for builder usage (1.5Mb file; gzip FITS is 119kb)
    do_this: bool, optional
      Set to True to actually do this. Default=False
    outfil: str, optional
      Name of output file.  Defaults to a given value
    fits:  bool, optional
      Generate a FITS file?  Default=True
    '''
    if not do_this:
        print('mktab_morton03: It is very unlikely you want to do this')
        print('mktab_morton03: Returning...')
        return

    # Read Morton2003 ASCII file
    m03 = parse_morton03(orig=True)

    # Write
    if fits:
        if outfil is None:
            outfil = lt_path + '/data/lines/morton03_table2.fits'
        m03.write(outfil,overwrite=True)
    else:
        if outfil is None:
            outfil = lt_path + '/data/lines/morton03_table2.vot'
        m03.write(outfil, format='votable')
    print('mktab_morton03: Wrote {:s}'.format(outfil))
    #
    print('mktab_morton03: Now compressing...')
    subprocess.call(['gzip', '-f', outfil])

def mktab_morton00(do_this=False, outfil=None):
    '''Used to generate a FITS Table for the Morton2000 paper
    Only intended for builder usage 
    do_this: bool, optional
      Set to True to actually do this. Default=False
    outfil: str, optional
      Name of output file.  Defaults to a given value
    '''
    if not do_this:
        print('mktab_morton00: It is very unlikely you want to do this')
        print('mktab_morton00: Returning...')
        return

    # Read Morton2003
    m00 = parse_morton00(orig=True)

    # Write
    if outfil is None:
        outfil = lt_path + '/data/lines/morton00_table2.fits'
    m00.write(outfil,overwrite=True)
    print('mktab_morton00: Wrote {:s}'.format(outfil))
    #
    print('mktab_morton03: Now compressing...')
    subprocess.call(['gzip', '-f', outfil])

def update_fval(table, verbose=False):
    '''Update f-values only
    Parameters:
    -----------
    table: QTable
      Data to be updated
    '''
    # Howk 2000
    howk00_fil = lt_path + '/data/lines/howk00_table1.ascii'
    howk00 = ascii.read(howk00_fil)

    # Dress up
    howk00['wrest'].unit = u.AA

    fval = []
    fsig = []
    for row in howk00:
        ipos1 = row['fval_sig'].find('(')
        ipos2 = row['fval_sig'].find(')')
        # 
        fval.append(float(row['fval_sig'][0:ipos1]))
        fsig.append(float(row['fval_sig'][ipos1+1:ipos2]))
    # Add columns
    howk00.add_column(Column(np.array(fval), name='f'))
    howk00.add_column(Column(np.array(fsig), name='fsig')) # Error in last decimal

    # Now, finally, update
    for row in howk00:
        mt = np.where( (np.abs(table['wrest']-row['wrest']*u.AA) < 1e-3*u.AA) & 
            (table['Z'] == 26) & (table['ion'] == 2))[0]
        if len(mt) == 0:
            if verbose:
                print('update_fval: Line {:g} not in your table.'.format(row['wrest']))
        elif len(mt) == 1:
            table['f'][mt[0]] = row['f']
        else:
            raise ValueError('Uh oh')

def update_wrest(table, verbose=True):
    '''Update wrest values (and Ej,Ek)
    Parameters:
    -----------
    table: QTable
      Data to be updated
    '''

    '''
    # TiII line (Morton 2003 vs Weise 2001) 
    # Went back to Morton 2003.  If you go to Weise, you have
    #  to expunge the Verner94 row
    mt = np.where( (np.abs(table['wrest']-1910.9538*u.AA) < 1e-3*u.AA))[0] 
    table['wrest'][mt[0]] = 1910.938 * u.AA
    table['Ek'][mt[0]] = 52330.33 / u.cm
    '''

#
def roman_to_number(val):
    '''Convert simple Roman numerals to Arabic

    Parameters:
    -------------
    val: str or unicoce
      Roman numeral for conversion
    Returns:
    ------------
    Number
    '''
    r_to_n = dict(I=1, II=2, III=3, IV=4, V=5, VI=6, 
        VII=7, VIII=8, IX=9, X=10)
    try:
        num = r_to_n[val.strip()]
    except KeyError:
        print(val)
        pdb.set_trace()

    return num



###########################################
