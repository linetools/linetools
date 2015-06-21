"""
Module for parsing Line List data
  Includes the Dict Definition for the Data
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp, glob

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
# def read_verner94
# def parse_morton03
# def mkvot_morton03
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
        'ion': 0,             # Ionic state
        'mol': ' '*10,            # Molecular name (H2, HD, CO, C13O)
        'Ref': ' '*50,            # Referencs
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
def read_verner94():
    # Read
    verner94 = lt_path + '/data/lines/verner94_tables.vot'
    print(
        'linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
            verner94))
    vot = vo_parse(verner94)
    # Table
    tbl_6 = vot.get_table_by_index(5).to_table(use_names_over_ids=True)

    # My table
    ldict, data = line_data(nrows=len(tbl_6))

    # Fill
    data['wrest'] = tbl_6['lambda']
    data['f'] = tbl_6['Fik']
    data['gj'] = tbl_6['Gi']
    data['gk'] = tbl_6['Gk']
    data['ion'] = tbl_6['Z'] - tbl_6['N'] + 1
    for ii,row in enumerate(tbl_6):
        data[ii]['name'] = (
            row['Species'][0:2].strip() + row['Species'][2:].strip() + 
            ' {:d}'.format(int(row['lambda'])))
        #xdb.set_trace()

    #  Finish
    data['group'] = 1
    data['Ref'] = 'Verner1994'
    data['mol'] = ''

    # Return
    return data

#
def parse_morton03(orig=False):
    '''Parse tables from Morton 2003, ApJS, 149, 205 
    '''
    # Look for VOT
    votf = lt_path + '/data/lines/morton03_table2.vot'
    morton03_tab2 = glob.glob(votf)
    if (len(morton03_tab2) > 0) & (not orig):
        print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
            morton03_tab2[0]))
        vot = vo_parse(morton03_tab2[0])
        data = vot.get_first_table().to_table(use_names_over_ids=True)
    else:

        ## Read Table 2
        morton03_tab2 = lt_path + '/data/lines/morton03_table2.dat'
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
        ionv = []
        for kk,line in enumerate(lines):
            if 'Z = ' in line:
                # Grab Z
                ipos = line.find('Z = ')
                elmZ.append(int(line[ipos+4:ipos+7]))
                ipos2 = line.find('= ')
                elmc.append(line[ipos2+2:ipos].strip())
                #xdb.set_trace()
                # Line index
                elmi.append(kk)
            if 'IP = ' in line:
                # Grab Z
                ipos = line.find(' ')
                ipos2 = line[ipos+1:].find(' ')
                ionv.append(line[ipos+1:ipos+ipos2+1].strip())
                # Line index
                ioni.append(kk)

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
                    # ion
                    gdi = np.where( (kk > np.array(ioni)) & (kk < np.roll(np.array(ioni),-1)))[0]
                    if len(gdi) != 1:
                        xdb.set_trace()
                        raise ValueError('Uh oh ion')
                    tbl[count]['ion'] = roman_to_number(ionv[gdi])
                    # Name
                    tbl[count]['name'] = elmc[gdZ]+ionv[gdi]+' {:d}'.format(
                        int(tbl[count]['wrest']))
                    #xdb.set_trace()
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
                    # Ex
                    #all_dict[count]['Ex'] = 0.  # Zero out units (for Table)

                    # Check
                    #print(tbl[count])
                    #xdb.set_trace()

                    # Increment
                    count += 1
        # Trim
        data = tbl[0:count]

        # Finish up
        #data['A'].unit = 1/u.s
        #data['gamma'].unit = 1/u.s
        #data['wrest'].unit = u.AA
        #data['Ex'].unit = 1/u.cm
        #data['Ej'].unit = 1/u.cm
        #data['Ek'].unit = 1/u.cm

        # Last
        data['group'] = 1
        data['Ref'] = 'Morton2003'
        data['mol'] = ''

    # Return
    return data

def mkvot_morton03(do_this=False, outfil=None):
    '''Used to generate a VO Table for the Morton2003 paper
    Only intended for builder usage (1.5Mb file; gzip FITS is 119kb)
    '''
    if not do_this:
        print('mkvot_morton03: It is very unlikely you want to do this')
        print('mkvot_morton03: Returning...')
        return
    # Read Morton2003
    m03 = parse_morton03()

    # Write
    if outfil is None:
        outfil = lt_path + '/data/lines/morton03_table2.vot'
    m03.write(outfil, format='votable')


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
    r_to_n = dict(I=0, II=1, III=2, IV=3, V=4, VI=5, 
        VII=6, VIII=7, IX=8, X=9)

    return r_to_n[val.strip()]



###########################################
