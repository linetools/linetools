"""
Module for parsing Line List data
  Includes the Dict Definition for the Data
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp

from astropy import units as u
from astropy import constants as const
from astropy.io import fits
from astropy.table import QTable, Column, Table

from xastropy.xutils import xdebug as xdb
lt_path = imp.find_module('linetools')[1]

#
def line_data():
    ''' Defines the dict for spectral line Data

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
        'name': '',           # Name
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
        'mol': '',            # Molecular name (H2, HD, CO, C13O)
        'Ref': '',            # References
        'group': 0            # Flag for grouping
        }

    return ldict

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
def parse_morton03():
    '''Parse tables from Morton 2003, ApJS, 149, 205 
    '''
    ## Read Table 2
    morton03_tab2 = lt_path + '/data/lines/morton03_table2.dat'
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
            elmc.append(line[ipos-5:ipos-1].strip())
            # Line index
            elmi.append(kk)
        if 'IP = ' in line:
            # Grab Z
            ipos = line.find(' ')
            ionv.append(line[ipos+1:ipos+5].strip())
            # Line index
            ioni.append(kk)

    ## Initialize Dicts
    ldict = line_data()
    ldict['Ref'] = 'Morton2003'
    all_dict = [ldict]*len(lines)

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
                all_dict[count]['wrest'] = float(line[19:28]) * u.AA
                # Z
                gdZ = np.where( (kk > np.array(elmi)) & (kk < np.roll(np.array(elmi),1)))[0]
                if len(gdZ) != 1:
                    raise ValueError('Uh oh')
                all_dict[count]['Z'] = elmZ[gdZ]
                # ion
                gdi = np.where( (kk > np.array(ioni)) & (kk < np.roll(np.array(ioni),1)))[0]
                if len(gdi) != 1:
                    raise ValueError('Uh oh')
                all_dict[count]['ion'] = roman_to_number(ionv[gdi])
                # Name
                all_dict[count]['name'] = elmc[gdZ]+ionv[gdi]+' {:d}'.format(
                    int(all_dict[count]['wrest'].value))
                # f
                all_dict[count]['f'] = float(line[79:89])
                # Ej, Ek
                all_dict[count]['Ej'] = float(line[29:38]) / u.cm
                all_dict[count]['Ek'] = float(line[40:50]) / u.cm
                # A
                all_dict[count]['A'] = float(line[59:68]) / u.s
                # gamma
                try:
                    all_dict[count]['gamma'] = float(line[69:79]) / u.s
                except ValueError:
                    pass
                # gl, gu
                all_dict[count]['gj'] = int(line[52:54])
                all_dict[count]['gk'] = int(line[56:58])           

                # Check
                print(all_dict[count])
                xdb.set_trace()

                # Increment
                count += 1
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
    if val.strip() == 'I':
        return 1
    elif val.strip() == 'II':
        return 2
    elif val.strip() == 'III':
        return 3
    elif val.strip() == 'IV':
        return 4
    elif val.strip() == 'V':
        return 5
    elif val.strip() == 'VI':
        return 6
    else:
        raise ValueError('Not setup for {:s}'.format(val)) 

###########################################
