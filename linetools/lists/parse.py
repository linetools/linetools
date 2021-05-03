""" Tools for parsing Line List data
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, imp, glob, pdb, gzip, sys
if not sys.version_info[0] > 2:
    import codecs
    open = codecs.open

from pkg_resources import resource_filename

from astropy import units as u
from astropy.units.quantity import Quantity
from astropy.io import fits, ascii
from astropy.table import Column, Table, vstack

from linetools.abund import roman, ions
from linetools.abund.elements import ELEMENTS

lt_path = imp.find_module('linetools')[1]


# TODO
# Ingest AGN lines
# Add Ej, Ek, Ex for emission lines (specially Balmer, Paschen and Brackett)

#
def line_data(nrows=1):
    """ Defines the dict (and/or Table) for spectral line Data

    Parameters
    ----------
    nrows : int, optional
      Number of rows in Table [default = 1]

    Notes
    -----
    Group definition:
       *    0: None
       *    1: "All" ISM (intended to be all atomic lines ever observed)
       *    2: Strong ISM
       *    4: HI Lyman series
       *    8: H2
       *   16: CO
       *   32: EUV
       *   64: Galaxy Emission
       *  128: Galaxy Absorption
       *  256: AGN
       *  512: ??
       * 1024: User1 (Reserved)
       * 2048: User2 (Reserved)
    """
    ldict = {
        'name': ' '*20,       # Name
        'wrest': 0.*u.AA,     # Rest Wavelength (Quantity)
        'f':  0.,             # Oscillator strength
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
            clm = Column( ([ldict[key].value]*nrows), name=key)
            clm.unit = ldict[key].unit
        else:
            clm = Column( [ldict[key]]*nrows, name=key)
        # Append
        clms.append(clm)

    # make it a masked Table so we can deal with Galaxy
    # emission and ISM absorption simultaneously by masking
    # out what does not make sense in one case or the other
    tbl = Table(clms, masked=True)

    return ldict, tbl


def read_sets(infil=None):
    """ Read sets file

    Parameters
    ----------
    infil : str, optional
      Set file
    """
    if infil is None:
        fils = glob.glob(lt_path+'/lists/sets/llist_v*')
        fils.sort()
        infil = fils[-1] # Should grab the lateset
    # Read
    print('read_sets: Using set file -- \n  {:s}'.format(infil))
    set_data = ascii.read(infil, format='fixed_width')

    # Return
    return set_data


def read_euv():
    """ read additional EUV lines

    Returns
    -------
    Table of EUV lines
    """
    EUV_fil = lt_path + '/data/lines/EUV_lines.ascii'
    print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(EUV_fil))
    data = Table.read(EUV_fil, format='ascii', guess=False, comment=';',delimiter='|')

    # Units
    data['wrest'].unit = u.AA

    # Return
    return data


def read_H2():
    """ Simple def to read H2 data

    Returns
    -------
    Table of H2 lines

    References
    ----------
    * Abgrall et al. 1993, A&AS, 101, 323
    * Abgrall et al. 1993, A&AS, 101, 273

    Kindly provide by co-author E. Roueff to JC Howk to JXP

    """
    H2_fil = lt_path + '/data/lines/H2_resonance.ascii'
    print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(H2_fil))
    data = Table.read(H2_fil, format='ascii', guess=False, comment=';')

    # Remove
    data.remove_column('log(w*f)')
    data.remove_column('#')

    # Rename
    data.rename_column('gu', 'gk')

    # Units
    data['wrest'].unit = u.AA
    data['gamma'].unit = 1./u.s

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

    # Reference
    data['Ref'] = 'Abgrall93'

    # Return
    return data


def read_CO():
    """ Simple def to read CO UV data

    Generated by JXP with some great effort.  See GRB 080607 paper.

    Returns
    -------
    Table of CO lines
    """
    CO_fil = lt_path + '/data/lines/CO_UV.ascii'
    print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(CO_fil))
    data = ascii.read(CO_fil)

    # Units

    # Remove column
    data.remove_column('lgwf')
    data.remove_column('id')

    # Rename some columns
    data.rename_column('Jp', 'Jj')
    data.rename_column('Jpp', 'Jk')
    data.rename_column('np', 'nj')
    data.rename_column('npp', 'nk')
    data.rename_column('iso', 'Am') # Isotope
    data.rename_column('wave', 'wrest')
    data.rename_column('gu', 'gk')
    data.rename_column('label', 'name')

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
    # Reference
    data['Ref'] = 'Prochaska09' # GRB paper

    # Return
    return data



def read_verner94():
    """ Read Verner1994 Table
    """
    # Read
    verner94 = lt_path + '/data/lines/verner94_tab6.fits'
    print(
        'linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
            verner94))
    tbl_6 = Table.read(verner94)

    # Deal with bad unit
    tbl_6['lambda'].unit = u.AA
    tbl_6 = Table(tbl_6)

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
            ' {:d}'.format(int(row['lambda'])))
    # name
    names = []
    for row in data:
        ionnm = ions.ion_to_name((row['Z'], row['ion']))
        names.append('{:s} {:d}'.format(ionnm, int(row['wrest'])))
    data['name'] = names
    #  Finish
    data['group'] = 1
    data['Ref'] = 'Verner1994'
    data['mol'] = ''

    # Return
    return data


def read_forbidden():
    """ read galaxy emission lines (forbidden)

    There may be more here: https://github.com/moustakas/impro/blob/master/pro/hiiregions/im_getmatrix.pro

    Returns
    -------
    Table of forbidden lines
    """
    forb_fil = lt_path + '/data/lines/galaxy_forbidden.ascii'
    print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(forb_fil))
    aux = Table.read(forb_fil, format='ascii')

    # My table
    ldict, data = line_data(nrows=len(aux))

    # load values using convention names
    data['wrest'] = aux['wave']
    data['wrest'].unit = u.AA
    data['Z'] = aux['Z']
    data['ion'] = aux['ion']
    for ii, row in enumerate(data):
        row['name'] = aux['name'][ii].replace('_', ' ')
    data['Ref'] = 'DESI_NIST_JM'

    # mask the galaxy data using default mask_keys
    data = mask_gal(data)

    # Return
    return data


def read_recomb():
    """ read galaxy emission lines (recombination)

    Returns
    -------
    Table of recombination lines
    """
    recomb_fil = lt_path + '/data/lines/galaxy_recomb.ascii'
    print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(recomb_fil))
    aux = Table.read(recomb_fil, format='ascii')

    # My table
    ldict, data = line_data(nrows=len(aux))

    # load values using convention names
    data['wrest'] = aux['wave']
    data['wrest'].unit = u.AA
    data['Z'] = aux['Z']
    data['ion'] = aux['ion']
    for ii, row in enumerate(data):
        row['name'] = aux[ii]['name'].replace('_', ' ')
    data['Ref'] = 'DESI_NIST_JM'

    # mask the galaxy data using default mask_keys
    data = mask_gal(data)

    # Return
    return data

def read_galabs():
    """ read galaxy absorption lines

    Returns
    -------
    Table of recombination lines
    """
    galabs_fil = lt_path + '/data/lines/galaxy_abs.ascii'
    print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(galabs_fil))
    aux = Table.read(galabs_fil, format='ascii')

    # My table -- insures proper format
    ldict, data = line_data(nrows=len(aux))

    # load values using convention names
    data['wrest'] = aux['wave']
    data['wrest'].unit = u.AA
    data['Z'] = aux['Z']
    data['ion'] = aux['ion']
    for ii, row in enumerate(data):
        row['name'] = aux[ii]['name'].replace('_', ' ')
    data['Ref'] = 'JXP_DK_Unknown'

    # mask the galaxy data using default mask_keys
    data = mask_gal(data)

    # Return
    return data


def read_cashman17():
    """ read lines from Cashman et al. 2017
    
    Returns 
    -----------
    Modified table from Cashman et al. 2017
    """
    cashman17_fil = lt_path + '/data/lines/cashman17.ascii'
    print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(cashman17_fil))
    data = Table.read(cashman17_fil, format='ascii', guess=False, comment=';',delimiter=',')

    # Remove
    data.remove_column('lower_LSJ_level')
    data.remove_column('lower_LSJ_level')
    data.remove_column('wobs')
    data.remove_column('log(g*f)')
    data.remove_column('grade')

    # Units
    data['wrest'].unit = u.AA

    # Return
    return data


def mask_gal(data, mask_keys=None):
    """Masks linelist attributes for all galaxy lines

    Parameters
    ----------
    data : Table (masked)
        The original table to mask columns for
    mask_keys : list of strings, optional
        List of column names to be masked if given
        Otherwise it uses the default:

    Returns
    -------
    data_masked : Table (masked)
        The masked version of `data`
    """
    # check input
    if not isinstance(data, (Table)):
        raise RuntimeError('The input table has to be astropy Table')

    if data.masked is not True:
        raise RuntimeError('The input Table has to be masked.')

    # set default keys to mask
    if mask_keys is None:
        mask_keys = ['A', 'el', 'nj', 'nk','group','Ek','f','mol',
                     'Ej','Am','Ex','Jj','Jk','gk','gj','gamma']

    for key in mask_keys:
        data[key].mask = True

    return data


def parse_verner96(orig=False, write=False):
    """Parse tables from Verner, Verner, & Ferland (1996, Atomic Data and Nuclear Data Tables, Vol. 64, p.1)

    Parameters
    ----------
    orig : bool, optional
      Use original code to parse the ASCII file
      Else, read from a FITS file

    Returns
    -------
    data : Table
      Atomic data
    """
    # Look for FITS
    fitsf = lt_path + '/data/lines/verner96_tab1.fits.gz'
    verner96_tab1 = glob.glob(fitsf)

    if (len(verner96_tab1) > 0) & (not orig):
        print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
            verner96_tab1[0]))
        data = Table(Table.read(verner96_tab1[0]), masked=True)
    else:
        # File
        verner96_tab1 = lt_path + '/data/lines/verner96_tab1.txt'
        # Read
        with open(verner96_tab1) as f:
            lines = f.readlines()
        # Grab the 'good' ones
        gdlines = [iline.strip() for iline in lines if len(iline.strip()) > 113]
        ldict, data = line_data(nrows=len(gdlines))
        # Loop
        for kk, line in enumerate(gdlines):
            # Z, ion
            data[kk]['Z'] = ELEMENTS[line[0:2].strip()].number
            data[kk]['ion'] = int(line[2:4].strip())
            # wrest
            data[kk]['wrest'] = float(line[47:56].strip())
            # name
            ionnm = ions.ion_to_name((data[kk]['Z'], data[kk]['ion']))
            data[kk]['name'] = '{:s} {:d}'.format(ionnm,
                    int(data[kk]['wrest']))
            # Ej, Ek
            data[kk]['Ej'] = float(line[59:73].strip())
            data[kk]['Ek'] = float(line[73:89].strip())
            # gj, gk
            data[kk]['gj'] = int(line[89:92].strip())
            data[kk]['gk'] = int(line[92:95].strip())
            # Ak
            data[kk]['A'] = float(line[95:103].strip())
            # f
            data[kk]['f'] = float(line[104:112].strip())
        # Update
        data['Ref'] = 'Verner1996'

        # Write
        if write:
            outfil = lt_path + '/data/lines/verner96_tab1.fits'
            data.write(outfil,overwrite=True)
            print('parse_verner96: Wrote {:s}'.format(outfil))
            # Compress and delete
            print('Now compressing...')
            with open(outfil) as src:
                with gzip.open(outfil+'.gz', 'wb') as dst:
                    dst.writelines(src)
            os.unlink(outfil)

    # Return
    return data


def parse_morton00(orig=False):
    """Parse tables from Morton 2000, ApJS, 130, 403

    Parameters
    ----------
    orig : bool, optional
      Use original code to parse the ASCII file

    Returns
    -------
    data : Table
      Atomic data
    """
    # Look for FITS
    fitsf = lt_path + '/data/lines/morton00_table2.fits.gz'
    morton00_tab2 = glob.glob(fitsf)

    if (len(morton00_tab2) > 0) & (not orig):
        print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
            morton00_tab2[0]))
        data = Table(Table.read(morton00_tab2[0]), masked=True)
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
    """Parse tables from Morton 2003, ApJS, 149, 205

    Parameters
    ----------
    orig : bool, optional
      Use original code to parse the ASCII file
    tab_fil : str, optional
      Filename to use.  Default = /data/lines/morton03_table2.dat 
    HIcombine : bool, optional
      Combine doublet for HI [True]

    Returns
    -------
    data : Table
      Atomic data
    """
    # Look for FITS
    fitsf = lt_path + '/data/lines/morton03_table2.fits.gz'
    morton03_tab2 = glob.glob(fitsf)

    if (len(morton03_tab2) > 0) & (not orig):
        print('linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
            morton03_tab2[0]))
        data = Table(Table.read(morton03_tab2[0]), masked=True)
    else:
        ## Read Table 2
        if tab_fil is None:
            morton03_tab2 = lt_path + '/data/lines/morton03_table2.dat'
        else:
            morton03_tab2 = tab_fil
        print(
            'linetools.lists.parse: Reading linelist --- \n   {:s}'.format(
                morton03_tab2))
        f = open(morton03_tab2, 'r', encoding="ISO-8859-1")
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
            try: # Deals with bad Byte in Morton00
                tmp = ('Z = ' in line) & ('A =' in line)  
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
                # Kludge!!
                if 'HOLMIUM' in line:
                    elmc[-1] = 'Ho'

            # ISOTOPE and ION
            try: # Deals with bad Byte in Morton00
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
                    # Skip ArI !
                    if 'Ar I' in line:
                        pass
                    else:
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
                        else:
                            gdi = gdi[0]
                    if ioni[gdi] in isoi: # Isotope
                        continue
                    # Ion
                    tbl[count]['ion'] = roman.fromRoman(ionv[gdi])

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
                    else:
                        gdZ = gdZ[0]
                    tbl[count]['Z'] = elmZ[gdZ]
                    # Name
                    try:
                        tbl[count]['name'] = elmc[gdZ]+ionv[gdi]+' {:d}'.format(
                            int(tbl['wrest'][count]))
                    except UnicodeEncodeError:
                        pdb.set_trace()
                    # Isotope (Atomic number)
                    if ioni[gdi] == Dline:
                        tbl[count]['Am'] = 2
                        tbl[count]['name'] = 'D'+ionv[gdi]+' {:d}'.format(
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
    #pdb.set_trace()
    return data

def mktab_morton03(do_this=False, outfil=None, fits=True):
    """Used to generate a VO or FITS Table for the Morton2003 paper

    Only intended for builder usage (1.5Mb file; gzip FITS is 119kb)

    Parameters
    ----------
    do_this : bool, optional
      Set to True to actually do this. Default=False
    outfil : str, optional
      Name of output file.  Defaults to a given value
    fits :  bool, optional
      Generate a FITS file?  Default=True
    """
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
        m03.write(outfil, format='votable', overwrite=True)
    print('mktab_morton03: Wrote {:s}'.format(outfil))
    # Compress and delete
    print('mktab_morton03: Now compressing...')
    with open(outfil) as src:
        with gzip.open(outfil+'.gz', 'wb') as dst:
            dst.writelines(src)
    os.unlink(outfil)


def mktab_morton00(do_this=False, outfil=None):
    """Used to generate a FITS Table for the Morton2000 paper

    Only intended for builder usage

    Parameters
    ----------
    do_this : bool, optional
      Set to True to actually do this. Default=False
    outfil : str, optional
      Name of output file.  Defaults to a given value
    """
    import subprocess
    if not do_this:
        print('mktab_morton00: It is very unlikely you want to do this')
        print('mktab_morton00: Returning...')
        return

    # Read Morton2003
    m00 = parse_morton00(orig=True)

    # Write
    if outfil is None:
        outfil = lt_path + '/data/lines/morton00_table2.fits'
    m00.write(outfil, overwrite=True)
    print('mktab_morton00: Wrote {:s}'.format(outfil))
    print('mktab_morton03: Now compressing...')
    _ = subprocess.call(['gzip', '-f', outfil])


def grab_galaxy_linelists(do_this=False):
    """ Pulls galaxy emission line lists from DESI project

    Specifically, desisim
    Writes to hard-drive
    Only run if you are building

    Parameters
    ----------
    do_this : bool, optional
      Set to True to actually do this. Default=False

    """
    if not do_this:
        print('mktab_morton00: It is very unlikely you want to do this')
        print('mktab_morton00: Returning...')
        return

    try:
        # For Python 3.0 and later
        from urllib.request import urlopen
    except ImportError:
        # Fall back to Python 2's urllib2
        from urllib2 import urlopen

    # Forbidden
    url = 'https://raw.githubusercontent.com/desihub/desisim/master/data/forbidden_lines.dat'
    f = urlopen(url)
    tab_fil = lt_path+'/data/lines/galaxy_forbidden.ascii'
    print('Writing {:s}'.format(tab_fil))
    with open(tab_fil, "wb") as code:
        code.write(f.read())

    # Recombination
    url = 'https://raw.githubusercontent.com/desihub/desisim/master/data/recombination_lines.dat'
    f = urlopen(url)
    tab_fil = lt_path+'/data/lines/galaxy_recomb.ascii'
    print('Writing {:s}'.format(tab_fil))
    with open(tab_fil, "wb") as code:
        code.write(f.read())


def update_fval(table, verbose=False):
    """Update f-values from the literature

    Primarily for modifying lines in the ISM lists (e.g. Morton2003)

    Parameters
    ----------
    table : Table
      Data to be updated
    verbose : bool, optional

    Returns
    -------
    table : Table
      Updated table.
      Note: This return is required to handle the vstack, i.e.
      as opposed to modifying the input table in place.
    """
    # Shectman et al. 1998, ApJ, 504, 921 
    #   Morton2003 cites this but uses a different f-value
    imn = np.argmin(np.abs(table['wrest']-1526.707))
    table['f'][imn] = 0.127
    table['Ref'][imn] = 'Shectman1998'

    # Howk 2000 (using Weise 2002 as in Morton for FeII 1142,1143,1144)
    howk00_fil = lt_path + '/data/lines/howk00_table1.ascii'
    howk00 = ascii.read(howk00_fil, comment='#')

    # Jenkins & Tripp (2006, ApJ, 637, 548)
    imn = np.argmin(np.abs(table['wrest']-1317.217))
    table['f'][imn] = 0.0571
    table['Ref'][imn] = 'JT2006'


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
        mt = np.where( (np.abs(table['wrest'].data*table['wrest'].unit-row['wrest']*u.AA) < 1e-3*u.AA) &
            (table['Z'] == 26) & (table['ion'] == 2))[0]
        if len(mt) == 0:
            if verbose:
                print('update_fval: Line {:g} not in your table.'.format(row['wrest']))
        else:
            table['f'][mt] = row['f']

    ## ##
    # Lines without f-value but of interest

    # AsII
    mn = np.min(np.abs(table['wrest']-1355.934))  # In Morton2000
    if mn > 0.05:  # Ang
        _, new_row = line_data()
        new_row['Z'] = 33
        new_row['ion'] = 2
        new_row['wrest'] = 1355.934
        new_row['name'] = 'AsII 1355'
        new_row['f'].mask = True
        # Stack
        table = vstack([table, new_row])

    return table


def update_gamma(table):
    """Update/add-in gamma values

    Parameters
    ----------
    table : Table
      Data to be updated
    verbose : bool, optional
    """
    # HI - Morton doesn't give these for the combined HI lines (sensible)
    #  Nor does he give them for the lines beyond Ly-d (not sure why)
    try:
        HI = np.where((table['Z']==1) & (table['ion']==1))[0]
    except KeyError: # Molecules
        pass
    else:
        if len(HI) > 0:
            # Same kludge as in atom.dat of VPFIT for higher order lines
            table['gamma'][HI] = table['A'][HI]
            # More accurate for stronger lines (pulled from Morton) [all in units of s^-1]
            gdict = {1215.670: 6.265E+08, 1025.7222: 1.897E+08, # From Morton
                972.5367: 8.127E+07, 949.7430: 4.204E+07, 937.8034: 2.450E+07}
            for key in gdict.keys():
                mt = np.where( (np.abs(table['wrest']-key) < 1e-4))[0]
                if len(mt) > 0:
                    table['gamma'][mt[0]] = gdict[key]


def update_wrest(table, verbose=True):
    """Update wrest values (and Ej,Ek)

    Parameters
    ----------
    table : Table
      Data to be updated
    verbose : bool, optional
    """
    pass
    '''
    # TiII line (Morton 2003 vs Weise 2001) 
    # Went back to Morton 2003.  If you go to Weise, you have
    #  to expunge the Verner94 row
    mt = np.where( (np.abs(table['wrest']-1910.9538*u.AA) < 1e-3*u.AA))[0] 
    table['wrest'][mt[0]] = 1910.938 * u.AA
    table['Ek'][mt[0]] = 52330.33 / u.cm
    '''
#
def load_datasets(datasets, tol=1e-3):
    """Load up all the inidividual linelist datasets

    Parameters
    ----------
    datasets : list of func
      Routines to call for generating the dataset
    Returns
    -------

    """
    flag_fval = True  # Update f-values?
    flag_wrest = True  # Update wavelengths?
    flag_gamma = True  # Update gamma values (recommended)

    full_table = None
    all_func = []
    tmp, _ = line_data(1)
    tkeys = list(tmp.keys())
    # Loop on data sets
    for func in datasets:
        # Query if read already
        if func not in all_func:
            # Read
            table = func()

            # Check keys
            for key in table.keys():
                if key not in tkeys:
                    pdb.set_trace()

            # Add extras
            if full_table is None:
                full_table = table
            else:
                '''
                # Unique values
                wrest = full_table['wrest']
                newi = []
                for jj, row in enumerate(table):
                    mt = np.abs(row['wrest'] - wrest) < tol
                    if mt.sum() == 0:
                        newi.append(jj)
                # Append
                '''
                full_table = vstack([full_table, table])
            # Save to avoid repeating
            all_func.append(func)

    # Save
    _fulltable = full_table.copy()

    # Update wavelength values
    if flag_wrest:
        update_wrest(_fulltable)

    # Update f-values (Howk00)
    if flag_fval:
        _fulltable = update_fval(_fulltable)

    # Update gamma-values (Mainly HI)
    if flag_gamma:
        update_gamma(_fulltable)

    # Finish
    return _fulltable


def _write_ref_table(outfile=None):
    """ Write a reference table enabling faster I/O for *all* line lists
    For developer use only.

    Note that after running this, you need to manually copy the table
    produced to linetools/data/lines/ISM_table.fits inside the github
    repository, and then check it in.
    """
    import warnings
    import datetime
    import getpass
    date = str(datetime.date.today().strftime('%Y-%b-%d'))
    user = getpass.getuser()
    if outfile is None:
        outfile = resource_filename('linetools', 'data/lines/linelist.ascii')

    # Define datasets: In order of Priority
    datasets = [parse_morton03, parse_morton00, parse_verner96,
                read_verner94, read_euv]  # Morton 2003, Morton 00, Verner 96, Verner 94
    datasets += [read_H2, read_CO] # H2 (Abrigail), CO (JXP)
    datasets += [read_forbidden, read_recomb, read_galabs] # Galaxy lines

    # Load
    full_table = load_datasets(datasets)

    # Check the sets file
    #set_data = read_sets()
    #pdb.set_trace()

    # Meta
    full_table.meta['Creator'] = user
    full_table.meta['CreationDate'] = date
    # Write
    warnings.warn("About to overwrite: {:s}".format(outfile))
    warnings.warn("Proceed only if you know what you are doing!")
    pdb.set_trace()
    full_table.write(outfile, format='ascii.ecsv', overwrite=True)


if __name__ == '__main__':
    # Generate tables
    #mktab_morton00(do_this=True)
    # ISM
    #_write_ref_ISM_table()
    # Full table
    _write_ref_table()

