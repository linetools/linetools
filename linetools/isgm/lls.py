""" Subclasses for LLS AbsSystem and AbsSurvey
"""

from __future__ import print_function, absolute_import, division, unicode_literals

import os, imp, glob
import numpy as np
import warnings

from astropy import units as u
from astropy.table import Table, Column

from linetools.spectralline import AbsLine
from linetools.analysis import absline as ltaa
from linetools.lists.linelist import LineList
from linetools.isgm.abssystem import AbsSystem, AbsSubSystem
from linetools.isgm import utils as ltiu
from linetools.isgm.abssurvey import AbslineSurvey
from linetools.abund import ions as ltai

class LLSSystem(AbsSystem):
    """
    Class for an LLS absorption system

    Attributes
    ----------
    tau_ll : float
      Opacity at the Lyman limit
    """

    @classmethod
    def from_datfile(cls, dat_file, tree=None, **kwargs):
        """ Read from dat_file (historical JXP format)

        Parameters
        ----------
        dat_file : str
          dat file
        tree : str, optional
          Path to data files
        kwargs :
          Passed to __init__

        Returns
        -------
        _datdict : dict
          Fills this attribute
        """
        if tree is None:
            tree = ''
        # Read datfile
        datdict = ltiu.read_dat_file(tree+dat_file)
        # Parse
        coord, zabs, name, NHI, sigNHI, clm_fil = ltiu.parse_datdict(datdict)
        kwargs['NHI'] = NHI
        kwargs['sig_NHI'] = sigNHI
        # Generate with type
        vlim = None
        slf = cls(coord, zabs, vlim, **kwargs)

        # Fill files
        slf.tree = tree
        slf.dat_file = slf.tree+dat_file

        # Parse datdict
        #   Includes Sub systems
        slf._datdict = datdict
        slf.parse_dat_file()

        return slf

    def __init__(self, radec, zabs, vlim, **kwargs):
        """Standard init

        NHI keyword is required

        Parameters
        ----------
        radec : tuple or coordinate
            RA/Dec of the sightline or astropy.coordinate
        zabs : float
          Absorption redshift
        vlim : Quantity array (2)
          Velocity limits of the system
          Defaulted to +/- 500 km/s if None (see Prochaska et al. 2016 HDLLS)
        NHI= : float, required despite being a keyword
          log10 of HI column density
        **kwargs : keywords
          passed to AbsSystem.__init__
        """
        # NHI
        try:
            NHI = kwargs['NHI']
        except KeyError:
            raise ValueError("NHI must be specified for LLSSystem")
        else:
            kwargs.pop('NHI')
        # vlim
        if vlim is None:
            vlim = [-500.,500.]*u.km/u.s
        # Generate with type
        AbsSystem.__init__(self, 'LLS', radec, zabs, vlim, NHI=NHI, **kwargs)

        # Set tau_LL
        self.tau_LL = (10.**self.NHI)*6.3391597e-18  # Should replace with photocross

        # Other
        self.zpeak = None  # Optical depth weighted redshift
        self.MH = 0.

        # Subsystems
        self.nsub = 0
        self.subsys = {}

    def parse_dat_file(self, vlim=[-300.,300]*u.km/u.s):
        """ Parse the datdict read from the .dat file

        Parameters
        ----------
        vlim : Quantity array (2), optional
          Velocity limits of the subsystems
          Should be pulled from the .clm files
        """

        # LLS keys
        self.bgsrc = self._datdict['QSO name']
        self.zem = float(self._datdict['QSO zem'])  # Was zqso
        self.MH = float(self._datdict['[M/H] ave'])
        self.nsub = int(self._datdict['N subsys'])
        self.cldyfil = self._datdict['Cloudy Grid File']

        # LLS Subsystems
        if self.nsub > 0:
            lbls= map(chr, range(65, 91))
            # Dict
            keys = (['zabs','NHI','NHIsig','NH','NHsig','log x','sigx','b','bsig','Abund file',
                     'U','Usig','flg_low','flg_alpha','[alpha/H]','sig[a/H]',
                     'flg_Fe','[Fe/H]','sig[Fe/H]','VPFIT file'])
            att = (['zabs','NHI','NHIsig','NH','NHsig','logx','sigx','bval','bsig','clm_file',
                     'U','Usig','flg_low','flg_alpha','alpha_H','sig_a_H',
                     'flg_Fe','Fe_H','sig_Fe_H','VPFIT_file'])
            values = ([0., 0., np.zeros(2), 0., np.zeros(2), 0., np.zeros(2), 0., 0.,
                    '', 0., np.zeros(2), 0, 0, 0., 0., 0, 0., 0., ''])
            null_dict = dict(zip(keys,values))
            # Loop on subsystems
            for i in range(self.nsub):
                # Generate
                zabs = float(self._datdict[lbls[i] + ' zabs'])
                self.subsys[lbls[i]] = AbsSubSystem(self, zabs, vlim, lbls[i])
                self.subsys[lbls[i]]._datdict = {}
                # Fill in dict
                for ii, key in enumerate(keys):
                    try:
                        tmpc = self._datdict[lbls[i]+' '+key]
                    except:
                        raise ValueError('lls_utils: Key "{:s}" not found in {:s}'
                                         .format(lbls[i]+key,self.dat_file))
                    else:  # Convert
                        val = null_dict[key]
                        if val.__class__ == np.ndarray:
                            self.subsys[lbls[i]]._datdict[att[ii]] = np.array(map(float,tmpc.split()))
                        else:  # Single value
                            self.subsys[lbls[i]]._datdict[att[ii]] = (map(type(val),[tmpc]))[0]
                # Set a few special ones as attributes
                self.subsys[lbls[i]].NHI = self.subsys[lbls[i]]._datdict['NHI']
                self.subsys[lbls[i]].sig_NHI = self.subsys[lbls[i]]._datdict['NHIsig']

    def get_ions(self, use_clmfile=False, idict=None, update_zvlim=True, linelist=None):
        """Parse the ions for each Subsystem

        And put them together for the full system
        Fills ._ionN with a QTable

        Parameters
        ----------
        idict : dict, optional
          dict containing the IonClms info
        use_clmfile : bool, optional
          Parse ions from a .clm file (JXP historical)
        update_zvlim : bool, optional
          Update zvlim from lines in .clm (as applicable)
        linelist : LineList
        """
        if idict is not None:
            # Manipulate for astropy Table
            #  Could probably use add_row or dict instantiation
            table = None
            for ion in idict.keys():
                Zion = ltai.name_ion(ion)
                if table is None:
                    tkeys = idict[ion].keys()
                    lst = [[idict[ion][tkey]] for tkey in tkeys]
                    table = Table(lst, names=tkeys)
                    # Extra columns
                    if 'Z' not in tkeys:
                        table.add_column(Column([Zion[0]], name='Z'))
                        table.add_column(Column([Zion[1]], name='ion'))
                else:
                    tdict = idict[ion]
                    tkeys = idict[ion].keys()
                    if 'Z' not in tkeys:
                        tdict['Z'] = Zion[0]
                        tdict['ion'] = Zion[1]
                    # Add
                    table.add_row(tdict)
            # Finish
            try:  # Historical keys
                table.rename_column('clm', 'logN')
            except:
                pass
            else:
                table.rename_column('sig_clm', 'sig_logN')
                table.rename_column('flg_clm', 'flag_N')
            self._ionN = table
        elif use_clmfile:
            # Subsystems
            if self.nsub > 0:  # This speeds things up (but is rarely used)
                linelist = LineList('ISM')
            for lbl in self.subsys.keys():
                clm_fil = self.tree+self.subsys[lbl]._datdict['clm_file']
                # Parse .clm file
                self.subsys[lbl]._clmdict = ltiu.read_clmfile(clm_fil, linelist=linelist)
                # Build components from lines
                abslines = []
                vmin,vmax = 9999., -9999.
                for wrest in self.subsys[lbl]._clmdict['lines']:
                    vmin = min(vmin, self.subsys[lbl]._clmdict['lines'][wrest].analy['vlim'][0].value)
                    vmax = max(vmax, self.subsys[lbl]._clmdict['lines'][wrest].analy['vlim'][1].value)
                    self.subsys[lbl]._clmdict['lines'][wrest].attrib['coord'] = self.coord
                    abslines.append(self.subsys[lbl]._clmdict['lines'][wrest])
                components = ltiu.build_components_from_abslines(abslines)
                # Update z, vlim
                if update_zvlim:
                    self.subsys[lbl].zabs = self.subsys[lbl]._clmdict['zsys']
                    self.subsys[lbl].vlim = [vmin, vmax]*u.km/u.s
                # Read .ion file and fill in components
                ion_fil = self.tree+self.subsys[lbl]._clmdict['ion_fil']
                self.subsys[lbl]._indiv_ionclms = ltiu.read_ion_file(ion_fil, components)
                # Parse .all file
                all_file = ion_fil.split('.ion')[0]+'.all'
                self.subsys[lbl].all_file=all_file #MF: useful to have
                _ = ltiu.read_all_file(all_file, components=components)
                # Build table
                self.subsys[lbl]._ionN = ltiu.iontable_from_components(components,ztbl=self.subsys[lbl].zabs)
                # Add to AbsSystem
                for comp in components:
                    self.add_component(comp)

            # Combine
            if self.nsub == 1:
                self._ionN = self.subsys['A']._ionN
                self._clmdict = self.subsys['A']._clmdict
                #xdb.set_trace()
            elif self.nsub == 0:
                raise ValueError('lls_utils.get_ions: Cannot have 0 subsystems..')
            else:
                self._ionN = self.subsys['A']._ionN
                self._clmdict = self.subsys['A']._clmdict
                warnings.warn('lls_utils.get_ions: Need to update multiple subsystems!! Taking A.')
        else:
            raise ValueError("Need an option in get_ions")

    def fill_lls_lines(self, bval=20.*u.km/u.s):
        """
        Generate an HI line list for an LLS.
        Goes into self.lls_lines 

        Parameters
        ----------
        bval : float, optional
          Doppler parameter in km/s
        """
        from linetools.lists import linelist as lll

        # May be replaced by component class (as NT desires)
        HIlines = lll.LineList('HI')

        self.lls_lines = []
        for lline in HIlines._data:
            aline = AbsLine(lline['wrest'], linelist=HIlines)
            # Attributes
            aline.attrib['N'] = 10**self.NHI / u.cm**2
            aline.attrib['b'] = bval
            aline.attrib['z'] = self.zabs
            # Could set RA and DEC too
            aline.attrib['RA'] = self.coord.ra
            aline.attrib['DEC'] = self.coord.dec
            self.lls_lines.append(aline)

    def flux_model(self, spec, smooth=0):
        """ Generate a LLS model given an input spectrum

        Parameters
        ----------
        spec :  Spectrum1D
        smooth : int, optional
          Number of pixels to smooth by

        Returns
        -------
        model : XSpectrum1D
          Output model is passed back as a Spectrum 
        """
        from linetools.analysis import voigt as lav

        # Energies in LLS rest-frame
        wv_rest = spec.dispersion / (self.zabs+1)
        energy = wv_rest.to(u.eV, equivalencies=u.spectral())

        # Get photo_cross and calculate tau
        tau_LL = (10.**self.NHI / u.cm**2) * ltaa.photo_cross(1,1,energy)

        # Check for lines
        if 'lls_lines' not in self.__dict__.keys():
            self.fill_lls_lines()

        tau_Lyman = lav.voigt_from_abslines(spec.dispersion, self.lls_lines, ret='tau')

        # Combine
        tau_model = tau_LL + tau_Lyman

        # Kludge around the limit
        pix_LL = np.argmin(np.fabs( wv_rest- 911.3*u.AA))
        pix_kludge = np.where((wv_rest > 911.5*u.AA) & (wv_rest < 912.8*u.AA))[0]
        tau_model[pix_kludge] = tau_model[pix_LL]
        
        # Fill in flux
        model = spec.copy()
        model.flux = np.exp(-1. * tau_model).value

        # Smooth?
        if smooth > 0:
            model.gauss_smooth(smooth)

        # Return
        return model

    def get_zpeak(self):
        """ Measure zpeak from an ionic transition
        """
        if self.ions is None:
            print('get_zpeak: Need to fill ions with get_ions first.')
            return

        # Ions for analysis
        low_ions = [ (14,2), (6,2), (13,2), (26,2), (13,3)]  # SiII,CII,AlII,FeII,AlIII
        high_ions= [(14,4), (6,4)]  # SiIV, CIV

        for tt in range(4):
            if tt == 0:
                ions = low_ions
                iflg = 1 # Standard
            elif tt == 1:
                ions = low_ions
                iflg = 2 # Saturated
            elif tt == 2:
                ions = high_ions
                iflg = 1 # Standard
            elif tt == 3:
                ions = high_ions
                iflg = 2 # Standard
            else:
                raise ValueError('Bad value')

            # Search 
            for ion in ions:
                try:
                    t = self.ions[ion]
                except KeyError:
                    continue
                # Measurement?
                if t['flag_N'] == iflg:
                # Identify the transition
                    gdi = np.where( (self.ions.trans['Z'] == ion[0]) &
                                (self.ions.trans['ion'] == ion[1]) &
                                (self.ions.trans['flag_N'] <= iflg) )[0]
                    # Take the first one
                    gdt = self.ions.trans[gdi[0]]
                    wrest = gdt['wrest']
                    flgs = self.clm_analy.clm_lines[wrest].analy['FLAGS']
                    spec_file = self.clm_analy.fits_files[flgs[1] % 64]
                    # Generate an Abs_Line with spectrum
                    line = AbsLine(wrest, z=self.clm_analy.zsys, spec_file=spec_file)
                    # vpeak
                    from linetools import utils as ltu
                    vpeak = line.vpeak()
                    self.zpeak = ltu.z_from_v(self.clm_analy.zsys, vpeak)
                    if tt == 3:
                        print('zpeak WARNING: Using saturated high-ions!!')
                    break
            else:
                continue
            break

        # Error catching
        if self.zpeak is None:
            # Skip primordial LLS
            print('lls.zpeak: No transition in {:s}'.format(self.clm_analy.clm_fil))
            return (0,0), 0.
        # Return
        return ion, vpeak

    # Output
    def __repr__(self):
        return ('[{:s}: {:s} {:s}, zabs={:g}, NHI={:g}, tau_LL={:g}, [Z/H]={:g} dex]'.format(
                self.__class__.__name__,
                 self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                 self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                 self.zabs, self.NHI, self.tau_LL, self.MH))

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'LLS'

class LLSSurvey(AbslineSurvey):
    """
    An LLS Survey class
    """

    @classmethod
    def load_HDLLS(cls):
        """ Default sample of LLS (HD-LLS, DR1)

        Return
        ------
        lls_survey
        """
        import urllib2
        import imp
        lt_path = imp.find_module('linetools')[1]
        # Pull from Internet (as necessary)
        summ_fil = glob.glob(lt_path+"/data/LLS/HD-LLS_DR1.fits")
        if len(summ_fil) > 0:
            summ_fil = summ_fil[0]
        else:
            url = 'http://www.ucolick.org/~xavier/HD-LLS/DR1/HD-LLS_DR1.fits'
            print('LLSSurvey: Grabbing summary file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            summ_fil = lt_path+"/data/HD-LLS_DR1.fits"
            with open(summ_fil, "wb") as code:
                code.write(f.read())

        # Ions
        ions_fil = glob.glob(lt_path+"/data/LLS/HD-LLS_ions.json")
        if len(ions_fil) > 0:
            ions_fil = ions_fil[0]
        else:
            url = 'http://www.ucolick.org/~xavier/HD-LLS/DR1/HD-LLS_ions.json'
            print('LLSSurvey: Grabbing JSON ion file from {:s}'.format(url))
            f = urllib2.urlopen(url)
            ions_fil = lt_path+"/data/HD-LLS_ions.json"
            with open(ions_fil, "wb") as code:
                code.write(f.read())

        # Read
        lls_survey = cls.from_sfits(summ_fil)
        # Load ions
        lls_survey.fill_ions(jfile=ions_fil)
        # Set data path (may be None)
        for lls in lls_survey._abs_sys:
            lls.spec_path = os.getenv('HDLLS_DATA')

        return lls_survey

    def __init__(self, **kwargs):
        AbslineSurvey.__init__(self, 'LLS', **kwargs)

    def cut_nhi_quality(self, sig_cut=0.4):
        """ Cut the LLS on NHI quality.

        Could put this in Absline_Survey

        Parameters
        ----------
        sig_cut : float, optional
            Limit to include as quality

        Returns
        -------
        gdNHI : ndarray
        bdNHI : ndarray
            Indices for those LLS that are good/bad
            gdNHI is a numpy array of indices
            bdNHI is a boolean array
        """
        # Cut
        gdNHI = np.where( (self.sigNHI[:, 0] < sig_cut)
                        & (self.sigNHI[:, 1] < sig_cut))[0]
        # Mask
        bdNHI = (self.NHI == self.NHI)
        bdNHI[gdNHI] = False

        # Return
        return gdNHI, bdNHI

def tau_multi_lls(wave, all_lls, **kwargs):
    """Calculate opacities on an input observed wavelength grid

    Parameters
    ----------
    wave : Quantity array
      Wavelengths
    all_lls : list
      List of LLS Class
    **kwargs : dict
      extra keywords go to lav.voigt_from_abslines

    Returns
    -------
    tau : ndarray
      Optical depth values at input wavelengths
    """
    from linetools.analysis import voigt as lav
    #
    all_tau_model = np.zeros(len(wave))
    # Loop on LLS
    for lls in all_lls:
        # LL
        wv_rest = wave / (lls.zabs+1)
        energy = wv_rest.to(u.eV, equivalencies=u.spectral())
        # Get photo_cross and calculate tau
        tau_LL = (10.**lls.NHI / u.cm**2) * ltaa.photo_cross(1,1,energy)

        # Lyman
        tau_Lyman = lav.voigt_from_abslines(wave, lls.lls_lines, ret='tau', **kwargs)
        tau_model = tau_LL + tau_Lyman

        # Kludge around the limit
        pix_LL = np.argmin(np.fabs( wv_rest- 911.3*u.AA ))
        pix_kludge = np.where((wv_rest > 911.5*u.AA) & (wv_rest < 912.8*u.AA))[0]
        tau_model[pix_kludge] = tau_model[pix_LL]
        # Add
        all_tau_model += tau_model
    # Return
    return all_tau_model

