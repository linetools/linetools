""" Class for absorption line component
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import pdb
import numpy as np
import warnings

from astropy import constants as const
from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Column

from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.analysis import absline as ltaa
from linetools.analysis import plots as ltap
from linetools.spectralline import AbsLine, SpectralLine
from linetools.abund import ions
from linetools import utils as ltu
from linetools.lists.linelist import LineList

# Global import for speed
c_kms = const.c.to('km/s').value

# Class for Components
class AbsComponent(object):
    """
    Class for an absorption component

    Attributes
    ----------
    name : str
        Name of the component, e.g. `Si II`
    coord : SkyCoord
        Sky coordinate
    Zion : tuple 
        Atomic number, ion -- (int,int)
        e.g. (8,1) for OI
        Note: (-1, -1) is special and is meant for molecules (e.g. H2)
              This notation will most likely be changed in the future.
    zcomp : float
        Component redshift
    vlim : Quantity array
        Velocity limits of the component
        e.g.  [-300,300]*u.km/u.s
    A : int
        Atomic mass -- used to distinguish isotopes
    Ej : Quantity
        Energy of lower level (1/cm)
    comment : str
        A comment, default is ``
    """
    @classmethod
    def from_abslines(cls, abslines, stars=None, **kwargs):
        """Instantiate from a list of AbsLine objects

        Parameters
        ----------
        abslines : list 
          List of AbsLine objects
        stars : str, optional
          Asterisks to append to the ion name (e.g. fine-structure, CII*)
        """
        # Check
        if not isinstance(abslines, list):
            raise IOError("Need a list of AbsLine objects")
        if not all(isinstance(x, AbsLine) for x in abslines):
            raise IOError("List needs to contain only AbsLine objects")

        # Instantiate with the first line
        init_line = abslines[0]
        slf = cls( init_line.attrib['coord'], (init_line.data['Z'],init_line.data['ion']),
                   init_line.z, init_line.limits.vlim,
                   Ej=init_line.data['Ej'], stars=stars)
        slf._abslines.append(init_line)
        # Append with component checking
        if len(abslines) > 1:
            for absline in abslines[1:]:
                slf.add_absline(absline, **kwargs)
        # Return
        return slf

    @classmethod
    def from_component(cls, component, **kwargs):
        """ Instantiate from an AbsComponent object

        Uses RA/DEC, Zion, Ej, A, z, vlim

        Parameters
        ----------
        component : AbsComponent
           An AbsComponent object

        Returns
        -------
        AbsComponent
        """
        # Check
        if not isinstance(component, AbsComponent):
            raise IOError('Need an AbsComponent object')
        # Return
        return cls(component.coord, component.Zion, component.zcomp, component.vlim, Ej=component.Ej,
                   A=component.A, name=component.name, **kwargs)

    @classmethod
    def from_dict(cls, idict, coord=None, **kwargs):
        """ Instantiate from a dict

        Parameters
        ----------
        idict : dict

        Returns
        -------

        """
        if coord is not None:
            radec = coord
        else:
            radec = SkyCoord(ra=idict['RA']*u.deg, dec=idict['DEC']*u.deg)
        # Init
        #slf = cls(radec, tuple(idict['Zion']), idict['zcomp'], Quantity(idict['vlim'], unit='km/s'),
        slf = cls(radec, tuple(idict['Zion']), idict['zcomp'], idict['vlim']*u.km/u.s,
                  Ej=idict['Ej']/u.cm, A=idict['A'],
                  Ntup = tuple([idict[key] for key in ['flag_N', 'logN', 'sig_logN']]),
                  comment=idict['comment'], name=idict['Name'])
        # Add lines
        for key in idict['lines'].keys():
            iline = SpectralLine.from_dict(idict['lines'][key], coord=coord, **kwargs)
            slf.add_absline(iline, **kwargs)
        # Return
        return slf

    def __init__(self, radec, Zion, zcomp, vlim, Ej=0./u.cm, A=None,
                 Ntup=None, comment='', name=None, stars=None):
        """  Initiator

        Parameters
        ----------
        radec : tuple or SkyCoord
            (RA,DEC) in deg or astropy.coordinate.SkyCoord
        Zion : tuple
            Atomic number, ion -- (int,int)
            e.g. (8,1) for OI
            Note: (-1, -1) is special and is meant for moleculer (e.g. H2)
                  This notation will most likely change in the future.
        zcomp : float
            Absorption component redshift
        vlim : Quantity array
            Velocity limits of the component w/r to `z`
            e.g.  [-300,300]*u.km/u.s
        A : int, optional
            Atomic mass -- used to distinguish isotopes
        Ntup : tuple
            (int,float,float)
            (flag_N,logN,sig_logN)
            flag_N : Flag describing N measurement  (0: no info; 1: detection; 2: saturated; 3: non-detection)
            logN : log10 N column density
            sig_logN : Error in log10 N
        Ej : Quantity, optional
            Energy of lower level (1/cm)
        stars : str, optional
            asterisks to add to name, e.g. '**' for CI**
            Required if name=None and Ej>0.
        comment : str, optional
            A comment, default is ``
        """

        # Required
        self.coord = ltu.radec_to_coord(radec)
        self.Zion = Zion
        self.zcomp = zcomp
        self.vlim = vlim

        # Optional
        self.A = A
        self.Ej = Ej
        self.comment = comment
        if Ntup is not None:
            self.flag_N = Ntup[0]
            self.logN = Ntup[1]
            self.sig_logN = Ntup[2]
            _, _ = ltaa.linear_clm(self)  # Set linear quantities
        else:
            self.flag_N = 0
            self.logN = 0.
            self.sig_logN = 0.

        # Name
        if (name is None) and (self.Zion != (-1, -1)):
            iname = ions.ion_name(self.Zion, nspace=0)
            if self.Ej.value > 0:  # Need to put *'s in name
                try:
                    iname += stars
                except:
                    raise IOError("Need to provide 'stars' parameter.")
            self.name = '{:s}_z{:0.5f}'.format(iname, self.zcomp)
        elif (name is None) and (self.Zion == (-1, -1)):
            self.name = 'mol_z{:0.5f}'.format(self.zcomp)
        else:
            self.name = name

        # Potential for attributes
        self.attrib = dict()

        # Other
        self._abslines = []

    def add_absline(self, absline, tol=0.1*u.arcsec, chk_vel=True,
                    chk_sep=True, vtoler=1., **kwargs):
        """Add an AbsLine object to the component if it satisfies
        all of the rules.

        For velocities, we demand that the new line has a velocity
        range that is fully encompassed by the component.

        Parameters
        ----------
        absline : AbsLine
        tol : Angle, optional
          Tolerance on matching coordinates.  Only used if chk_sep=True
        chk_vel : bool, optional
          Perform velocity test (can often be skipped)
          Insist the bounds of the AbsLine are within 1km/s of the Component
             (allows for round-off error)
        chk_sep : bool, optional
          Perform coordinate check (expensive)
        vtoler : float
          Tolerance for velocity in km/s (must be positive)
        """
        if vtoler < 0:
            raise ValueError('vtoler must be positive!')

        # Perform easy checks
        if chk_sep:
            testc = bool(self.coord.separation(absline.attrib['coord']) < tol)
        else:
            testc = True

        if self.Zion == (-1,-1):  #(-1,-1) represents molecules
            testZ = True
            testi = True
            testE = True
        else: # atoms
            testZ = self.Zion[0] == absline.data['Z']
            testi = self.Zion[1] == absline.data['ion']
            testE = bool(self.Ej == absline.data['Ej'])

        # Now redshift/velocity
        if chk_vel:
            dz_toler = (1 + self.zcomp) * vtoler / c_kms  # Avoid Quantity for speed
            zlim_line = absline.limits.zlim  # absline.z + (1 + absline.z) * absline.limits.vlim.to('km/s').value / c_kms
            zlim_comp = self.zcomp + (1+self.zcomp) * self.vlim.to('km/s').value / c_kms
            testv = (zlim_line[0] >= (zlim_comp[0] - dz_toler)) & (
                zlim_line[1] <= (zlim_comp[1] + dz_toler))
        else:
            testv = True
        # Combine
        test = testc & testZ & testi & testE & testv
        # Isotope
        if self.A is not None:
            raise ValueError('Not ready for this yet.')
        # Append?
        if test:
            self._abslines.append(absline)
        else:
            warnings.warn("Failed add_absline test")
            print('Input absline with wrest={:g} at z={:.3f} does not match component rules. Not appending'.format(absline.wrest,
                                                                                                                   absline.z))
            if not testv:
                print("Absline velocities lie beyond component\n Set chk_vel=False to skip this test.")
            if not testc:
                print("Absline coordinates do not match.  Best to set them")

    def add_abslines_from_linelist(self, llist='ISM', init_name=None, wvlim=None, min_Wr=None, **kwargs):
        """
        It adds associated AbsLines satisfying some conditions (see parameters below).

        Parameters
        ----------
        llist : str, optional
            Name of the linetools.lists.linelist.LineList
            object where to look for the transition names.
            Default is 'ISM', which means the function looks
            within `list = LineList('ISM')`.
        init_name : str, optional
            Name of the initial transition used to define the AbsComponent
        wvlims : Quantity array, optional
            Observed wavelength limits for AbsLines to be added.
            e.g. [1200, 2000]*u.AA.
        min_Wr : Quantity, optional
            Minimum rest-frame equivalent with for AbsLines to be added.
            This is calculated in the very low optical depth regime tau0<<1,
            where Wr is independent of Doppler parameter or gamma (see eq. 9.15 of
            Draine 2011). Still, a column density attribute for the AbsComponent
            is needed.

        Returns
        -------
        Adds AbsLine objects to the AbsComponent._abslines list.

        Notes
        -----
        **kwargs are passed to AbsLine.add_absline() method.

        """
        # get the transitions from LineList
        llist = LineList(llist)
        if init_name is None:  # we have to guess it
            if (self.Zion) == (-1, -1):  # molecules
                # init_name must be in self.attrib (this is a patch)
                init_name = self.attrib['init_name']
            else:  # atoms
                init_name = ions.ion_name(self.Zion, nspace=0)
        transitions = llist.all_transitions(init_name)

        # unify output to be always QTable
        if isinstance(transitions, dict):
            transitions = llist.from_dict_to_qtable(transitions)

        # check wvlims
        if wvlim is not None:
            cond = (transitions['wrest']*(1+self.zcomp) >= wvlim[0]) & \
                   (transitions['wrest']*(1+self.zcomp) <= wvlim[1])
            transitions = transitions[cond]

        # check outputs
        if len(transitions) == 0:
            warnings.warn("No transitions satisfying the criteria found. Doing nothing.")
            return

        # loop over the transitions when more than one found
        for transition in transitions:
            iline = AbsLine(transition['name'], z=self.zcomp, linelist=llist)
            iline.limits.set(self.vlim)
            iline.attrib['coord'] = self.coord
            iline.attrib['logN'] = self.logN
            iline.attrib['sig_logN'] = self.sig_logN
            iline.attrib['flag_N'] = self.flag_N
            iline.attrib['N'] = 10**iline.attrib['logN'] / (u.cm * u.cm)
            iline.attrib['sig_N'] = 10**iline.attrib['sig_logN'] / (u.cm * u.cm)

            for key in self.attrib.keys():
                iline.attrib[key] = self.attrib[key]

            if min_Wr is not None:
                # check logN is defined
                logN = self.logN
                if logN == 0:
                    warnings.warn("AbsComponent does not have logN defined. Appending AbsLines "
                                  "regardless of min_Wr.")
                else:
                    N = 10**logN / (u.cm*u.cm)
                    Wr_iline = iline.get_Wr_from_N(N=N)  # valid for the tau0<<1 regime.
                    if Wr_iline < min_Wr:  # do not append
                        continue
            # add the absline
            self.add_absline(iline)

    def build_table(self):
        """Generate an astropy QTable out of the abs lines
        Returns
        -------
        comp_tbl : QTable
        """
        if len(self._abslines) == 0:
            return
        comp_tbl = QTable()
        comp_tbl.add_column(Column([iline.wrest.to(u.AA).value for iline in self._abslines]*u.AA, name='wrest'))
        comp_tbl.add_column(Column([iline.z for iline in self._abslines], name='z'))
        for attrib in ['flag_N', 'logN', 'sig_logN']:
            comp_tbl.add_column(Column([iline.attrib[attrib] for iline in self._abslines], name=attrib))
        # Return
        return comp_tbl

    def cog(self, redo_EW=False, show_plot=False, **kwargs):
        """Perform a COG analysis on the component

        Parameters
        ----------
        redo_EW : bool, optional
          Re-analyze each line for its EW
        show_plot : bool, optional
          Generate plot and show

        Returns
        -------
        logN : float
          COG column density
        b : Quantity
          COG Doppler parameter (km/s)
        """
        from linetools.analysis import cog as ltcog
        reload(ltcog)
        # Redo EWs?
        if redo_EW:
            for aline in self._abslines:
                aline.measure_restew(**kwargs)
        # COG setup
        wrest = np.array([aline.wrest.to('AA').value for aline in self._abslines])*u.AA
        f = np.array([aline.data['f'] for aline in self._abslines])
        EW = np.array([aline.attrib['EW'].to('AA').value for aline in self._abslines])*u.AA
        sig_EW = np.array([aline.attrib['sig_EW'].to('AA').value for aline in self._abslines])*u.AA
        # COG analysis
        COG_dict = ltcog.single_cog_analysis(wrest, f, EW, sig_EW=sig_EW)
        # COG plot
        if show_plot:
            ltcog.cog_plot(COG_dict)
        # Return
        return COG_dict

    def plot_Na(self, show=True, **kwargs):
        """Plot apparent column density Na profiles
        """
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import matplotlib as mpl
        try:  # Nicer view, especially in notebook
            import seaborn as sns
            sns.set(context="notebook", font_scale=2)
        except ImportError:
            pass
        mpl.rcParams['font.family'] = 'stixgeneral'
        mpl.rcParams['font.size'] = 15.
        # Check for spec
        gdiline = []
        for iline in self._abslines:
            if isinstance(iline.analy['spec'], XSpectrum1D):
                gdiline.append(iline)
        nplt = len(gdiline)
        if nplt == 0:
            print("Load spectra into the absline.analy['spec']")
            return
        atom_cst = (const.m_e.cgs*const.c.cgs / (np.pi * (const.e.esu**2).cgs)).to(u.AA*u.s/(u.km*u.cm**2))
        # Setup plot
        plt.clf()
        ax = plt.gca()

        fw_sv = 0.*u.AA
        ymax = 0.
        for qq, iline in enumerate(gdiline):
            # Calculate
            velo = iline.analy['spec'].relative_vel((1+iline.z)*iline.wrest)
            cst = atom_cst/(iline.data['f']*iline.wrest)  # / (u.km/u.s) / u.cm * (u.AA/u.cm)
            Na = np.log(1./np.maximum(iline.analy['spec'].flux, iline.analy['spec'].sig)) * cst

            # Figure out ymnx
            pixmnx = (velo > self.vlim[0]) & (velo < self.vlim[1])
            if iline.data['f']*iline.wrest > fw_sv:
                ymax = max(np.max(Na[pixmnx].value), ymax)
                fw_sv = iline.data['f']*iline.wrest
            # Plot
            ax.plot(velo, Na, '-', linestyle='steps-mid', label=iline.data['name'])
            # ax.plot(velo, iline.analy['spec'].sig, 'r:')
        # Axes
        ax.set_xlim(self.vlim.value)
        ax.set_ylim(-0.2*ymax, 5*ymax)
        # ax.set_ylim(ymnx)
        ax.minorticks_on()
        ax.set_xlabel('Relative Velocity (km/s)')
        ax.set_ylabel(r'Apparent Column (cm$^{-2}$ per km/s)')
        # Legend
        legend = ax.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                           handletextpad=0.3, fontsize='large')

        plt.tight_layout(pad=0.2, h_pad=0., w_pad=0.1)
        if show:
            plt.show()
        plt.close()

    def reset_vlim_from_abslines(self, verbose=False):
        """ Resets the vlim value using the AbsLines

        Parameters
        ----------

        """
        for aline in self._abslines:
            if aline.analy['vlim'][0] < self.vlim[0]:
                if verbose:
                    print('Resetting vlim0 from {}'.format(aline))
                self.vlim[0] = aline.analy['vlim'][0]
            if aline.analy['vlim'][1] > self.vlim[1]:
                if verbose:
                    print('Resetting vlim1 from {}'.format(aline))
                self.vlim[1] = aline.analy['vlim'][1]

    def synthesize_colm(self, overwrite=False, redo_aodm=False, **kwargs):
        """Synthesize column density measurements of the component.
        Default is to use the current AbsLine values, but the user can
        request that those be re-calculated with AODM.

        Parameters
        ----------
        overwrite : bool, optional
          Clobber any previous measurement
        redo_aodm : bool, optional
          Redo the individual column density measurements (likely AODM)

        Returns
        -------
        None
          Fills the component attributes instead
        """
        # Check
        if (self.flag_N != 0) and (not overwrite):
            raise IOError("Column densities already set.  Use overwrite=True to redo.")
        # Redo?
        if redo_aodm:
            for aline in self._abslines:
                aline.measure_aodm(**kwargs)
        # Collate
        self.flag_N = 0
        for aline in self._abslines:
            if aline.attrib['flag_N'] == 0:  # No value
                warnings.warn("Absline {} has flag=0.  Hopefully you expected that".format(str(aline)))
                continue
            # Check N is filled
            if np.allclose(aline.attrib['N'].value, 0.):
                raise ValueError("Need to set N in attrib.  \n Consider linear_clm in linetools.analysis.absline")
            if aline.attrib['flag_N'] == 1:  # Good value?
                if self.flag_N == 1:  # Weighted mean
                    # Original
                    weight = 1. / self.sig_N**2
                    mu = self.N * weight
                    # Update
                    weight += 1./aline.attrib['sig_N']**2
                    self.N = (mu + aline.attrib['N']/aline.attrib['sig_N']**2) / weight
                    self.sig_N = np.sqrt(1./weight)
                else:  # Fill
                    self.N = aline.attrib['N']
                    self.sig_N = aline.attrib['sig_N']
                    self.flag_N = 1
            elif aline.attrib['flag_N'] == 2:  # Lower limit
                if self.flag_N in [0, 3]:
                    self.N = aline.attrib['N']
                    self.sig_N = aline.attrib['sig_N']
                    self.flag_N = 2
                elif self.flag_N == 2:
                    if aline.attrib['N'] > self.N:
                        self.N = aline.attrib['N']
                        self.sig_N = aline.attrib['sig_N']
                elif self.flag_N == 1:
                    pass
            elif aline.attrib['flag_N'] == 3:  # Upper limit
                if self.flag_N == 0:
                    self.N = aline.attrib['N']
                    self.sig_N = aline.attrib['sig_N']
                    self.flag_N = 3
                elif self.flag_N in [1, 2]:
                    pass
                elif self.flag_N == 3:
                    if aline.attrib['N'] < self.N:
                        self.N = aline.attrib['N']
                        self.sig_N = aline.attrib['sig_N']
            elif aline.attrib['flag_N'] == 0:  # No value
                warnings.warn("Absline {} has flag=0.  Hopefully you expected that")
            else:
                raise ValueError("Bad flag_N value")
        # Log values
        if self.flag_N > 0:
            self.logN, self.sig_logN = ltaa.log_clm(self)

    def repr_vpfit(self, b=10.*u.km/u.s, tie_strs=('', '', ''), fix_strs=('', '', '')):
        """
        String representation for VPFIT (line fitting software) in its fort.26 format

        Parameters
        ----------
        b : Quantity, optional
            Doppler parameter of the component. Default is 10*u.km/u.s
        tie_strs : tuple of strings, optional
            Strings to be used for tying parameters (z,b,logN),
            respectively.  These are all converted to lower case
            format, following VPFIT convention.
        fix_strs : tuple of strings, optional
            Strings to be used for fixing parameters (z,b,logN),
            respectively.  These are all converted to upper case
            format, following VPFIT convention.  These will take
            precedence over tie_strs if different than ''.

        Returns
        -------
        repr_vpfit : str

        """
        # get Doppler parameter to km/s
        b = b.to('km/s').value

        # Ion name
        name = ions.ion_name(self.Zion, nspace=1)
        name = name.replace(' ', '')

        # Deal with fix and tie parameters
        # Check format first
        for i, x_strs in enumerate([tie_strs, fix_strs]):
            if (not isinstance(x_strs, tuple)) or (not all(isinstance(s, (str, basestring)) for s in x_strs)):
                if i == 0:
                    raise TypeError('`tie_strs` must be a tuple of strings.')
                elif i == 1:
                    raise TypeError('`fix_strs` must be a tuple of strings.')
            if len(x_strs) != 3:
                raise SyntaxError('`tie_strs` and `fix_strs` must have len() == 3')

        # reformat for VPFIT standard
        fix_strs = np.array([s.upper() for s in fix_strs])
        tie_strs = np.array([s.lower() for s in tie_strs])
        # preference to fix_strs over tie_strs
        strs = np.where(fix_strs != '', fix_strs, tie_strs)

        # create the line string
        s = '{:s} {:.5f}{:s} {:.5f} {:.2f}{:s} {:.2f} {:.2f}{:s} {:.2f}'.format(name, self.zcomp, strs[0], 0, b,
                                                                                strs[1], 0, self.logN, strs[2], 0)
        if len(self.comment) > 0:
            s += '! {:s}'.format(self.comment)
        s += '\n'
        return s

    def repr_alis(self, T_kin=1e4*u.K, bturb=0.*u.km/u.s,
                  tie_strs=('', '', '', ''), fix_strs=('', '', '', '')):
        """
        String representation for ALIS (line fitting software)

        Parameters
        ----------
        T_kin : Quantity, optional
            Kinetic temperature. Default 1e4*u.K
        bturb : Quantity, optional
            Turbulent Doppler parameter. Default 0.*u.km/u.s
        tie_strs : tuple of strings, optional
            Strings to be used for tying parameters
            (logN,z,bturb,T_kin), respectively.  These are all
            converted to lower case format, following ALIS convention.
        fix_strs : tuple of strings, optional
            Strings to be used for fixing parameters
            (logN,z,bturb,T_kin), respectively.  These are all
            converted to upper case format, following ALIS convention.
            These will take precedence over tie_strs if different from
            ''.

        Returns
        -------
        repr_alis : str

        """

        # Convert to the units ALIS wants
        T_kin = T_kin.to('K').value
        bturb = bturb.to('km/s').value

        # A patch for nucleons; todo: come up with a better way to do this using ELEMENTS?
        if self.Zion[0] == 1:
            nucleons = 1
        elif self.Zion[0] > 1:
            nucleons = 2 * self.Zion[0]

        # name
        name = ions.ion_name(self.Zion, nspace=1)
        name = '{}'.format(nucleons)+name.replace(' ', '_')

        # Deal with fix and tie parameters
        # Check format first
        for i, x_strs in enumerate([tie_strs, fix_strs]):
            if (not isinstance(x_strs, tuple)) or (not all(isinstance(s, (str, basestring)) for s in x_strs)):
                if i == 0:
                    raise TypeError('`tie_strs` must be a tuple of strings.')
                elif i == 1:
                    raise TypeError('`fix_strs` must be a tuple of strings.')
            if len(x_strs) != 4:
                raise SyntaxError('`tie_strs` and `fix_strs` must have len()== 4')

        # reformat for ALIS standard
        fix_strs = np.array([s.upper() for s in fix_strs])
        tie_strs = np.array([s.lower() for s in tie_strs])
        # preference to fix_strs over tie_strs
        strs = np.where(fix_strs != '', fix_strs, tie_strs)

        s = 'voigt   ion={:s} {:.2f}{:s} redshift={:.5f}{:s} {:.1f}{:s} {:.1E}{:s}'.format(name, self.logN, strs[0],
                                                                                           self.zcomp, strs[1], bturb,
                                                                                           strs[2], T_kin, strs[3])

        if len(self.comment) > 0:
            s += '# {:s}'.format(self.comment)
        s += '\n'
        return s

    def repr_joebvp(self, specfile, flags=(2,2,2), b_default=10*u.km/u.s):
        """
        String representation for JOEBVP (line fitting software).

        Parameters
        ----------
        specfile : str
            Name of the spectrum file
        flags : tuple of ints, optional
            Flags (nflag, bflag, vflag). See JOEBVP input for details
            about these flags.
        b_default : Quantity, optional
            Doppler parameter value adopted in case an absorption
            line within the component has not set this attribute
            Default is 10 km/s.

        Returns
        -------
        repr_joebvp : str
            May contain multiple "\n" (1 per absline within component)

        """
        # Reference:
        # specfile|restwave|zsys|col|bval|vel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|trans
        s = ''
        for aline in self._abslines:
            s += '{:s}|{:.5f}|'.format(specfile, aline.wrest.to('AA').value)
            logN = aline.attrib['logN']
            b_val = aline.attrib['b'].to('km/s').value
            if b_val == 0:  # set the default
                b_val = b_default.to('km/s').value
            s += '{:.8f}|{:.4f}|{:.4f}|0.|'.format(self.zcomp, logN, b_val)  # `vel` is set to 0. because z is zcomp
            s += '{}|{}|{}|'.format(int(flags[0]), int(flags[1]), int(flags[2]))
            vlim = aline.limits.vlim.to('km/s').value
            wvlim = aline.limits.wvlim.to('AA').value
            s += '{:.4f}|{:.4f}|{:.5f}|{:.5f}|'.format(vlim[0], vlim[1], wvlim[0], wvlim[1])
            s += '{:s}'.format(aline.data['ion_name'])

            if len(self.comment) > 0:
                s += '# {:s}'.format(self.comment)
            s += '\n'
        return s

    def stack_plot(self, return_fig=False, **kwargs):
        """Show a stack plot of the component, if spec are loaded
        Assumes the data are normalized.

        Parameters
        ----------
        return_fig : bool, optional
            If True, return stack plot as plt.Figure() instance for further manipulation

        Returns
        -------
        fig : matplotlib Figure, optional
            Figure instance containing stack plot with subplots, axes, etc.
        """
        if return_fig:
            fig = ltap.stack_plot(self._abslines, vlim=self.vlim, return_fig=True, **kwargs)
            return fig
        else:
            ltap.stack_plot(self._abslines, vlim=self.vlim, **kwargs)

    def to_dict(self):
        """ Convert component data to a dict
        Returns
        -------
        cdict : dict
        """
        cdict = dict(Zion=self.Zion, zcomp=self.zcomp, vlim=self.vlim.to('km/s').value,
                     Name=self.name,
                     RA=self.coord.ra.value, DEC=self.coord.dec.value,
                     A=self.A, Ej=self.Ej.to('1/cm').value, comment=self.comment,
                     flag_N=self.flag_N, logN=self.logN, sig_logN=self.sig_logN)
        cdict['class'] = self.__class__.__name__
        # AbsLines
        cdict['lines'] = {}
        for iline in self._abslines:
            cdict['lines'][iline.wrest.value] = iline.to_dict()
        # Polish
        cdict = ltu.jsonify(cdict)
        # Return
        return cdict

    def copy(self):
        """ Generate a copy of itself
        Returns
        -------
        abscomp : AbsComponent

        """
        # Instantiate with required attributes
        abscomp = AbsComponent(self.coord, self.Zion, self.zcomp, self.vlim)
        # Add in the rest
        attrs = vars(self).keys()
        for attr in attrs:
            if attr == '_abslines':
                for iline in self._abslines:
                    abscomp._abslines.append(iline.copy())
            else:
                setattr(abscomp, attr, getattr(self, attr))
        # Return
        return abscomp

    def __getitem__(self, attrib):
        """Passback attribute, if it exists

        Useful for columns

        Parameters
        ----------
        attrib : str
        """
        return getattr(self, attrib)

    def __repr__(self):
        txt = '<{:s}: {:s} {:s}, Name={:s}, Zion=({:d},{:d}), Ej={:g}, z={:g}, vlim={:g},{:g}'.format(
            self.__class__.__name__, self.coord.ra.to_string(unit=u.hour,sep=':', pad=True),
                self.coord.dec.to_string(sep=':',pad=True,alwayssign=True), self.name, self.Zion[0], self.Zion[1], self.Ej, self.zcomp, self.vlim[0], self.vlim[1])

        # Column?
        if self.flag_N > 0:
            txt = txt + ', logN={:g}'.format(self.logN)
            txt = txt + ', sig_logN={:g}'.format(self.sig_logN)
            txt = txt + ', flag_N={:d}'.format(self.flag_N)

        # Finish
        txt = txt + '>'
        return (txt)

