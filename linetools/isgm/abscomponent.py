""" Classes for absorption line component
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
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Column

from specutils import Spectrum1D

from linetools.analysis import absline as ltaa
from linetools.spectralline import AbsLine
from linetools.abund import ions

#import xastropy.atomic as xatom
#from xastropy.stats import basic as xsb
#from xastropy.xutils import xdebug as xdb

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
    def from_abslines(cls, abslines):
        """Instantiate from a list of AbsLine objects

        Parameters
        ----------
        abslines : list 
          List of AbsLine objects
        """
        # Check
        if not isinstance(abslines,list):
            raise IOError("Need a list of AbsLine objects")
        if not all(isinstance(x,AbsLine) for x in abslines):
            raise IOError("List needs to contain only AbsLine objects")

        # Instantiate with the first line
        init_line = abslines[0]
        slf = cls( init_line.attrib['coord'], (init_line.data['Z'],init_line.data['ion']), init_line.attrib['z'], init_line.analy['vlim'], Ej=init_line.data['Ej'])
        slf._abslines.append(init_line)
        # Append with component checking
        if len(abslines) > 1:
            for absline in abslines[1:]:
                slf.add_absline(absline)
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
        if not isinstance(component,AbsComponent):
            raise IOError('Need an AbsComponent object')
        # Return
        return cls(component.coord, component.Zion, component.zcomp, component.vlim, Ej=component.Ej, A=component.A, **kwargs)

    def __init__(self, radec, Zion, z, vlim, Ej=0./u.cm, A=None, Ntup=None,comment=''):
        """  Initiator

        Parameters
        ----------
        radec : tuple or SkyCoord
            (RA,DEC) in deg or astropy.coordinate
        Zion : tuple 
            Atomic number, ion -- (int,int)
            e.g. (8,1) for OI
        z : float
            Absorption component redshift
        vlim : Quantity array
            Velocity limits of the component w/r to `z`
            e.g.  [-300,300]*u.km/u.s
        A : int, optional
            Atomic mass -- used to distinguish isotopes
        Ntup : tuple
            (int,float,float)
            (flag_N,logN,sig_N)
            flag_N : Flag describing N measurement
            logN : log10 N column density
            sig_logN : Error in log10 N
        Ej : Quantity, optional
            Energy of lower level (1/cm)
        comment : str, optional
            A comment, default is ``
        """

        # Required
        if isinstance(radec,(tuple)):
            self.coord = SkyCoord(ra=radec[0], dec=radec[1])
        elif isinstance(radec,SkyCoord):
            self.coord = radec
        self.Zion = Zion
        self.zcomp = z
        self.vlim = vlim

        # Optional
        self.A = A
        self.Ej = Ej
        self.comment = comment
        if Ntup is not None:
            self.flag_N = Ntup[0]
            self.logN = Ntup[1]
            self.sig_logN = Ntup[2]
            _,_ = ltaa.linear_clm(self) # Set linear quantities
        else:
            self.flag_N = 0
            self.logN = 0.
            self.sig_logN = 0.

        # Other
        self.name = ions.ion_name(self.Zion, nspace=1)
        self._abslines = []

    def add_absline(self,absline,tol=0.1*u.arcsec):
        """Add an AbsLine object to the component if it satisfies
        all of the rules.

        For velocities, we demand that the new line has a velocity
        range that is fully encompassed by the component.

        Parameters
        ----------
        absline : AbsLine
        tol : Angle, optional
          Tolerance on matching coordinates
        """
        # Perform easy checks
        test = bool(self.coord.separation(absline.attrib['coord']) < tol)
        test = test & (self.Zion[0] == absline.data['Z'])
        test = test & (self.Zion[1] == absline.data['ion'])
        test = test & bool(self.Ej == absline.data['Ej'])
        # Now redshift/velocity
        zlim_line = (1+absline.attrib['z'])*absline.analy['vlim']/const.c.to('km/s')
        zlim_comp = (1+self.zcomp)*self.vlim/const.c.to('km/s')
        test = test & (zlim_line[0]>=zlim_comp[0]) & (zlim_line[1]<=zlim_comp[1])
        # Isotope
        if self.A is not None:
            raise ValueError('Not ready for this yet')
        # Append?
        if test:
            self._abslines.append(absline)
        else:
            warnings.warn('Input absline with wrest={:g} does not match component rules. Not appending'.format(absline.wrest))

    def build_table(self):
        """Generate an astropy QTable out of the component.
        Returns
        -------
        comp_tbl : QTable
        """
        if len(self._abslines) == 0:
            return
        comp_tbl = QTable()
        comp_tbl.add_column(Column([iline.wrest.to(u.AA).value for iline in self._abslines]*u.AA,name='wrest'))
        for attrib in ['z', 'flag_N', 'logN', 'sig_logN']:
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
        sig_EW=np.array([aline.attrib['sig_EW'].to('AA').value for aline in self._abslines])*u.AA
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
        try: # Nicer view, especially in notebook
            import seaborn as sns; sns.set(context="notebook",font_scale=2)
        except ImportError:
            pass
        mpl.rcParams['font.family'] = 'stixgeneral'
        mpl.rcParams['font.size'] = 15.
        # Check for spec
        gdiline = []
        for iline in self._abslines:
            if isinstance(iline.analy['spec'],Spectrum1D):
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
        for qq,iline in enumerate(gdiline):
            # Calculate
            velo = iline.analy['spec'].relative_vel((1+iline.attrib['z'])*iline.wrest)
            cst = atom_cst/(iline.data['f']*iline.wrest) #/ (u.km/u.s) / u.cm * (u.AA/u.cm)
            Na = np.log(1./np.maximum(iline.analy['spec'].flux,
                                      iline.analy['spec'].sig))*cst
            # Figure out ymnx
            pixmnx = (velo > self.vlim[0]) & (velo < self.vlim[1])
            if iline.data['f']*iline.wrest > fw_sv:
                ymax = max(np.max(Na[pixmnx].value),ymax)
                fw_sv = iline.data['f']*iline.wrest
            # Plot
            ax.plot(velo, Na, '-', linestyle='steps-mid', label=iline.data['name'])
            #ax.plot(velo, iline.analy['spec'].sig, 'r:')
        # Axes
        ax.set_xlim(self.vlim.value)
        ax.set_ylim(-0.2*ymax, 5*ymax)
        #ax.set_ylim(ymnx)
        ax.minorticks_on()
        ax.set_xlabel('Relative Velocity (km/s)')
        ax.set_ylabel(r'Apparent Column (cm$^{-2}$ per km/s)')
        # Legend
        legend = ax.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                           handletextpad=0.3, fontsize='large')

        plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
        if show:
            plt.show()
        plt.close()

    def synthesize_colm(self, overwrite=False, redo_aodm=False, **kwargs):
        """Synthesize column density measurements of the component.
        Default is to use the current AbsLine values, but the user can
        request that those be re-calculated with AODM

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
            raise IOError("Column densities already set.  Use clobber=True to redo.")
        # Redo?
        if redo_aodm:
            for aline in self._abslines:
                aline.measure_aodm(**kwargs)
        # Collate
        self.flag_N = 0
        for aline in self._abslines:
            if np.allclose(aline.attrib['N'].value,0.):
                raise ValueError("Need to set N in attrib.  \n Consider linear_clm in linetools.analysis.absline")
            if aline.attrib['flag_N'] == 1: # Good value?
                if self.flag_N == 1: # Weighted mean
                    # Original
                    weight = 1. / self.sig_N**2
                    mu= self.N * weight
                    # Update
                    weight += 1./aline.attrib['sig_N']**2
                    self.N = (mu + aline.attrib['N']/aline.attrib['sig_N']**2) / weight
                    self.sig_N = np.sqrt(1./weight)
                else: # Fill
                    self.N = aline.attrib['N']
                    self.sig_N = aline.attrib['sig_N']
                    self.flag_N = 1
            elif aline.attrib['flag_N'] == 2: # Lower limit
                if self.flag_N in [0,3]:
                    self.N = aline.attrib['N']
                    self.sig_N = 99.
                    self.flag_N = 2
                elif self.flag_N == 2:
                    self.N = max(self.N,aline.attrib['N'])
                    self.sig_N = 99.
                elif self.flag_N == 1:
                    pass
            elif aline.attrib['flag_N'] == 3: # Upper limit
                if self.flag_N == 0:
                    self.N = aline.attrib['N']
                    self.sig_N = aline.attrib['sig_N']
                    self.flag_N = 3
                elif self.flag_N in [1,2]:
                    pass
                elif self.flag_N == 3:
                    if aline.attrib['N'] < self.N:
                        self.N = aline.attrib['N']
                        self.sig_N = aline.attrib['sig_N']
            else:
                raise ValueError("Bad flag_N value")
        # Log values
        self.logN, self.sig_logN = ltaa.log_clm(self)

    def stack_plot(self, nrow=6, show=True):
        """Show a stack plot of the component, if spec are loaded
        Assumes the data are normalized.

        Parameters
        ----------
        nrow : int, optional  
          Maximum number of rows per column
        show : bool, optional
          Show the plot?
        """
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
        import matplotlib as mpl
        mpl.rcParams['font.family'] = 'stixgeneral'
        mpl.rcParams['font.size'] = 15.
        # Check for spec
        gdiline = []
        for iline in self._abslines:
            if isinstance(iline.analy['spec'],Spectrum1D):
                gdiline.append(iline)
        nplt = len(gdiline)
        if nplt == 0:
            print("Load spectra into the absline.analy['spec']")
            return
        # Setup plot
        nrow = min(nplt,nrow)
        ncol = nplt // nrow + (nplt % nrow > 0)
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)
        ymnx = (-0.1,1.1)

        for qq,iline in enumerate(gdiline):
            ax = plt.subplot(gs[qq%nrow, qq//nrow])
            # Plot
            velo = iline.analy['spec'].relative_vel((1+iline.attrib['z'])*iline.wrest)
            ax.plot(velo, iline.analy['spec'].flux, 'k-', linestyle='steps-mid')
            ax.plot(velo, iline.analy['spec'].sig, 'r:')
            # Lines
            ax.plot([0]*2, ymnx, 'g--')
            # Axes
            ax.set_xlim(self.vlim.value)
            ax.set_ylim(ymnx)
            ax.minorticks_on()
            if ((qq+1) % nrow == 0) or ((qq+1) == nplt):
                ax.set_xlabel('Relative Velocity (km/s)')
            else:
                ax.get_xaxis().set_ticks([])
            # Label
            ax.text(0.1, 0.1, iline.data['name'], transform=ax.transAxes, ha='left', va='center', fontsize='x-large')#, bbox={'facecolor':'white'})

        plt.tight_layout(pad=0.2,h_pad=0.,w_pad=0.1)
        if show:
            plt.show()
        plt.close()


    def repr_vpfit(self, b=10.):
        """
        String representation for VPFIT in its fort.26 format

        Parameters
        ----------
        b : float, optional
            Doppler parameter of the component

        Returns
        -------
        repr_vpfit : str
        """
        name = self.name.replace(' ','')
        try:
            logN = self.logN
        except:
            logN = 0.

        s = '{} {:.5f} {:.5f} {:.2f} {:.2f} {:.2f} {:.2f}'.format(name,self.zcomp,0,b,0,logN,0)
        if len(self.comment)>0:
            s += '! {}'.format(self.comment)
        return s


    def __getitem__(self, attrib):
        """Passback attribute, if it exists

        Useful for columns

        Parameters
        ----------
        attrib : str
        """
        return getattr(self,attrib)

    def __repr__(self):
        txt = '<{:s}: {:s} {:s}, Name={:s}, Zion=({:d},{:d}), Ej={:g}, z={:g}, vlim={:g},{:g}'.format(
            self.__class__.__name__, self.coord.ra.to_string(unit=u.hour,sep=':', pad=True), self.coord.dec.to_string(sep=':',pad=True,alwayssign=True), self.name, self.Zion[0], self.Zion[1], self.Ej, self.zcomp, self.vlim[0], self.vlim[1])

        # Column?
        if self.flag_N > 0:
            txt = txt + ', logN={:g}'.format(self.logN)
            txt = txt + ', sig_N={:g}'.format(self.sig_logN)
            txt = txt + ', flag_N={:d}'.format(self.flag_N)

        # Finish
        txt = txt + '>'
        return (txt)

