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

#import xastropy.atomic as xatom
#from xastropy.stats import basic as xsb
#from xastropy.xutils import xdebug as xdb

# Class for Components
class AbsComponent(object):
    """
    Class for an absorption component

    Attributes
    ----------
    coord : SkyCoord
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
            raise IOError('Need a list of AbsLine objects')
        if not all(isinstance(x,AbsLine) for x in abslines):
            raise IOError('List needs to contain AbsLine objects')
        # Instantiate with the first line
        init_line = abslines[0]
        slf = cls( init_line.attrib['coord'],
           (init_line.data['Z'],init_line.data['ion']),
           init_line.attrib['z'], init_line.analy['vlim']) 
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

    def __init__(self, radec, Zion, z, vlim, Ej=0./u.cm, A=None, Ntup=None):
        """  Initiator

        Parameters
        ----------
        radec : tuple or SkyCoord
          (RA,DEC) in deg or astropy.coordinate
        Zion : tuple 
          Atomic number, ion -- (int,int) 
             e.g. (8,1) for OI
        z : float
          Absorption redshift
        vlim : Quantity array
          Velocity limits of the component
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
        if Ntup is not None:
            self.flag_N = Ntup[0]
            self.logN = Ntup[1]
            self.sig_logN = Ntup[2]
            _,_ = ltaa.linear_clm(self) # Set linear quantities
        else:
            self.flag_N = 0

        # Other
        self._abslines = []

    def add_absline(self,absline,toler=0.1*u.arcsec):
        """Add an AbsLine object to the component if it satisfies
        all of the rules.

        For velocities, we demand that the new line has a velocity
        range that is fully encompassed by the component.

        Parameters
        ----------
        absline : AbsLine
        toler : Angle, optional
          Tolerance on matching coordinates
        """
        # Perform easy checks
        test = bool(self.coord.separation(absline.attrib['coord']) < toler)
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

    def __getitem__(self, attrib):
        """Passback attribute, if it exists

        Useful for columns

        Parameters
        ----------
        attrib : str
        """
        return getattr(self,attrib)

    def __repr__(self):
        txt = '[{:s}: {:s} {:s}, Zion=({:d},{:d}), z={:g}, vlim={:g},{:g}'.format(
            self.__class__.__name__,
            self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
            self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
            self.Zion[0], self.Zion[1], self.zcomp,
            self.vlim[0],self.vlim[1])
        # Column?
        try:
            txt = txt+', logN={:g}'.format(self.logN)
        except AttributeError:
            pass
        else:
            txt = txt+', sig_N={:g}'.format(self.sig_logN)
        # Finish
        txt = txt + ']'
        return (txt)

def chk_components(components, chk_match=False, chk_A_none=False, toler=0.2*u.arcsec):
    """ Performs checks on a list of components

    Parameters
    ----------
    components : list
      list of AbsComponent objects
    chk_match : bool, optional
      if True, require that the components match in RA/DEC, Zion, Ej, A but not velocity
    chk_A_none : bool, optional
      if True, require that A *not* be set
    toler : Angle, optional
      Tolerance on matching coordinates
    """
    tests = True
    # List
    if not isinstance(components,list):
        tests = False
        raise IOError('Need a list of AbsComponent objects')
    # Object
    if not all(isinstance(x,AbsComponent) for x in components):
        tests = False
        raise IOError('List needs to contain only AbsComponent objects')
    # A None
    if chk_A_none:
        if not all(x.A is None for x in components):
            tests = False
            raise IOError('Not ready for components with A set')
    # Matching?
    if chk_match:
        match = True
        comp0 = components[0]
        for comp in components[1:]:
            # RA/DEC
            match = match & bool(comp0.coord.separation(comp.coord) < toler)
            # Zion
            match = match & (comp0.Zion == comp.Zion)
            # Ej
            match = match & np.testing.assert_almost_equal(comp0.Ej.to('1/cm').value,comp.Ej.to('1/cm').value)
            # A
            match = match & (comp0.A == comp.A)
        tests = tests & match
    # Return
    return tests

def iontable_from_components(components):
    """Generate a QTable from a list of components

    Method does *not* perform logic on redshifts or vlim.
    Includes rules for adding components of like ion
    Not ready for varying atomic mass (e.g. Deuterium)

    Parameters
    ----------
    components : list
      list of AbsComponent objects
    """
    # Checks
    assert chk_components(components,chk_A_none=True)

    # Set z from mean
    ztbl = np.mean(comp.zcomp for comp in components)

    # Construct the QTable
    iontbl = QTable()

    # Identify unique Zion, Ej (not ready for A)
    uZiE = np.array([comp.Zion[0]*1000000+comp.Zion[1]*10000+
                      comp.Ej.to('1/u.cm').value for comp in components])
    uniZi, auidx = np.unique(uZiE, return_index=True)

    # Loop
    for uidx in auidx:
        # Synthesize components with like Zion, Ej
        mtZiE = uZiE == uZiE[idx]
        synthcomp = synthesize_comps(components[idx[mtZiE]],zcomp=ztbl)
        # Add a row to QTable


    # Add zlim to metadata

    # Return

def synthesize_components(components, zcomp=None):
    """Synthesize a list of components into one

    Requires consistent RA/DEC, Zion, Ej, (A; future)
    Is agnostic about z+vlim
    Melds column densities
    Melds velocities with a small buffer (10 km/s)

    Note: Could make this a way to instantiate AbsComponent

    Parameters
    ----------
    components : list
      list of AbsComponent objects
    zcomp : float, optional
      Input z to reference the synthesized component
      If not input, the mean of the input components is used
    """
    reload(ltaa)
    # Checks
    assert chk_components(components,chk_A_none=True,chk_match=True)


    # Final component
    synth_comp = AbsComponent.from_component(components[0])
    synth_comp.logN = 0.
    synth_comp.sig_logN = 0.

    # Meld column densities
    for comp in components[1:]:
        synth_comp.flag_N, synth_comp.logN, synth_comp.sig_logN = ltaa.sum_logN(synth_comp,comp)

    # Meld z, vlim
    # zcomp
    if zcomp is None:
        zcomp = np.mean(comp.zcomp for comp in components)
    synth_comp.zcomp = zcomp
    # Set vlim by min/max  [Using non-relativistic + buffer]
    v = [(comp.zcomp-zcomp)/(1+zcomp)*3e5 for comp in components]
    synth_comp.vlim = [np.min(v)-10.,np.max(v)+10.]*u.km/u.s
