"""
#;+ 
#; NAME:
#; spectralline
#;    Version 1.0
#;
#; PURPOSE:
#;    Module for SpectralLine class
#;   23-Jun-2015 by JXP
#;-
#;------------------------------------------------------------------------------
"""
from __future__ import print_function, absolute_import, division, unicode_literals

# Python 2 & 3 compatibility
try:
    basestring
except NameError:
    basestring = str

import pdb

from astropy import constants as const
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import QTable, Column

from specutils import Spectrum1D

#import xastropy.atomic as xatom
#from xastropy.stats import basic as xsb
#from xastropy.xutils import xdebug as xdb

# Class for Components
class AbsComponent(object):
    """Class for a spectral line.  Emission or absorption 
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

    # Initialize with wavelength
    def __init__(self, radec, Zion, z, vlim, Ej=0./u.cm, A=None, Ntup=None):
        """  Initiator

        Parameters
        ----------
        radec : tuple or SkyCoord
          (RA,DEC) in deg
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
          (flgN,logN,sigN)
          flgN : Flag describing N measurement
          logN : log10 N column density
          sigN : Error in log10 N
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
            self.flgN = Ntup[0]
            self.logN = Ntup[1]
            self.sigN = Ntup[2]

        # Other
        self._abslines = []

    def add_absline(self,absline):
        """Add an AbsLine object to the component if it satisfies
        all of the rules.
        For velocities, we demands that the new line has a velocity
          range that is fully encompassed by the component.

        Parameters
        ----------
        absline : AbsLine
        """
        # Perform easy checks
        test = self.coord.separation(absline.attrib['coord']) < 0.1*u.arcsec
        test = test & self.Zion[0] == absline.data['Z']
        test = test & self.Zion[1] == absline.data['ion']
        test = test & bool(self.Ej == absline.data['Ej'])
        # Now redshift/velocity
        zlim_line = (1+absline.attrib['z'])*absline.analy['vlim']/const.c.to('km/s')
        zlim_comp = (1+self.zcomp)*self.vlim/const.c.to('km/s')
        test = test & (zlim_line[0]>=zlim_comp[0]) & (zlim_line[1]>=zlim_comp[1]) 
        
        # Isotope
        if self.A is not None:
            raise ValueError('Not ready for this yet')
        # Append?
        if test:
            self._abslines.append(absline)
        else:
            print('Input absline with wrest={:g} does not match component rules'.format(absline.wrest))
            print('Not appending')

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
        for attrib in ['z', 'flagN', 'N', 'Nsig']:
            comp_tbl.add_column(Column([iline.attrib[attrib] for iline in self._abslines], name=attrib))
        # Return
        return comp_tbl

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


    # Output
    def __repr__(self):
        return ('[AbsComponent: {:s} {:s}, Zion=({:d},{:d}), z={:g}]'.format(
                self.coord.ra.to_string(unit=u.hour,sep=':',pad=True),
                self.coord.dec.to_string(sep=':',pad=True,alwayssign=True),
                self.Zion[0], self.Zion[1], self.zcomp))

