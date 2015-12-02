""" Subclasses for DLA AbsSystem and AbsSurvey
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import pdb

from astropy import units as u

from linetools.isgm.abssystem import AbsSystem
from linetools.isgm.abssurvey import AbslineSurvey
from linetools.isgm import utils as ltiu


class DLASystem(AbsSystem):
    """
    Class for a DLA absorption system
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
        slf._datdict = datdict

        # QSO keys
        slf.qso = slf._datdict['QSO name']
        slf.zqso = float(slf._datdict['QSO zem'])
        # Name
        slf.name = '{:s}_z{:0.3f}'.format(slf.qso,zabs)

        # Abund
        slf.flg_ZH = float(slf._datdict['flg_mtl'])
        slf.ZH = float(slf._datdict['[M/H]'])
        slf.sig_ZH = float(slf._datdict['sig([M/H])'])

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
          Defaulted to +/- 500 km/s if None
        NHI : float, required despite being a keyword
          log10 of HI column density
          must be 20.3 or higher
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
        if NHI < 20.3:
            raise ValueError("This is not a DLA!  Try an LLS (or SLLS)")
        # vlim
        if vlim is None:
            vlim = [-500., 500.]*u.km/u.s
        # Generate with type
        AbsSystem.__init__(self, 'DLA', radec, zabs, vlim, NHI=NHI, **kwargs)

        # Other
        self.ZH = 0.

    def get_ions(self, use_clmfile=False, update_zvlim=True, linelist=None):
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
        reload(ltiu)
        if use_clmfile:
            clm_fil = self.tree+self._datdict['Abund file']
            # Read
            self._clmdict = ltiu.read_clmfile(clm_fil, linelist=linelist)
            # Build components
            components = ltiu.build_components_from_abslines([], clmdict=self._clmdict, coord=self.coord)
            # Read .ion file and fill in components
            ion_fil = self.tree+self._clmdict['ion_fil']
            self._indiv_ionclms = ltiu.read_ion_file(ion_fil, components)
            # Parse .all file
            all_file = ion_fil.split('.ion')[0]+'.all'
            self.all_file=all_file #MF: useful to have
            _ = ltiu.read_all_file(all_file, components=components)
            # Build table
            self._ionN = ltiu.iontable_from_components(components, ztbl=self.zabs)
            # Add to AbsSystem
            for comp in components:
                self.add_component(comp)
        else:
            raise IOError("Not ready for this")

    # Output
    def __repr__(self):
        return ('[{:s}: {:s} {:s}, {:g}, NHI={:g}, Z/H={:g}]'.format(
                self.__class__.__name__,
                 self.coord.ra.to_string(unit=u.hour, sep=':', pad=True),
                 self.coord.dec.to_string(sep=':', pad=True),
                 self.zabs, self.NHI, self.ZH))

    def print_abs_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'DLA'

# #######################################################################
# #######################################################################
# #######################################################################
# Class for DLA Survey
class DLASurvey(AbslineSurvey):
    """An DLA Survey class

    Attributes:
        
    """
    @classmethod
    def default_sample(cls):
        """
        Returns
        -------
        dlasurvey : AbslineSurvey
        """
        # Default sample of DLA:  Neeleman
        if os.getenv('DLA') is None:
            print('Need to grab the DLA tree from JXP')
            return None
        dlasurvey = cls.from_flist('Lists/Neeleman13.lst', tree=os.environ.get('DLA'))
        dlasurvey.ref = 'Neeleman+13'

        # Return
        return dlasurvey

    def __init__(self, **kwargs):
        # Generate with type
        AbslineSurvey.__init__(self, 'DLA', **kwargs)














## #################################
## #################################
## TESTING
## #################################
if __name__ == '__main__':

    flg_test = 0
    #flg_test = 1  # ions
    #
    #flg_test += 2**9 # DLA Survey NHI
    flg_test += 2**10 # DLA Survey ions

    # Test Absorption System
    print('-------------------------')
    tmp1 = DLA_System(dat_file='Data/PH957.z2309.dat',
                      tree=os.environ.get('DLA'))
    print(tmp1)

    # Test ions
    if (flg_test % 2**1) >= 2**0:
        print('-------------------------')
        clm_fil = tmp1.tree+tmp1.datdict['Abund file']
        tmp1.get_ions()
        # Print
        print('Si II: ')
        print(tmp1.ions[(14,2)])
        print(tmp1.ions.trans[15]) # CIV 1550



    # #############################
    # DLA Survey
    if (flg_test % 2**10) >= 2**9:
        print('-------------------------')
        dla = DLA_Survey('Lists/metal_MAR_all.lst', tree=os.environ.get('DLA'))
        dla.fill_ions()
        xdb.xhist(dla.NHI, binsz=0.10)

    # DLA Survey ions
    if (flg_test % 2**11) >= 2**10:
        dla = DLA_Survey('Lists/metal_MAR_all.lst', tree=os.environ.get('DLA'))
        dla.fill_ions()
        xdb.xhist(dla.ions((6,4),skip_null=True)['clm'], binsz=0.3,
                  xlabel=r'$\log_{10} N({\rm C}^{+3})$')
        xdb.xhist(dla.ions((14,2),skip_null=True)['clm'], binsz=0.3,
                  xlabel=r'$\log_{10} N({\rm Si}^{+})$')
