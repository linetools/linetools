"""
Module for dealing with LSFs of various astronomical instruments.
"""
from __future__ import print_function

import numpy as np
from astropy.io import fits, ascii
from astropy.units import Quantity
from astropy.table import Table, Column
import astropy.units as u
import glob, imp

lt_path = imp.find_module('linetools')[1]

class LSF(object):

    """Class to deal with line-spread-functions (LSFs) from
    various different astornomical spectrographs.

    Note: only implemented for HST/COS at the moment.

    Parameters
    ----------
    wv_array : Quantity numpy.ndarray, shape(N,)
        Wavelength array for which the LSF kernel is defined. The 
        central wavelength value of `wv_array` define the wavelength
        at which the LSF is defined, while the limits of `wv_array` 
        define the extent of the kernel.
    instr_config : dict
        A dictionary with the instrument configuration details relevant
        to the required LSF. Mandatory keywords of the dict are: ['name'], 
        all of which must be either string or None. 
        Note: There must be extra relevant keywords specific to each instrument.



    Attributes
    ----------

    """

    def __init__(self, wv_array, instr_config):
        #mandatory keys for characterizing a spectrograph mode
        self.mandatory_dict_keys = ['name']
                
        #Check correct format
        if not ((isinstance(wv_array, np.ndarray)) or (isinstance(wv_array, Quantity))):
            raise SyntaxError('`wv_array` must be Quantity numpy.ndarray')
        elif len(wv_array.shape) != 1:
            raise SyntaxError('`wv_array` must be of shape(N,), i.e. 1-dimensional array')
        elif not isinstance(instr_config,dict):
            raise SyntaxError('`instr_config` must be a dictionary.')
        elif not all([key in instr_config.keys() for key in self.mandatory_dict_keys]):
            raise SyntaxError('`instr_config` must have the following mandatory keys {}:'.format(self.mandatory_dict_keys))

        #Initialize
        self.wv_array = wv_array
        self.wv_min = np.min(wv_array)
        self.wv_max = np.max(wv_array)
        self.wv0 = 0.5 * (self.wv_max - self.wv_min)
        self.instr_config = instr_config
        self.name = instr_config['name']

        #load data
        self.rel_pix , self._data = self.load_data()

    def load_data(self):
        """Load the right data according to `instr_config`"""
        if self.name not in ['COS']:
            raise ValueError('Not ready for this instrument: {}'.format(self.name))

        if self.name == 'COS':
            #define pixel scales; values obtained from STScI
            #these values must be consistent with the given LSFs
            pixel_scale = {'G130M': 9.97 / 1000. * u.AA,
                        'G160M': 12.23 / 1000. * u.AA,
                        'G140L': 80.3 / 1000. * u.AA,
                        'G230L': 390. / 1000. * u.AA,
                        'G185M': 37. / 1000. * u.AA,
                        'G225M': 33. / 1000. * u.AA,
                        'G285M': 40. / 1000. * u.AA}
            #define chanel based on grating name
            channel = {'G130M':  'FUV',
                        'G160M': 'FUV',
                        'G140L': 'FUV',
                        'G230L': 'NUV',
                        'G185M': 'NUV',
                        'G225M': 'NUV',
                        'G285M': 'NUV'}

            try:
                grating = self.instr_config['grating']
            except:
                raise SyntaxError('`grating` keyword missing in `instr_config` dictionary.')
            
            if grating not in channel.keys():
                raise ValueError('Not ready for this HST/COS grating: {}'.format(grating))

            if channel[grating] == 'NUV': #there is only 1 LSF file for NUV data
                file_name = 'nuv_all_lp1.txt'
            #COS
            elif channel[grating] == 'FUV':
                # Use the ones corrected by scattering when possible
                # (currently, these are only available for lifetime-position 1)
                # check: http://www.stsci.edu/hst/cos/performance/spectral_resolution
                try:
                    life_position = self.instr_config['life_position']
                except:
                    raise SyntaxError('`life_position` keyword missing in `instr_config` dictionary.')

                if life_position not in ['1','2']:
                    raise ValueError('HST/COS `life_position` should be either `1` or `2` (strings)') 

                if life_position == '1':
                    if grating == 'G140L': #use theoretical values 
                        file_name = 'fuv_G140L_lp1.txt'
                        
                    elif grating == 'G130M': #use empirical values corrected by scattering
                        file_name = 'fuv_G130M_lp1_empir.txt'

                    elif grating == 'G160M': #use empirical values corrected by scattering
                        file_name = 'fuv_G160M_lp1_empir.txt'
                
                elif life_position == '2':
                    try:
                        cen_wave = self.instr_config['cen_wave']
                    except:
                        raise SyntaxError('`cen_wave` keyword missing in `instr_config` dictionary. This should point to the central wavelength of the grating in Angstroms.')
                    #adjust format in cases where cen_wave is of the form: str(1230A)
                    if cen_wave.endswith('A'): #adjust format
                        cen_wave = cen_wave[:-1]
                    
                    #filenames in this case have a well defined naming convention
                    file_name = 'fuv_{}_{}_lp2.txt'.format(grating,cen_wave)
            
            else: #this should never happen
                raise ValueError('Not ready for the given HST/COS channel; only `NUV` and `FUV` channels allowed.')
            
            #still in COS
            #point to the right file
            file_name = lt_path + '/data/lsf/{}/{}'.format(self.name,file_name)
            
            #get column names
            f = open(file_name,'r')
            line = f.readline() #first line of file
            f.close()
            #get rid of '\n'
            line = line.split('\n')[0]
            #by construction first column should be separated by `,`
            col_names = line.split(',')
            col_names[0] = 'rel_pix'
            
            rel_pix = pixel_scale[grating] #read from dictionary defined above
            #read data
            data = ascii.read(file_name,data_start=1,names=col_names)
            
            return rel_pix , data