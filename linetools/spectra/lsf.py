"""
Module for dealing with LSFs of various astronomical instruments.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
from scipy.interpolate import interp1d
from astropy.io import fits, ascii
from astropy.units import Quantity
import astropy.units as u
from astropy.table import Table, QTable, Column
import glob, imp
from linetools.analysis.interp import interp_Akima
import warnings

lt_path = imp.find_module('linetools')[1]

class LSF(object):
    """Class to deal with line-spread-functions (LSFs) from
    various different astronomical spectrographs.

    Note: only implemented for HST/COS and HST/STIS at the moment.

    Parameters
    ----------
    instr_config : dict
        A dictionary with the instrument configuration details relevant
        to the required LSF. Mandatory keywords of the dict are: ['name'], 
        all of which must be either string or None. 
        Note: There must be extra relevant keywords specific to each instrument.

    """

    def __init__(self, instr_config):
        # mandatory keys for characterizing a spectrograph mode
        self.mandatory_dict_keys = ['name']
                
        # Check correct format
        if not isinstance(instr_config, dict):
            raise TypeError('`instr_config` must be a dictionary.')
        elif not all([key in instr_config.keys() for key in self.mandatory_dict_keys]):
            raise SyntaxError('`instr_config` must have the following mandatory keys {}:'.format(self.mandatory_dict_keys))

        # Initialize basics
        self.instr_config = instr_config
        self.name = instr_config['name']
        if self.name not in ['COS', 'STIS']:
            raise NotImplementedError('Not ready for this instrument: {}'.format(self.name))
        
        # initialize specific to given instrument name
        # only implemented for HST/COS  and HST/STIS so far
        if self.name == 'COS':
            self.pixel_scale, self._data = self.load_COS_data()
        elif self.name == 'STIS':
            self.pixel_scale, self._data = self.load_STIS_data()
        # IMPORTANT: make sure that LSFs are given in linear wavelength scales !!!

        #reformat self._data
        self.check_and_reformat_data()
        
        #other relevant values to initialize?

    def get_lsf(self, wv_array, kind='Akima'):
        """ Given a wavelength array `wv_array`, it returns
        the LSF kernel at the central wavelength of the array, 
        using the same pixel scale and extent of `wv_array`. 

        Method for non-Gaussian kernels: First, tabulated LSFs
        are linearly interpolated to the center of `wv_array`
        (see LSF.interpolate_to_wv0() for details); then, the
        LSF is interpolated to match the `wv_array`
        scale and extent using Akima (or cubic) interpolation (see
        LSF.interpolate_to_wv_array() for details).

        Parameters
        ----------
        wv_array : Quantity numpy.ndarray, shape(N,)
            Wavelength array for which the LSF kernel is defined. The 
            central wavelength value of `wv_array` define the wavelength
            at which the LSF is defined, while the limits of `wv_array` 
            define the extent of the kernel.
        kind : str, optional
            Specifies the kind of interpolation as a string either 
            ('cubic', 'Akima')

        Returns
        -------
        lsf_array : numpy.ndarray, shape(N,)
            The lsf kernel.

        """
        lsf_array = self.interpolate_to_wv_array(wv_array, kind=kind)
        return lsf_array['kernel'].data

    def check_and_reformat_data(self):
        """Any re-formating of self._data should happen here.

        At the moment this function does the following:
         - Make sure that the number of relative pixels of the LSF is odd integer
         - Impose the middle value to define the 0 relative pixel
         - Make sure tables with 'rel_pix' given in fraction of pixels (e.g. COS NUV) work properly
         - Normalize tabulated LSFs
        """

        self._data['rel_pix'] = self.check_and_reformat_relpix(self._data['rel_pix'].data)

        #normalize given LSFs
        for col_name in self._data.keys()[1:]:
            self._data[col_name] /= np.sum(self._data[col_name])

    def check_and_reformat_relpix(self, relpix):
        """Performs checks and reformating on a relative pixel array.

        At the moment this function does the following:
         - Make sure that the number of relative pixels of the LSF is odd integer
         - Impose the middle value to define the 0 relative pixel

        Parameters
        ----------
        relpix : np.array
            Relative pixel array

        Returns
        -------
        new_relpix array that conforms to the LSF class standards.
        """

        # odd integer for total number of relative pixels
        n_pix = len(relpix)
        assert n_pix % 2 != 0, ValueError('LSF tables must be given as odd integers!')

        # make sure relpix == 0 is in the middle of the array
        n_half = (n_pix - 1) / 2
        mid_value = relpix[n_half]
        new_relpix = relpix - mid_value  # 0 in the middle

        return new_relpix

    def load_COS_data(self):
        """Load the right data according to `instr_config` for HST/COS 
        instrument"""

        # define pixel scales; values obtained from STScI
        # these values must be consistent with the given LSFs
        pixel_scale_dict = {'G130M': 9.97 / 1000. * u.AA,
                    'G160M': 12.23 / 1000. * u.AA,
                    'G140L': 80.3 / 1000. * u.AA,
                    'G230L': 390. / 1000. * u.AA,
                    'G185M': 37. / 1000. * u.AA,
                    'G225M': 33. / 1000. * u.AA,
                    'G285M': 40. / 1000. * u.AA}
        # define channel based on grating name
        channel_dict = {'G130M':  'FUV',
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
        
        if grating not in channel_dict.keys():
            raise NotImplementedError('Not ready for this HST/COS grating: {}'.format(grating))

        if channel_dict[grating] == 'NUV': #there is only 1 LSF file for NUV data
            file_name = 'nuv_all_lp1.txt'
        # COS
        elif channel_dict[grating] == 'FUV':
            # Use the ones corrected by scattering when possible
            # (currently, these are only available for lifetime-position 1)
            # check: http://www.stsci.edu/hst/cos/performance/spectral_resolution
            try:
                life_position = self.instr_config['life_position']
            except:
                raise SyntaxError('`life_position` keyword missing in `instr_config` dictionary.')

            if life_position not in ['1','2','3']:
                raise ValueError('HST/COS `life_position` should be either `1` or `2` or `3` (strings)')

            if life_position == '1':
                if grating == 'G140L': #use theoretical values 
                    file_name = 'fuv_G140L_lp1.txt'
                    
                elif grating == 'G130M': #use empirical values corrected by scattering
                    file_name = 'fuv_G130M_lp1_empir.txt'

                elif grating == 'G160M': #use empirical values corrected by scattering
                    file_name = 'fuv_G160M_lp1_empir.txt'
            
            elif life_position in ['2','3']:
                try:
                    cen_wave = self.instr_config['cen_wave']
                except:
                    raise SyntaxError('`cen_wave` keyword missing in `instr_config` dictionary. This should provide the central wavelength of the grating in Angstroms as a string.')
                #adjust format in cases where cen_wave is of the form: str(1230A)
                if cen_wave.endswith('A'): #adjust format
                    cen_wave = cen_wave[:-1]
                
                #filenames in this case have a well defined naming convention, and strict format.
                if life_position == '2':
                    file_name = 'fuv_{}_{}_lp2.txt'.format(grating,cen_wave)
                elif life_position == '3':
                    file_name = 'fuv_{}_{}_lp3.txt'.format(grating,cen_wave)
                else: # this should never happen
                    raise NotImplementedError('Unexpected error: please contact linetools developers!')

        else: # Wrong COS channel
            raise NotImplementedError('Not ready for the given HST/COS channel; only `NUV` and `FUV` channels allowed.')
        
        # point to the right file
        file_name = lt_path + '/data/lsf/{}/{}'.format(self.name,file_name)
        
        # get column names
        f = open(file_name,'r')
        line = f.readline()  # first line of file
        f.close()
        # get rid of '\n' in first line
        line = line.split('\n')[0]
        # by construction first column should be separated by `,`
        col_names = line.split(',')
        col_names[0] = 'rel_pix'
        
        pixel_scale = pixel_scale_dict[grating]  # read from dictionary defined above
        # read data
        data = ascii.read(file_name, data_start=1, names=col_names)
        
        return pixel_scale, data

    def load_STIS_data(self):
        """Load the right data according to `instr_config` for HST/STIS
        instrument"""

        # define pixel scales; values obtained from STScI
        # these values must be consistent with the given LSFs
        pixel_scale_dict = {
            'G140L': 0.60 * u.AA,
            'G140M': 0.05 * u.AA,
            'G230L': 1.58 * u.AA,
            'G230M': 0.09 * u.AA,
            'E140M': 'lambda/91700',
            'E140H': 'lambda/228000',
            'E230M': 'lambda/60000',
            'E230H': 'lambda/228000',
            'G230LB': 1.35 * u.AA,
            'G230MB': 0.15 * u.AA,
            'G430L': 2.73 * u.AA,
            'G430M': 0.28 * u.AA,
            'G750L': 4.92 * u.AA,
            'G750M': 0.56 * u.AA
        }
        # define channel based on grating name
        channel_dict = {
            'G140L': 'FUV-MAMA',
            'G140M': 'FUV-MAMA',
            'G230L': 'NUV-MAMA',
            'G230M': 'NUV_MAMA',
            'E140M': 'FUV-MAMA',
            'E140H': 'FUV-MAMA',
            'E230M': 'NUV-MAMA',
            'E230H': 'NUV-MAMA',
            'G230LB': 'CCD',
            'G230MB': 'CCD',
            'G430L': 'CCD',
            'G430M': 'CCD',
            'G750L': 'CCD',
            'G750M': 'CCD'
        }
        # also need slits
        available_slits = {
            'G140L': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'G140M': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'G230L': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'G230M': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'E140H': ['0.1x0.03', '0.2x0.09', '0.2x0.2', '6x0.2'],
            'E140M': ['0.1x0.03', '0.2x0.06', '0.2x0.2', '6x0.2'],
            'E230H': ['0.1x0.03', '0.1x0.09', '0.1x0.2', '6x0.2'],
            'E230M': ['0.1x0.03', '0.2x0.06', '0.2x0.2', '6x0.2'],
            'G430L': ['52x0.1', '52x0.2', '52x0.5', '52x2.0'],
            'G750L': ['52x0.1', '52x0.2', '52x0.5', '52x2.0']
        }

        try:
            grating = self.instr_config['grating']
        except:
            raise SyntaxError('`grating` keyword missing in `instr_config` dictionary.')
        if grating not in channel_dict.keys():
            raise NotImplementedError('Not ready for this HST/STIS grating: {}. '
                                      'Available gratings for HST/STIS are: {}'.format(grating, channel_dict.keys()))
        if grating in ['E140H', 'E140M', 'E230M', 'E230H', 'G750M', 'G430M', 'G430L']:
            raise NotImplementedError('{} not implemented yet; coming soon...'.format(grating))

        # We also need to know the slit width
        try:
            slit = self.instr_config['slit']
        except:
            raise SyntaxError('`slit` keyword missing in `instr_config` dictionary.')
        if slit not in available_slits[grating]:
            raise NotImplementedError('Not ready for this HST/STIS slit: {}. '
                                      'Available slits for HST/STIS grating {} are: {}'.format(slit, grating, available_slits[grating]))

        # now we need to read the right files
        lsf_files = glob.glob(lt_path + '/data/lsf/STIS/stis_LSF_{}_????.txt'.format(grating))
        # figure relevant wavelengths from file names
        wa_names = [fname.split('/')[-1].split('_')[-1].split('.')[0] for fname in lsf_files]
        # sort them
        sorted_inds = np.argsort(wa_names)
        lsf_files = np.array(lsf_files)[sorted_inds]
        wa_names = np.array(wa_names)[sorted_inds]

        # read the relevant kernels; they may have different rel_pix values depending on wave
        kernels_dict = dict()
        for ii, file_name in enumerate(lsf_files):
            # get the column names
            f = open(file_name,'r')
            lines = f.readlines()  # get all lines
            f.close()
            col_names = lines[1] # column names are in second line
            # get rid of '\n'
            col_names = col_names.split('\n')[0]
            # split by blank space(s) and remove the first element
            col_names = col_names.split()[1:]
            # rename new first column
            col_names[0] = 'rel_pix'

            # get original data
            data_aux = ascii.read(file_name, data_start=2, names=col_names)
            # reformat rel_pix
            data_aux['rel_pix'] = self.check_and_reformat_relpix(data_aux['rel_pix'])
            # create column with absolute wavelength based on pixel scales
            wave_aux = float(wa_names[ii])*u.AA + data_aux['rel_pix']*pixel_scale_dict[grating]
            data_aux['wv'] = wave_aux  # not used for now, but may be useful with a different approach
            # create column with normalized kernel for relevant slit
            kernel = data_aux[slit] / np.sum(data_aux[slit])
            data_aux['kernel'] = kernel

            # store only rel_pix, wv, kernel as Table
            kernels_dict[wa_names[ii]] = data_aux['rel_pix', 'wv', 'kernel']

        # at this point the kernels_dict have the original information as given by the STScI
        # tables, in their various formats...

        # project to a single rel_pix scale using interpolation; for simplicity use a
        # custom rel_pix because STScI sometimes provide them as non-constant pixel fractions.
        # This new rel_pix grid will apply to all the wavelengths

        rel_pix = kernels_dict[wa_names[0]]['rel_pix']  # use the first one as reference
        rel_pix = np.linspace(-1*np.max(rel_pix), np.max(rel_pix), len(rel_pix))  # impose them linear and symmetric
        rel_pix = self.check_and_reformat_relpix(rel_pix)  # again just in case something odd happened
        data_table = Table()
        data_table['rel_pix'] = rel_pix
        for wa_name in wa_names:
            kernel_aux = interp_Akima(data_table['rel_pix'],
                                      kernels_dict[wa_name]['rel_pix'], kernels_dict[wa_name]['kernel'])
            data_table['{}A'.format(wa_name)] = kernel_aux

        pixel_scale = pixel_scale_dict[grating]  # read from dictionary defined above
        # import pdb; pdb.set_trace()
        # todo: work out a cleverer approach to this whole issue of having different rel_pix, pixel_scales, etc
        return pixel_scale, data_table

    def interpolate_to_wv0(self, wv0):
        """Retrieves a unique LSF valid at wavelength wv0

        This is done by linearly interpolating from tabulated values
        at different wavelengths. These tabulated values (stored
        internally in self._data) are usually given as calibration
        products by instrument developers and should be loaded by
        self.load_XX_data() in the initialization stage of LSF(),
        where XX is the name of the instrument)

        Parameters
        ----------
        wv0 : Quantity 
            Wavelength at which an LSF solution is required

        Returns
        -------
        lsf_table : Table
            The interpolated lsf at wv0. This table has two 
            columns: 'wv' and 'kernel'
        """
        # get wa0 to Angstroms
        wv0 = wv0.to('AA').value

        # transform to wavelength in float() form assuming Angstroms
        col_names = self._data.keys()
        col_waves = np.array([float(name.split('A')[0]) for name in col_names[1:]])

        # find out the closest 2 columns in self._data to wv0 (on each side); these kernels will be used for interpolation
        seps =  np.fabs(col_waves - wv0)
        closest_ind = np.argsort(seps)[0]  # the 1 closest
        if closest_ind == 0:  # lower edge
            ind_blue = closest_ind
        elif closest_ind == len(col_waves) - 1:  # upper edge
            ind_blue = closest_ind - 1
        elif (col_waves[closest_ind] - wv0) > 0:  # middle case 1
            ind_blue = closest_ind - 1
        else:  # middle case 2
            ind_blue = closest_ind
        ind_red = ind_blue + 1

        # create a smaller version of self._data with the 2 most relevant columns
        good_keys = col_names[1:]  # get rid of the first name, i.e. 'rel_pix'
        good_keys = good_keys[ind_blue : ind_red + 1]
        data_aux = self._data[good_keys]
        col_waves_aux = col_waves[ind_blue : ind_red + 1]
        lsf_vals = []
        for row in data_aux:
            aux_val = []
            for i in range(0, len(row)):
                aux_val += [row[i]]

            # we don't want to extrapolate wildly, but allow LSF instantiations for wv0 outside range of 'col_waves'
            if (wv0 >= col_waves_aux[0]) & (wv0 <= col_waves_aux[-1]):
                f = interp1d(col_waves_aux,aux_val,bounds_error=True,kind='linear')  # no need to extrapolate
                lsf_vals += [f(wv0)]

            elif (wv0 < col_waves[0]) & ((col_waves[0] - wv0) < np.abs(col_waves[1] - col_waves[0])):
                f = interp1d(col_waves_aux, aux_val, bounds_error=False,
                             fill_value=aux_val[0], kind='linear')  # assign shortest wv LSF definition
                lsf_vals += [f(wv0)]

                # warning
                if (col_waves[0] - wv0) > (np.abs(col_waves[1] - col_waves[0])/2.):
                    warnings.warn(
                        "LSF may result from extrapolation outside wavelength range characterized for current grating.")

            elif (wv0 > col_waves[-1]) & ((wv0 - col_waves[-1]) < np.abs(col_waves[-1] - col_waves[-2])):
                f = interp1d(col_waves_aux, aux_val, bounds_error=False,
                             fill_value=aux_val[-1], kind='linear') #assign longest wv LSF definition
                lsf_vals += [f(wv0)]

                # warning
                if (wv0 - col_waves[-1]) > (np.abs(col_waves[-1] - col_waves[-2])/2.):
                    warnings.warn(
                        "LSF may result from extrapolation outside wavelength range characterized for current grating.")
            else:
                raise ValueError("wv0 too far outside range of defined LSFs. Perhaps you've chosen the wrong grating?")
        lsf_vals = np.array(lsf_vals)

        # normalize
        lsf_vals /= np.sum(lsf_vals)

        # create Column to store the interpolated LSF
        # lsf_vals = Column(name='{:.0f}A'.format(wv0.value),data=lsf_vals)
        lsf_vals = Column(name='kernel', data=lsf_vals)

        # create column of relative pixel in absolute wavelength
        wv_array = [(self.pixel_scale * self._data['rel_pix'][i] + wv0*u.AA).value for i in range(len(self._data))]
        wv = Column(name='wv',data=wv_array, unit=u.AA)

        # create lsf Table
        lsf = Table()
        lsf.add_column(wv)
        lsf.add_column(lsf_vals)

        # return lsf Table()
        return lsf

    def interpolate_to_wv_array(self, wv_array, kind='Akima', debug=False):
        """ Interpolate an LSF to a wavelength array.
        
        Given `wv_array` this function interpolates an LSF
        to match both scale and extent of `wv_array` using the 
        Akima or cubic-spline interpolators (default is Akima). 
        Some checks are performed too.

        Parameters
        ----------
        wv_array : Quantity numpy.ndarray, shape(N,)
            Wavelength array for which the LSF kernel is defined. The 
            central wavelength value of `wv_array` define the wavelength
            at which the LSF is defined, while the limits of `wv_array` 
            define the extent of the kernel.
        kind : str, optional
            Specifies the kind of interpolation as a string either 
            ('cubic', 'Akima'); default is `Akima`.

        Returns
        -------
        lsf_table : Table
            The interpolated lsf using at the central wavelength of 
            `wv_array`, using the same pixel scale as `wv_array`. 
            This table has two columns: 'wv' and 'kernel'. (lst_table['wv'] 
            is equal to `wv_array` by construction.)

        """
        # Check correct format
        if not ((isinstance(wv_array, np.ndarray)) or (isinstance(wv_array, Quantity))):
            raise SyntaxError('`wv_array` must be Quantity numpy.ndarray')
        elif len(wv_array.shape) != 1:
            raise SyntaxError('`wv_array` must be of shape(N,), i.e. 1-dimensional array')
        if kind not in ['cubic','Akima','akima']:
            raise ValueError('Only `cubic` or `Akima` interpolation available.')

        # define useful quantities
        wv_min = np.min(wv_array)
        wv_max = np.max(wv_array)
        wv0 = 0.5 * (wv_max + wv_min)

        lsf_tab = self.interpolate_to_wv0(wv0)

        # make sure the wv_array is dense enough to sample the LSF kernel
        kernel_wvmin = np.min(lsf_tab['wv']) * u.AA
        kernel_wvmax = np.max(lsf_tab['wv']) * u.AA
        cond = (wv_array >= kernel_wvmin) & (wv_array <= kernel_wvmax)
        if np.sum(cond) < 10:  # this number is somewhat arbitrary but reasonable
            raise ValueError('The input `wv_array` is undersampling the LSF kernel! Try a finer grid.')

        # convert to Angstroms
        wv_array_AA = np.array([wv.to('AA').value for wv in wv_array])
        
        # interpolate to wv_array
        if kind == 'cubic':
            f = interp1d(lsf_tab['wv'], lsf_tab['kernel'], kind='cubic', bounds_error= False, fill_value=0)
            lsf_vals =  f(wv_array_AA)
        elif kind in ('Akima','akima'):
            # f = Akima1DInterpolator(lsf_tab['wv'],lsf_tab['kernel']) 
            # NT: I tried Akima interpolator from scipy.interpolate
            # and is not robust in extreme situations where the
            # wv_array is large compared to the kernel FWHM.
            # Let's try linetools.analysis.interp Akima version
            lsf_vals = interp_Akima(wv_array_AA,lsf_tab['wv'],lsf_tab['kernel'])

        # make sure the kernel is never negative
        cond = lsf_vals < 0
        if np.sum(cond) > 0:
            warnings.warn('The interpolated kernel has negative values; imposing them to be 0.')
            if debug:
                import matplotlib.pyplot as plt
                plt.plot(wv_array_AA, lsf_vals, 'k-')
                # import pdb; pdb.set_trace()
            lsf_vals = np.where(lsf_vals < 0, 0., lsf_vals)

        # normalize
        lsf_vals /= np.sum(lsf_vals)

        # re-define Table
        lsf_tab = Table()
        lsf_tab.add_column(Column(name='wv', data=wv_array))
        lsf_tab.add_column(Column(name='kernel', data=lsf_vals))

        return lsf_tab
