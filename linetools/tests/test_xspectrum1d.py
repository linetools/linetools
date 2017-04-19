# Module to run test on XSpectrum1D

from __future__ import print_function,absolute_import,division,unicode_literals
import numpy as np
import glob,os,sys,copy
from astropy.table import QTable,table
from astropy.io import ascii
from astropy import units as u
from astropy import constants as const
from matplotlib import pyplot as plt
from enigma.qpq import spec as qpqs

datadir = os.getenv('DROPBOX_DIR') + '/QSOPairs/data/'
spec_dict = qpqs.load_spec(datadir+'ESI_redux/SDSSJ103900.01+502652.8_F.fits.gz')

spec_dict['spec'].normalize(spec_dict['conti'])
spec_dict['spec'] = spec_dict['spec'].normalized_spec()

plt.plot(spec_dict['spec'].data['flux'].compressed())
plt.show()
