from astropy import units as u

from linetools.analysis import plots as ltap
from linetools.spectralline import AbsLine
from linetools.spectra import io as ltsio

import os

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), '../../spectra/tests/files')
    return os.path.join(data_dir, filename)

def test_stack_plot(show=False):
    abslin1 = AbsLine(1548.195*u.AA)
    abslin2 = AbsLine('CIV 1550')
    # no spectrum first
    ltap.stack_plot([abslin1], show=show)
    # Set spectrum
    spec = ltsio.readspec(data_path('UM184_nF.fits'))  # already normalized
    abslin1.analy['spec'] = spec
    abslin1.analy['wvlim'] = [6079.78, 6168.82]*u.AA
    abslin1.attrib['z'] = 2.92929
    ltap.stack_plot([abslin1], show=show)
    # second line
    abslin2.analy['spec'] = spec
    abslin2.analy['wvlim'] = [6079.78, 6168.82]*u.AA
    abslin2.attrib['z'] = 2.92929
    ltap.stack_plot([abslin1, abslin2], show=show)

