# Module to run tests on initializing AbsLine

# TEST_UNICODE_LITERALS

import numpy as np
import os, pdb
import pytest
from astropy import units as u

from linetools.spectralline import AbsLine, AbsComponent

def test_mk_abscomponent():
	# Init HI Lya
    lya = AbsLine(1215.6700*u.AA)
    # Init HI Lyb
    lyb = AbsLine('HI 1025')
    # Init HI Lyc, etc
    lyc = AbsLine('HI 972')
    lyd = AbsLine('HI 949')
    lye = AbsLine('HI 937')
    lyf = AbsLine('HI 930')
    
    #Init component in different ways
    component = AbsComponent(lya)
    component = AbsComponent([lya])
    component = AbsComponent([lya,lyb,lyc])
    
    #test adding AbsLines
    component = AbsComponent(lya)
    component.add_abslines(lyb)
    component.add_abslines([lyc,lyd])
    component.add_abslines([lye])
    
    #remove AbsLines
    component.remove_absline('HI 1025')
    component.add_abslines(lyb)
    component.remove_abslines(['HI 1025'])
    component.add_abslines(lyb)
    component.remove_abslines(['HI 1025','HI 1215'])