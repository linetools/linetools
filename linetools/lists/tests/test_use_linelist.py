# Module to run tests on Generating a LineList
#   Also tests some simple functionality

from __future__ import (print_function, absolute_import, division,
unicode_literals)

# TEST_UNICODE_LITERALS

import os, pdb
import pytest
import astropy.io.ascii as ascii
from astropy import units as u
from astropy.table import QTable
import numpy as np

from linetools.lists.linelist import LineList

#import pdb
#pdb.set_trace()

# ISM LineList
def test_lines_from_ion():
    ism = LineList('ISM')
    # 
    lines = ism[(6,2)]
    assert (1334.5323*u.AA in lines['wrest'])

def test_subset():
    ism = LineList('ISM')
    subset = np.array([1215.6700, 1608.4511])*u.AA
    #pytest.set_trace()
    ism = ism.subset_lines(subset)
    assert len(ism._data) == 2
    np.testing.assert_allclose(ism['FeII 1608']['wrest'], 1608.4511*u.AA, rtol=1e-7)

    # Now with names
    ism = LineList('ISM')
    subset = ['HI 1215', 'HI 1025', 'CIV 1548']
    ism = ism.subset_lines(subset)
    np.testing.assert_allclose(ism['HI 1215']['wrest'], 1215.6700*u.AA, rtol=1e-7)

def test_closest():
    ism = LineList('ISM')
    ism.closest=True
    # 
    line = ism[1250.584*u.AA]
    np.testing.assert_allclose(line['wrest'], 1250.578*u.AA, rtol=1e-7)

def test_all_transitions():
    error_msg = 'Something is wrong in all_transitions()'
    ism = LineList('ISM')
    #check simple case
    line = 'OVI'
    ovi_transitions = ism.all_transitions(line)
    assert len(ovi_transitions) == 2, error_msg
    #print(ovi_transitions['name'])
    #check unknown
    line = 'unknown'
    out = ism.all_transitions(line)
    assert type(out) == dict, error_msg
    #check case of single transition ion
    line = 'CIII'
    out = ism.all_transitions(line)
    assert type(out) == dict, error_msg
    #check case of transitions from excited levels
    line='FeII*'
    out = ism.all_transitions(line)
    assert len(out) == 27, "wrong line counts"
    print(out)
    # wrest
    out = ism.all_transitions(1215.6700*u.AA)
    assert len(out) == 30,"wrong line counts" # 30 Lyman series transitions
    #print('test_all_transitions() passed')
    h2 = LineList('H2')
    line = 'B19-0P(1)'
    out = h2.all_transitions(line)
    assert len(out) == 7


def test_strongest_transitions():
    error_msg = 'Something is wrong in strongest_transitions()'
    ism = LineList('ISM')
    wvlims = (1200,1800)*u.AA
    z = 0.5
    transitions = ism.strongest_transitions('HI',wvlims/(1+z),n_max=5)
    assert len(transitions) == 5,  error_msg
    assert transitions[0]['name'] == 'HI 1025' , error_msg
    assert isinstance(transitions,QTable), error_msg

    wvlims = (1500,1700)*u.AA
    z = 0.5
    transitions = ism.strongest_transitions('HI',wvlims/(1+z),n_max=5)
    assert isinstance(transitions,dict), error_msg #only Lyb should be available, so dict is expected
    assert transitions['name'] == 'HI 1025'

    wvlims = (1100,1200)*u.AA
    z = 0.0
    transitions = ism.strongest_transitions('HI',wvlims/(1+z),n_max=5)
    assert transitions is None, error_msg

def test_available_transitions():
    error_msg = 'Something is wrong in available_transitions()'
    ism = LineList('ISM')
    wvlims = (900,1800)*u.AA
    z = 0.1

    transitions = ism.available_transitions(wvlims/(1+z),n_max_tuple=5)
    assert transitions[2]['name'] == 'HI 972' , error_msg
    assert isinstance(transitions,QTable), error_msg

    transitions = ism.available_transitions(wvlims/(1+z),n_max_tuple=2)
    assert transitions[2]['name'] == 'CIII 977' , error_msg

    wvlims = (1200,1800)*u.AA
    z = 0.5
    transitions = ism.available_transitions(wvlims/(1+z), n_max_tuple=2)
    assert transitions[0]['name'] == 'HI 1025', error_msg 
    assert 'OVI 1031' in transitions['name'], error_msg 
    assert 'CIII 977' in transitions['name'], error_msg

    wvlims = (1000,3000)*u.AA
    z = 1.5
    transitions = ism.available_transitions(wvlims/(1+z),n_max_tuple=2)
    assert 'NeVIII 770' in transitions['name'], error_msg
    assert 'MgX 609' in transitions['name'], error_msg
    assert 'HI 1215' not in transitions['name'], error_msg

    wvlims = (1215.6,1217)*u.AA
    z = 0
    transitions = ism.available_transitions(wvlims/(1+z),n_max_tuple=2)
    assert isinstance(transitions,dict), error_msg

def test_sortdata():
    error_msg = 'Something is wrong in sortdata()'
    ism = LineList('ISM', sort_by='name')
    assert ism.name[0] == 'AlII 1670', error_msg
    ism.sortdata('name', reverse=True)
    assert ism.name[0] == 'ZrIII 1798', error_msg
    ism.sortdata(['abundance', 'rel_strength'], reverse=True)
    assert ism.name[0] == 'HI 1215', error_msg
    ism.sortdata(['rel_strength'])
    assert ism.name[0] == 'CI** 1123b', error_msg
