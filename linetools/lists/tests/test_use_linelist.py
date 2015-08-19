# Module to run tests on Generating a LineList
#   Also tests some simple functionality

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
	# wrest
	out = ism.all_transitions(1215.6700*u.AA)
	assert len(out) == 30 # 30 Lyman series transitions
	#print('test_all_transitions() passed')

def test_strongest_transitions():
	error_msg = 'Something is wrong in strongest_transitions()'
	ism = LineList('ISM')
	wvlims = (1200,1800)*u.AA
	z = 0.5
	transitions = ism.strongest_transitions('HI',wvlims,z,n_max=5)
	assert len(transitions) == 5, , error_msg
	assert transitions[0]['name'] == 'HI 1025', , error_msg
	assert isinstance(transitions,QTable), error_msg

	wvlims = (1500,1700)*u.AA
	z = 0.5
	transitions = ism.strongest_transitions('HI',wvlims,z,n_max=5)
	assert len(transitions) == 1, error_msg #Only Lyb should be available
	assert isinstance(transitions,dict), error_msg

	wvlims = (1100,1200)*u.AA
	z = 0.0
	transitions = ism.strongest_transitions('HI',wvlims,z,n_max=5)
	assert transitions is None, error_msg

def test_available_transitions():
	error_msg = 'Something is wrong in available_transitions()'
	ism = LineList('ISM')
	wvlims = (900,1800)*u.AA
	z = 0.1

	transitions = ism.available_transitions(wvlims,z,n_max=3,n_max_tuple=5)
	assert len(transitions) == 3, error_msg 
	assert transitions[2]['name'] == 'HI 972', , error_msg
	assert isinstance(transitions,QTable), error_msg

	transitions = ism.available_transitions(wvlims,z,n_max=10,n_max_tuple=2)
	assert len(transitions) == 10, error_msg 
	assert transitions[2]['name'] == 'CIII 977' , error_msg #Only Lya and Lyb as the first two

	wvlims = (1200,1800)*u.AA
	z = 0.5
	transitions = ism.available_transitions(wvlims,z,n_max=15,n_max_tuple=2)
	assert len(transitions) == 15, error_msg 
	assert transitions[0]['name'] == 'HI 1025', error_msg 
	assert 'OVI 1031' in transitions['name'], error_msg 
	assert 'CIII 977' in transitions['name'], error_msg

	wvlims = (1000,3000)*u.AA
	z = 1.5
	transitions = ism.available_transitions(wvlims,z,n_max=100,n_max_tuple=2)
	assert 'NeVIII 770' in transitions['name'], error_msg
	assert 'MgX 609' in transitions['name'], error_msg
	assert 'HI 1215' not in transitions['name'], error_msg

	wvlims = (1215.6,1217)*u.AA
	z = 0
	transitions = ism.available_transitions(wvlims,z,n_max=100,n_max_tuple=2)
	assert isinstance(transitions,dict), error_msg
	