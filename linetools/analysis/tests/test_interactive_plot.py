from __future__ import print_function, absolute_import, division, \
     unicode_literals

import linetools.analysis.interactive_plot as inter 

import matplotlib.pyplot as plt
import numpy as np
import os

def test_PlotWrapNav():
    
    wa = np.linspace(1, 10, 50)
    fl = np.ones_like(wa)

    fig = plt.gcf()
    ax = plt.gca()
    artists = {}
    artists['fl'] = ax.plot(wa, fl)[0]

    inter.PlotWrapNav(fig, ax, wa, fl, artists, printhelp=False)

def test_InteractiveCoFit():

    wa = np.linspace(1, 10, 50)
    fl = np.ones_like(wa)
    er = np.ones_like(wa) * 0.01
    contpoints = [(1,1), (5, 1.2), (9, 0.9)]
    if os.path.exists('_knots.jsn'):
        os.remove('_knots.jsn')
    inter.InteractiveCoFit(wa, fl, er, contpoints, co=None,
                     fig=None)

    
