""" Curve of Growth calculations
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import sys
import os
import warnings
import pdb

from scipy import integrate
from scipy.interpolate import interp1d

from astropy import units as u
from astropy.modeling import FittableModel, Parameter
from astropy.modeling import fitting

#from xastropy.xutils import xdebug as xdb

# Begin
def _ftau_intgrnd(x,tau0=0.1):
    return 1 - np.exp(-tau0 * np.exp(-x**2))

# Generate Ftau (could archive, but this is reasonably fast)
neval = 10000
lgt = np.linspace(-3, 9, neval)
all_tau0 = 10.**lgt
#
xFtau0 = np.zeros(neval)
for jj,tau0 in enumerate(all_tau0):
    xFtau0[jj], ferr = integrate.quad(_ftau_intgrnd, 0, np.inf, args=(tau0,))

# Now interpolate
intFtau0 = interp1d(all_tau0, xFtau0, bounds_error=False,fill_value=0.)

##############################
def cog_plot(COG_dict):
    """Generate a plot for COG solution

    Parameters
    ----------
    COG_dict : dict
      dict containing the COG inputs and solution from single_cog_analysis
    """
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'stixgeneral'
    mpl.rcParams['font.size'] = 17.
    # Plot
    plt.clf()
    ax = plt.gca()
    # Data with values
    gdv = COG_dict['redEW'] > 0.
    yerr=(COG_dict['sig_EW'][gdv]/COG_dict['wrest'][gdv])/COG_dict['redEW'][gdv]
    ax.errorbar(np.log10(COG_dict['f'][gdv]*COG_dict['wrest'][gdv].to('cm').value),
                np.log10(COG_dict['redEW'][gdv]),
                yerr=yerr, fmt='o')
    # Upper limit
    upper = COG_dict['redEW'] <= 0.
    ax.scatter(np.log10(COG_dict['f'][upper]*COG_dict['wrest'][upper].to('cm').value),
               np.log10(3*COG_dict['sig_EW'][upper]/COG_dict['wrest'][upper]), color='red', marker='v')
    # Model
    xval = np.log10(COG_dict['f']*COG_dict['wrest'].to('cm').value)
    xmod = np.linspace(np.min(xval), np.max(xval), 200)
    tau0 = 1.497e-15*(10**(xmod+8))*(10.**COG_dict['logN'])/COG_dict['b'].to('km/s').value
    Ftau0 = intFtau0(tau0)
    ymod = np.log10(2*COG_dict['b'].to('km/s').value*Ftau0/3e5)
    #pdb.set_trace()
    ax.plot(xmod,ymod,'g--')
    # Axes
    ax.set_xlabel(r'$\log_{10} \, (f \, \lambda)$')
    ax.set_ylabel(r'$\log_{10} \, (W / \lambda)$')
    # Label
    ax.text(0.1, 0.9, r'$\log N = ${:.2f}$\pm${:.2f}'.format(COG_dict['logN'],COG_dict['sig_logN']), transform=ax.transAxes, ha='left', va='center', fontsize='large')#, bbox={'facecolor':'white'})
    ax.text(0.1, 0.8, r'$b = ${:.1f}$\pm${:.2f}'.format(COG_dict['b'].value,COG_dict['sig_b']), transform=ax.transAxes, ha='left', va='center', fontsize='large')#, bbox={'facecolor':'white'})
    # Finish
    plt.show()


def single_cog_analysis(wrest, f, EW, sig_EW=None, guesses=None):
    """Perform COG analysis on a single component

    Parameters
    ----------
    wrest : Quantity array
      Rest wavelengths
    f : float array
      f-values
    EW : Quantity array
      Measured EWs
    sig_EW : Quantity array, optional
      Measured sig_EWs
    guesses : tuple of float,float
      Guesses for logN, b

    Returns
    -------
    COG_dict : dict
      dict containing inputs and solution, e.g.
       * logN : float
       * b : Quantity
           Doppler parameter (km/s)
    """
    if guesses is None:
        logN=14.
        b=10.
    else:
        logN = guesses[0]
        b = guesses[1]
    # Reduced EW
    redEW = (EW / wrest).value
    # Weights
    if sig_EW is not None:
        weights = (wrest/sig_EW)**2
    # COG model
    cog_model = single_cog_model(logN=logN, b=b)
    # Fitter
    fitter = fitting.LevMarLSQFitter()
    # Fit
    parm = fitter(cog_model, wrest.to('AA').value*f, redEW, weights=weights)
    # Error
    covar = fitter.fit_info['param_cov']
    sigma = np.zeros(covar.shape[0])
    for ii in range(sigma.size):
        sigma[ii] = np.sqrt(covar[ii,ii])
    # Generate COG dict
    COG_dict = dict(wrest=wrest,f=f,EW=EW,sig_EW=sig_EW,
                    redEW=redEW,
                    logN=parm.logN.value, sig_logN=sigma[0],
                    b=parm.b.value*u.km/u.s, sig_b=sigma[1]*u.km/u.s,
                    parm=parm)
    # Return
    return COG_dict

class single_cog_model(FittableModel):
    """Generate a single COG model

    Parameters
    ----------
    logN
    b

    input : wrest*f
    output : redEW
      reduced EWs
    """
    inputs = ('wrestxf',)
    outputs = ('redEW',)

    # Free parameters (generally)
    logN=Parameter()
    b=Parameter()  # Assumes km/s

    # Fixed parameters

    @staticmethod
    def evaluate(wrestf,logN,b):
        # F(tau0)
        tau0 = 1.497e-15*(wrestf)*(10.**logN)/b
        Ftau0 = intFtau0(tau0)
        # Finish
        redEW = 2*b*Ftau0/3e5
        return redEW
