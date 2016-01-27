.. _XSpectrum1D:

*****************
XSpectrum1D Class
*****************

.. index:: XSpectrum1D

Overview
========

`~linetools.spectra.xspectrum1d.XSpectrum1D` describes a 1-d spectrum,
which usually consists of a wavelength, flux and flux uncertainty
array. It subclasses the specutils Spectrum1D class, adding several
new attributes and methods.

Attributes
==========

Its main attributes are the `wavelength`, `flux` and
`uncertainty`. Let's create a spectrum using the
`~linetools.spectra.xspectrum1d.XSpectrum1D.from_tuple` method::

    >>> from linetools.spectra.xspectrum1d import XSpectrum1D
    >>> import numpy as np
    >>> wa = np.arange(3000, 7000.1, 0.5)
    >>> fl = np.ones_like(wa)
    >>> sig = np.ones_like(fl) * 0.1
    >>> sp = XSpectrum1D.from_tuple((wa, fl, sig))
    >>> sp.wavelength # doctest: +SKIP
    <Quantity [ 3000. , 3000.5, 3001. ,..., 6999. , 6999.5, 7000. ] Angstrom>
    >>> sp.flux
    <Quantity [ 1., 1., 1.,...,  1., 1., 1.]>
    >>> sp.uncertainty # doctest: +SKIP
    <astropy.nddata.nduncertainty.StdDevUncertainty at 0x10bec8e80>
   
Note that the wavelength and flux both have units. If you don't
specify a unit when you create an new XSpectrum1D instance, Angstroms
are assumed. In this case the flux is unitless. The one sigma
uncertainty is assumed to have the same units as the flux, and you can
access its underlying numpy array via ``sp.uncertainty.array``.

Methods
=======

Reading and Writing
-------------------

Read spectra from a file using ``sp.from_file``, which uses the same
syntax as `~linetools.spectra.io.readspec`. The easiest way to create
a new spectrum from data arrays is to use ``sp.from_tuple`` as shown
above.

To write a spectrum toa file, use either `sp.write_to_fits` or
`sp.write_to_ascii`. FITS files are preferable because they are
generally faster to read and write.

Plotting
--------

`sp.plot()` plots the spectrum, which you can then navigate around
using the same keys as `~lt_xspec` (as well as the usual matplotlib
navigation tools).

Rebinning
---------

`~linetools.spectra.xspectrum1d.XSpectrum1D.rebin` rebins the spectrum
to an arbitrary input wavelength array.  Flux is conserved.  If
*do_sig=True*, the error array is rebinned as well and a crude attempt
is made to conserve S/N.  Generally, neighboring pixels will be
correlated::

    >>> newspec = sp.rebin(new_wv, do_sig=True) # doctest: +SKIP


Continuum fitting
-----------------

`~linetools.spectra.xspectrum1d.XSpectrum1D.fit_continuum` enables you
to interactively fit a continuum to the spectrum. Currently it's
optimised to estimate the continuum for high-resolution quasar
spectra, but it should be applicable to any spectrum with a slowly
varying continuum level and narrow absorption features. Once a
continuum has been fitted, it can be accessed using under the `co`
attribute. The spectrum can also be normalised (i.e the flux is
divided by the continuum) with the
`~linetools.spectra.xspectrum1d.XSpectrum1D.normalize`
method. Finally, you can apply small variations to the continuum
anchor points with
`~linetools.spectra.xspectrum1d.XSpectrum1D.perturb_continuum` to see
how changes in the continuum level affect your analysis.


Other methods
-------------

You can join one XSpectrum1D instance with another overlapping
spectrum using `~linetools.spectra.xspectrum1d.XSpectrum1D.splice`.
`~linetools.spectra.xspectrum1d.XSpectrum1D.pix_minmax` finds the
pixel indices corresponding to a wavelength or velocity range, and
`~linetools.spectra.xspectrum1d.XSpectrum1D.add_noise` add noise to
the spectrum. For a complete list of all the available methods, see
the API: `~linetools.spectra.xspectrum1d.XSpectrum1D`.
  
