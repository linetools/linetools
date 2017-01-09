.. _xspec_multi:

*****************
Multi XSpectrum1D
*****************

.. index:: xspec_multi

Overview
========

`~linetools.spectra.xspectrum1d.XSpectrum1D` may contain
multiple spectra, typically taken with the same instrument
and configuration.  These docs describe several methods
related to multple spectra.

Generation
==========

There are currently two ways to generate a multi-spectrum
`~linetools.spectra.xspectrum1d.XSpectrum1D` object.

Instantiation
-------------

The first is to feed it a set of np.arrays each with dimension
(nspec,npix).::

    nspec, npix = 3, 100
    wave = np.outer(np.ones(nspec), np.arange(npix))
    flux = np.ones_like(wave)
    sig = np.ones_like(wave)
    #
    mspec = XSpectrum1D(wave, flux, sig)

Collate
-------

Alternatively, one can generate from a list of
`~linetools.spectra.xspectrum1d.XSpectrum1D` objects
using the collate() method.::

    from linetools.spectra import utils as ltsu
    mspec = ltsu.collate([spec1,spec2])

If any of the input `~linetools.spectra.xspectrum1d.XSpectrum1D` objects
are already multi-dimensional, these will all be ingested.

Note: if any of the input objects has a continuum, the resultant
object will provide a continuum for all with default value=0.
These 0. values are ignored even if the spectrum is later normalized.

Slice
=====

Return a sliced portion of the XSpectrum1D object.  Indices
can be repeated.  A couple of examples::

    two_spec = mspec[np.array([0,1])]
    one_spec = mspec[1]
    two_spec = mspec[0:2]

Note that the full set of headers are maintained, i.e.
these are not sliced.


Rebin to Rest
=============

By inputting an array of redshifts and a velocity
width, one can rebin the multi-spec to a common
rest-frame spectrum with constant dv pixels.::

    rest_spec = ltsu.rebin_to_rest(mspec, zarr, 100*u.km/u.s)

The output is a new multi-spec object with a common
rest-frame wavelength array.

Smash(stack)
============

Smash down a multi-spec object into a 1D spectrum.::

    stack = ltsu.smash_spectra(rest_spec, method='average')


