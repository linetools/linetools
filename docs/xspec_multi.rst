.. _XSpectrum1D:

***********
Multi Xspec
***********

.. index:: XSpec_multi

Overview
========

`~linetools.spectra.xspectrum1d.XSpectrum1D` may contain
multiple spectra, typically taken with the same instrument
and configuration.  These docs describe several methods
related to multple spectra.

Generation
==========

There are currently two ways to generate a multi-spectrum
`~linetools.spectra.xspectrum1d.XSpectrum1D` object.  The
first is to feed it a set of np.arrays each with dimension
(nspec,npix).::

    nspec, npix = 3, 100
    wave = np.outer(np.ones(nspec), np.arange(npix))
    flux = np.ones_like(wave)
    sig = np.ones_like(wave)
    #
    mspec = XSpectrum1D(wave, flux, sig)

Alternatively, one can generate from a list of
`~linetools.spectra.xspectrum1d.XSpectrum1D` objects.::

    from linetools.spectra import utils as ltsu
    mspec = ltsu.collate([spec1,spec2])


Generation
==========
