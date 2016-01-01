.. _XSpectrum1D:

*****************
XSpectrum1D Class
*****************

.. index:: XSpectrum1D

Notebooks
=========


Overview
========

This Class describes a spectrum. It subclasses the specutils
Spectrum1D class, adding several new attributes and methods.


Methods
=======

rebin
-----

This method rebins the spectrum to an arbitrary input wavelength array.
Flux is conserved.  If *do_sig=True*, the error array is rebinned as well
and a crude attempt is made to conserve S/N.  Generally, neighboring
pixels will be correlated. ::

    newspec = spec.rebin(new_wv, do_sig=True)

