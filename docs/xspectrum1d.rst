.. _XSpectrum1D:

*****************
XSpectrum1D Class
*****************

.. index:: XSpectrum1D

Overview
========

`~linetools.spectra.xspectrum1d.XSpectrum1D` contains and
manipulates 1-d spectra, each of
which usually consists of a wavelength, flux and flux uncertainty
array.  For absorption-line analysis, it also often contains a
continuum array.

The data are held in a masked numpy array which
may contain multiple spectra.
By default pixels on the edges of the
spectrum with an input error array having values <=0 are masked
at instantiation. It is important to appreciate this masking.
It does mean that you will not view, print, analyze, etc. pixels
that have been masked.

Attributes
==========

The main attributes of XSpectrum1D are `wavelength`, `flux` and
`sig`. Let's begin by creating a spectrum using the
`~linetools.spectra.xspectrum1d.XSpectrum1D.from_tuple` method::

    >>> from linetools.spectra.xspectrum1d import XSpectrum1D
    >>> import numpy as np
    >>> wa = np.arange(3000, 7000.1, 0.5)
    >>> fl = np.ones_like(wa)
    >>> sig = np.ones_like(fl) * 0.1
    >>> sp = XSpectrum1D.from_tuple((wa, fl, sig), verbose=False)
    >>> sp.wavelength # doctest: +SKIP
    <Quantity [ 3000. , 3000.5, 3001. ,..., 6999. , 6999.5, 7000. ] Angstrom>
    >>> sp.flux # doctest: +SKIP
    <Quantity [ 1., 1., 1.,...,  1., 1., 1.]>
    >>> sp.sig # doctest: +SKIP
    <Quantity [ 1., 1., 1.,...,  1., 1., 1.]>

Note that all three arrays have units. If you don't
specify a unit when you create an new XSpectrum1D instance, Angstroms
are assumed for wavelength and dimensionless_unscaled
for flux. The 1-sigma uncertainty is always assumed to have the
same units as the flux. All of these are specified in the sp.units dict.

If one loads multiple 1D spectra (e.g. a brick of data from DESI
or a set of spectra from
`specdb <https://github.com/specdb/specdb>`_),
the selected spectrum is given by the spec.select index.

All of the values are stored in the masked spec.data numpy array
with columns `wave`, `flux`, `sig`, and `co` (the latter is
for a continuum).

Init
====

Reading
-------

Read spectra from a file using ``XSpectrum1D.from_file``, which uses the same
syntax as `~linetools.spectra.io.readspec`.  See
below for a complete listing of permitted file formats.

The easiest way to create
a new spectrum from a set of data arrays for a single
spectrum is to use ``sp.from_tuple`` as shown above.
Here are a series of example calls to generate the class::

    sp = XSpectrum1D.from_file('PH957_f.fits')      # From a FITS file
    sp = XSpectrum1D.from_file('q0002m422.txt.gz')  # From an ASCII table
    sp = xspec1.copy()                              # From an XSpectrum1D object
    sp = XSpectrum1D.from_tuple((wa, fl, sig), verbose=False)



Masking
-------

The guts of XSpectrum1D is a ndarray array named data
which contains the wave, flux, sig, etc. values.  This
is a masked array which is convenient for many applications.
If you wish to view/analyze all pixels in your spectrum including
those with 0 or NAN sig values, then disable the mask when
creating the object (masking='None') or by using the unnmask() method::

    sp = XSpectrum1D.from_tuple((wa, fl, sig), masking='none')
    sp = XSpectrum1D.from_file('PH957_f.fits')
    sp.unmask()

Methods
=======

Writing
-------

There are a number of methods to write a file, e.g.
`sp.write_to_fits`. FITS files are preferable because they are
generally faster to read and write, require less space, and
are generally easier for other software to read.
Another option is an HDF5 file which better preserves the
data format of XSpectrum1D.  Here are some examples::

    sp.write_to_fits('QSO.fits')            # Standard FITS file
    sp.write('QSO.fits')                    # Same
    sp.write('QSO.fits', FITS_TABLE=True)   # Binary FITS table
    sp.write_to_hdf5('QSO.hdf5')            # HDF5 file
    sp.write('QSO.hdf5')                    # Same
    sp.write_to_ascii('QSO.ascii')          # ASCII (heaven forbid)
    sp.write('QSO.ascii')                   # Same


One can collate a list of XSpectrum1D objects into one with collate::

    sp1 = XSpectrum1D.from_file('PH957_f.fits')
    sp2 = XSpectrum1D.from_file('q0002m422.txt.gz')
    sp = linetools.spectra.utils.collate([sp1,sp2])


Plotting
--------

`sp.plot()` plots the spectrum, which you can then navigate around
using the same keys as `~lt_xspec` (as well as the usual matplotlib
navigation tools).
**Note**:  if you are using MacOSX then you will
probably need to change your *backend* from macosx to TkAgg
in the matplotlibrc file.

Rebinning
---------

`~linetools.spectra.xspectrum1d.XSpectrum1D.rebin` rebins the spectrum
to an arbitrary input wavelength array.  Flux is conserved.  If
*do_sig=True*, the error array is rebinned as well and a crude attempt
is made to conserve S/N.  Generally, neighboring pixels will be
correlated::

    newspec = sp.rebin(new_wv, do_sig=True)

If the XSpectrum1D object containts multiple spectra, you can rebin
all of them to the new wavelength array as well::

    newspec = sp.rebin(new_wv, do_sig=True, all=True)


Continuum fitting
-----------------

`~linetools.spectra.xspectrum1d.XSpectrum1D.fit_continuum` enables you
to interactively fit a continuum to the spectrum. Currently it's
optimised to estimate the continuum for high-resolution quasar
spectra, but it should be applicable to any spectrum with a slowly
varying continuum level and narrow absorption features. Once a
continuum has been fitted, it can be accessed using the `co`
attribute. The spectrum can also be normalised (i.e the flux values
returned by spec.flux are divided by the continuum) with the
`~linetools.spectra.xspectrum1d.XSpectrum1D.normalize`
method.  This also sets spec.normed to True.

Finally, you can apply small variations to the continuum
anchor points with
`~linetools.spectra.xspectrum1d.XSpectrum1D.perturb_continuum` to see
how changes in the continuum level affect your analysis.

Smoothing
---------

There are several algorithms included that smooth the
input spectrum and return a new XSpectrum1D.  These are
`~linetools.spectra.xspectrum1d.XSpectrum1D.box_smooth`,
`~linetools.spectra.xspectrum1d.XSpectrum1D.gauss_smooth`,
and
`~linetools.spectra.xspectrum1d.XSpectrum1D.ivar_smooth`.

Other methods
-------------

You can join one XSpectrum1D instance with another overlapping
spectrum using `~linetools.spectra.xspectrum1d.XSpectrum1D.splice`.
`~linetools.spectra.xspectrum1d.XSpectrum1D.pix_minmax` finds the
pixel indices corresponding to a wavelength or velocity range, and
`~linetools.spectra.xspectrum1d.XSpectrum1D.add_noise` adds noise to
the spectrum. We have also implemented a method that estimates a local
average signal-to-noise ratio at a given observed wavelength
(`~linetools.spectra.xspectrum1d.XSpectrum1D.get_local_s2n`), which is capable
of masking out pixels that are below a flux threshold (useful for excluding
strong absorption features from the calculation). For a complete list of
all the available methods, see the API: `~linetools.spectra.xspectrum1d.XSpectrum1D`.

Multi-spec methods
------------------

See :doc:`xspec_multi` for more.

File Formats Read
=================

Below is a table of the types of spectra files that can be read by
`~linetools.spectra.io.readspec`.  If your file cannot be read, please
open an issue on the `linetools issue tracker
<http://github.com/linetools/linetools/issues>`_.

========================================================== =================
Description                                                Instruments
========================================================== =================
simple 1D FITS files                                       ESI, HIRES, etc.
binary FITS table from LowRedux                            LRIS,Kast,etc.
multi-extension 1D FITS files from LowRedux                LRIS,Kast,etc.
binary FITS tables from many other sources                 COS, SDSS, etc.
multi-extension binary FITS tables from PYPIT              LRIS,Kast,etc.
brick files (2D images: flux, ivar; 1D image: wavelength)  DESI
`UVES_popler`_ output files                                UVES
========================================================== =================

.. _UVES_popler: http://astronomy.swin.edu.au/~mmurphy/UVES_popler/
