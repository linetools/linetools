0.2 (unreleased)
----------------

Updates
.......

- Modify XSpectrum1D to use masked numpy arrays
- Enable hdf5 I/O  [requires h5py]
- Added .header property to XSpectrum1D (reads from .meta)
- Added XSpectrum1D.write, a generic write wrapper

Bug fixes
.........

- Fix `XSpectrum1D.from_tuple` to allow an astropy table column to
  specify wavelengths and fluxes.


0.1 (2016-01-31)
----------------

First public release.
