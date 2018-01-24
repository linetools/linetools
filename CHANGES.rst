0.3 (Unreleased)
----------------

Updates
.......
- Added extra attributes to AbsComponents
- Added script to get HST/COS life-time position from date
- Added Cashman+2017 catalog in LineList
- Significant refactor of AbsComponent
- LineList "AGN" added

Bug fixes
.........


0.2 (2017-10-24)
----------------

Updates
.......
- Extra features to main Objects like XSpectrum1D, AbsComponent, AbsLine, LineList
- Added some extra emisison lines to Galaxy LineList
- Refactor from pyQt4 -> pyQt5
- Improvements to GUIs and scripts
- Added EmLine and EmSystem classes
- LineList.available_transitions() no longer has key argument n_max
- LineList: extra attributes for transitions added (`ion_name`, `log(w*f)`, `abundance`, `ion_correction`, `rel_strength`)
- ASCII tables with no header are required to be 4 columns or less for `io.readspec` to work
- Modify XSpectrum1D to use masked numpy arrays
- Enable hdf5 I/O [requires h5py]
- Added .header property to XSpectrum1D (reads from .meta)
- Added XSpectrum1D.write, a generic write wrapper
- Added XSpectrum1D.get_local_s2n, a method to calculate average signal-to-noise at a given wavelength
- Added xabssysgui GUI
- Added new LineList(e.g. Galaxy)
- Added EmLine child to SpectralLine
- Added LineLimits class
- Added SolarAbund class
- Added lt_radec and lt_line scripts
- Added LSF class to handle line-spread-functions. Currently implemented for HST/COS and most HST/STIS configurations.

Bug fixes
.........

- Fix `XSpectrum1D.from_tuple` to allow an astropy table column to
  specify wavelengths and fluxes.


0.1 (2016-01-31)
----------------

First public release.
