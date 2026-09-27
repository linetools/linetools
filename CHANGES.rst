0.3 (Unreleased)
----------------

Updates
.......
- Dropped support for Python 3.8, 3.9 and 3.10; ``requires-python`` is now
  ``>=3.11``.  **This is a breaking change for installers on older Pythons.**
- Raised dependency floors to the oldest releases with CPython 3.11 wheels
  (numpy>=1.24, astropy>=5.3, scipy>=1.10, h5py>=3.8, matplotlib>=3.6,
  PyYAML>=6.0, IPython>=8.0, qtpy>=2.0).  The previous floors were not
  installable on 3.11 and were never exercised by CI.
- Migrated all packaging metadata from ``setup.cfg`` into ``pyproject.toml``
  (PEP 621).  ``setup.py``, ``setup.cfg``, ``old_setup.py``, ``old_setup.cfg``
  and ``.travis.yml`` have been removed.
- ``linetools.__version__`` is now available (previously ``linetools``
  exported nothing).
- Added a ``LICENSE`` file at the repository root so the BSD-3-Clause license
  is detected by GitHub and packaging tools.
- Replaced the deprecated ``astropy.utils.isiterable`` with ``numpy.iterable``
  throughout.
- Rebuilt CI: tests now run on Python 3.11-3.14, on macOS as well as Linux,
  with an ``oldestdeps`` job that pins to the declared dependency floors and a
  job that builds the docs.
- Dropped the ``codecov`` test dependency and the Coveralls badge; coverage is
  no longer uploaded to an external service.
- ``testpaths`` now covers the whole ``linetools`` package.  A bare ``pytest``
  previously collected only ``linetools/tests`` (about 30 of 232 tests) and
  silently skipped every sub-package test directory.
- Filled in the copyright holder in ``LICENSE`` / ``licenses/LICENSE.rst``,
  which still carried the unedited ``Copyright (c) year, author`` placeholder
  from the astropy package template.
- Added extra attributes to AbsComponents
- Added script to get HST/COS life-time position from date
- Added Cashman+2017 catalog in LineList
- Significant refactor of AbsComponent
- LineList "AGN" added
- Refactor from PyQt5 -> PySide2
- Refactor from PySide2 -> QtPy

Bug fixes
.........
- Fixed a Python 3 ``TypeError`` when gzip-compressing generated line-list
  tables: ``linetools.lists.parse.mktab_morton03`` and
  ``parse_verner96(write=True)`` opened the source file in text mode and fed it
  to a binary gzip handle.  Both paths had been broken for every Python 3 user.
- ``linetools.utils.compare_two_files`` used ``~`` (bitwise inversion) on a
  bool, which gave the wrong result for any ``verbose`` value other than
  ``True`` and is deprecated in Python 3.16.
- Tests no longer write scratch files into the source tree or into the
  installed package directory; they use pytest's ``tmp_path`` fixture.
- Re-enabled ``test_morton03``, which had been skipped on every Python 3 since
  the Python 2 era and was masking the gzip bug above.
- change linestyle='steps-mid' to drawstyle='steps-mid' in Axes.plot calls


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
