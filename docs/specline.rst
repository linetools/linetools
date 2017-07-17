.. _SpectralLine:

******************
SpectralLine Class
******************

.. index:: SpectralLine

Notebooks
=========

.. toctree::
   :maxdepth: 1

Overview
========

This Class is designed to organize and analyze an spectral line,
either emission or absorption. This is an abstract base class,
i.e. one must instantiate one of its subclasses (e.g. AbsLine).

Sub Classes
===========

The primary children of SpectralLine are
:ref:`AbsLine` and EmissionLine (to be implemented).
See their documentation for a description of instantiation
and additional attributes.

Attributes
==========

The base attributes for the SpectralLine class are:

================ ================= ========== ========================================
Property         Variable          Type       Description
================ ================= ========== ========================================
Limits           limits            LineLimits The redshift and limits of the line in redshift, velocity
                                              (w/r to its redshift) and observed wavelength.
RA, Dec          attrib['coord']   Coord      astropy.coordinate
Velocity         attrib['v']       Quantity   line velocity relative to its redshift
Velocity sigma   attrib['sig_v']   Quantity   1 sigma uncertainty in the velocity
Equivalent Width attrib['EW']      Quantity   Equivalent width
EW sigma         attrib['sig_EW']  Quantity   1 sigma uncertainty in EW
EW flag          attrib['flag_EW'] int        Equivalent width flag
================ ================= ========== ========================================

Note that redshift is sufficiently important that it is contained
within its own object.  It is also accessible as a property, e.g.::

   z = sline.z

.. _spec-analysis:

Analysis
========

It is common that one wishes to associate a line with a spectrum
to perform a range of analyses.
This is accomplished through::

   spline.analy['spec'] = sp

where sp is an :ref:`XSpectrum1D` object.

Methods
=======

cut_spec
--------

Provide a spectrum has been associated to the line (see `Analysis`_):
then this method returns the portion of the spectrum surrounding
the line.  The limits are specified in the LineLimits class held
in the attribute *limits*,
usually either with observed wavelengths or velocities relative
to the line's redshift.  The code returns the flux, error array,
and a *dict* containing the wavelength and velocity arrays.
::

   spline.limits.set([-300., 300.]*u.km/u.s) # vlim
   fx, sig, wv_dict = spline.cut_spec()

ismatch
-------

Check whether two Lines are equal rest wavelength (to 1e-6 AA),
whether they have common RA/DEC to 0.1" (if provided),
and whether they are from the same ion::

   print specline.ismatch(another_line)


.. _measure-ew:

measure_ew
----------

Measure the observed Equivalent width for a SpectralLine.
Following absorption-line convention, absorption will
have a positive value and emission will have a negative value.

To perform the calculation, the line must be associated to
a spectrum (see `Analysis`_ above) and the LineLimits of the line
must have previously been specified.

When executed, the EW and sig_EW attibutes are filled::

   specline.measure_ew()
   print(specline.attrib['EW'])


measure_kin
-----------

Measure kinematic characteristics of an AbsLine.
To perform the calculation, the line must be associated to
a spectrum (see `Analysis`_) and vlim must
be specified.  When executed, the 'kin' attribute is filled
with a dict of measurements.  Default set of measurements
are the v90, fedg, and fmm statistics of Prochaska & Wolfe 1997::

   specline.measure_kin()


measure_restew
--------------

Measure the rest-frame Equivalent width of a SpectralLine.
See :ref:`measure-ew` for other details.

to_dict
-------

Convert the Class to a JSON-friendly dict that might
be easily written to the disk, e.g.::

   adict = specline.to_dict()
   with io.open(outfil, 'w', encoding='utf-8') as f:
      f.write(unicode(json.dumps(tdict, sort_keys=True,
         indent=4, separators=(',', ': '))))



Utilities
=========

There are several utilites related to spectral lines.
These are located in the line_utils module.

parse_speclines
---------------

Given a list of SpectralLines and desired property (key),
this method returns a list or array of the values.::

   from linetools import line_utils
   array_of_values = line_utils.parse_speclines(list_of_speclines, mk_array=True)

transtable_from_speclines
-------------------------

Given a list of SpectralLines, this method returns a Table
of a subset of the properties (e.g. wavelength, name, EW).::

   trans_tbl = line_utils.transtable_from_speclines(list_of_speclines)


