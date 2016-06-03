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

================ ================= ========= ========================================
Property         Variable          Type      Description
================ ================= ========= ========================================
RA, Dec          attrib['coord']   Coord     astropy.coordinate
Redshift         attrib['z']       float     Reference redshift
Redshift sigma   attrib['sig_z']   float     Reference redshift uncertainty
Velocity         attrib['v']       Quantity  line velocity relative to its redshift
Velocity sigma   attrib['sig_v']   Quantity  1 sigma uncertainty in the velocity
Equivalent Width attrib['EW']      Quantity  Equivalent width
EW sigma         attrib['sig_EW']  Quantity  1 sigma uncertainty in EW
EW flag          attrib['flag_EW'] int       Equivalent width flag
================ ================= ========= ========================================

.. _specanalysis

You can access these attributes with a getitem call, e.g.::

   specline['EW']

A search is performed through the object attributes (e.g. wrest),
the attrib dict, the analy dict, and then the data dict.

Analysis
========

It is common that one wishes to associate a line with a spectrum
to perform a range of analyses.
This is accomplished through::

   specline.analy['spec'] = spec

where spec is an :ref:`XSpectrum1D` object.

Methods
=======

cut_spec
--------

Provide a spectrum has been associated to the line (see `Analysis`_):
then this method returns the portion of the spectrum surrounding
the line.  The limits are specified by either analy['wvlim'] (in
observed wavelength) or analy['vlim'] with velocities relative
to the line's redshift.  The code returns the flux, error array,
and a *dict* containing the wavelength and velocity arrays.
::

   specline.analy['vlim'] = [-300., 300.]*u.km/u.s
   fx, sig, wv_dict = spline.cut_spec()

ismatch
-------

Check whether two Lines are equal rest wavelength (to 1e-6 AA),
whether they have common RA/DEC to 0.1" (if provided),
and whether they are from the same ion::

   print specline.ismatch(another_line)


measure_ew
----------

Measure the observed Equivalent width for a SpectralLine.
Following absorption-line convention, absorption will
have a positive value and emission will have a negative value.

To perform the calculation, the line must be associated to
a spectrum (see `Analysis_`) and either wvlim or vlim must
be specified.  When executed, the EW and sig_EW attibutes
are filled::

   specline.measure_ew()


measure_kin
-----------

Measure kinematic characteristics of an AbsLine.
To perform the calculation, the line must be associated to
a spectrum (see `Analysis_`) and vlim must
be specified.  When executed, the 'kin' attribute is filled
with a dict of measurements.  Default set of measurements
are the v90, fedg, and fmm statistics of Prochaska & Wolfe 1997::

   specline.measure_kin()


measure_restew
--------------

Measure the rest-frame Equivalent width of a SpectralLine.
See `measure_ew`_ for details.

to_dict
-------

Convert the Class to a JSON-friendly dict that might
be easily written to the disk, e.g.::

   adict = specline.to_dict()
   with io.open(outfil, 'w', encoding='utf-8') as f:
      f.write(json.dumps(tdict, sort_keys=True,
         indent=4, separators=(',', ': ')))


