.. _AbsComponent:

******************
AbsComponent Class
******************

.. index:: AbsComponent

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <AbsComponent_examples>
   Column Densities <AbsComponent_ColumnDensities>

Overview
========

This Class is designed to organize and analyze a set of
absorption lines.

By definition, an AbsComponent is a unique collection of
absorption lines specified by:

=============== ========   ============== ============================================
Property        Variable   Type           Description
=============== ========   ============== ============================================
RA, DEC         radec      tuple or coord RA,DEC in deg or astropy.coordinate
Z, ion          Zion       tuple          Atomic Number (Z) and ionization state (ion)
Redshift        z          float          absorption redshift
Velocity limits vlim       Quantity array -/+ velocity limits of the component
Energy level    Ej         Quantity       Energy of the level relative to ground
Isotope         A          int            Atomic Mass number (optional)
=============== ========   ============== ============================================


Instantiation
=============

The AbsComponent Class may be instantiated in a few ways.
The default sets the properties listed above::

	abscomp = AbsComponent((10.0*u.deg, 45*u.deg), (14,2), 1.0, [-300,300]*u.km/u.s)

More commonly, one will instantiate with one
`~linetools.spectralline.AbsLine` object::

    lya = AbsLine(1215.670*u.AA, z=2.92939)
    lya.limits.set([-300.,300.]*u.km/u.s)  # vlim
    abscomp1 = AbsComponent.from_abslines([lya])

or multiple::

    lyb = AbsLine(1025.7222*u.AA, z=lya.attrib['z'])
    lyb.limits.set([-300.,300.]*u.km/u.s)  # vlim
    abscomp = AbsComponent.from_abslines([lya,lyb])

One may also instantiate from a *dict*, usually read
from the hard-drive::

    abscomp = AbsComponent.from_dict(idict)

::::

Inspecting
==========

Here are a few simple methods to explore/inspect the class.

Generate a QTable
+++++++++++++++++

If the class contains one or more AbsLines, you may generate a
`~astropy.table.QTable` from their attributes and data::

    comp_tbl = abscomp.build_table()

Show a Stack Plot
+++++++++++++++++

If the AbsLine have spectra attached to them (in attrib['spec']),
a stack plot (aka velocity plot) is generated with::

    abscomp.stack_plot()

Apparent Column Densitities
+++++++++++++++++++++++++++

Show a plot of the apparent column density profiles, :math:`N_a`::

    abscomp.plot_Na()

::::

Analysis
========

Here are some methods related to analysis.

Synthesize Columns
++++++++++++++++++

If one inputs a set of AbsLine(s) with column density measurements,
the synthesize_colm method collates these.  Positive, unsaturated detections
are combined in a weighted mean whereas limits are compared
and the most sensitive one is adopted.::

    abscomp.synthesize_colm()

Here is the set of rules:

1.  If all measurements are upper limits, take the lowest value and flag as an upper limit (*flgN=3*).
2.  If all measurements are a mix of upper and lower limits, take the highest lower limit and flag as a lower limit (*flgN=2*).
3.  If one or more measurements are a proper detection, take the weighted mean of these and flag as a detection (*flgN=1*).

Curve of Growth
+++++++++++++++

A standard, single-component curve-of-growth (COG) analysis may be
performed on the set of AbsLines::

    COG_dict = abscomp.cog(show_plot=True)

The output dict includes:

========== ============== =====================================
Key        Type           Description
========== ============== =====================================
EW         Quantity array Input equivalent widths
sigEW      Quantity array Input error in equivalent widths
f          ndarray        Input f-values
wrest      Quantity array Input rest wavelengths
logN       float          Output fitted column density (log10)
sig_logN   float          Output error in fitted logN
b          Quantity       Output b-value (km/s)
sig_b      Quantity       Output error in b-value (km/s)
========== ============== =====================================

Misc
====

I/O
+++

One may generate a *dict* of the key properties of the AbsSystem
with the to_dict() method::

   cdict = component.to_dict()


Synthesize Components
+++++++++++++++++++++

This method combines a list of two or more components into a new one.
It checks first for consistent RA/DEC, Zion, and Ej.  It does
not place any constraints on z and vlim.  The column density of
the new component is the sum of the input ones (with rules for
limits).  And the redshift and vlim are set to encompass the
velocities of the input components.::

   from linetools.isgm import utils as ltiu
   synth_SiII = ltiu.synthesize_components([SiIIcomp1,SiIIcomp2])

See the :doc:`AbsComponent_examples` notebook for a complete example.

Generate Multiple Components
++++++++++++++++++++++++++++

This method generates multiple components from a list of
AbsLines.::

   comps = ltiu.build_components_from_abslines([lya,lyb,SiIIlines[0],SiIIlines[1]])

