.. highlight:: rest

.. _LineList:

******************
LineList Class
******************

.. index:: LineList

Overview
========

This Class is designed to organize and handle information about atomic
and/or molecular transition lines (e.g. HI Lya, CIV 1548, Hydrogen
Balmer series, etc.) observed in a variety of astrophysical
environments. It is currently implemented for absorption lines, but we
expect to also include common emission lines in the near future.

..
   (:ref:`AbsLine Class`).  add this back in when written

By definition, a LineList is a unique collection of transitions lines
specified by:

=============== =========   ============== ====================================================
Property        Variable    Type           Description
=============== =========   ============== ====================================================
LineList key    llst_keys   str or list    A key to define a subsample of transitions to load:
                                            * 'ISM'   :: "All" ISM lines (can be overwhelming!)
                                            * 'Strong':: Strong ISM lines
                                            * 'HI'    :: HI Lyman series
                                            * 'H2'    :: H2 (Lyman-Werner)
                                            * 'CO'    :: CO UV band-heads
                                            * 'EUV'   :: Extreme UV lines
=============== =========   ============== =====================================================


Instantiation
=============

The LineList Class may be instantiated in a few ways::

	ism = LineList(['ISM'])
  linetools.lists.parse: Reading linelist --- 
      /home/ntejos/python/linetools/linetools/data/lines/morton03_table2.fits.gz
  WARNING: UnitsWarning: The unit 'Angstrom' has been deprecated in the FITS standard. Suggested: nm (with data multiplied by 0.1). [astropy.units.format.utils]
  linetools.lists.parse: Reading linelist --- 
      /home/ntejos/python/linetools/linetools/data/lines/morton00_table2.fits.gz
  linetools.lists.parse: Reading linelist --- 
      /home/ntejos/python/linetools/linetools/data/lines/verner94_tab6.fits
  linetools.lists.parse: Reading linelist --- 
      /home/ntejos/python/linetools/linetools/data/lines/EUV_lines.ascii
  read_sets: Using set file -- 
      /home/ntejos/python/linetools/linetools/lists/sets/llist_v0.4.ascii



More commonly, one will instantiate with one AbsLine object::

    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    abscomp1 = AbsComponent.from_abslines([lya])

or multiple::

    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    abscomp = AbsComponent.from_abslines([lya,lyb])

::::

Inspecting
==========

Here are a few simple methods to explore/inspect the class.

Generate a QTable
+++++++++++++++++

If the class contains one or more AbsLines, you may generate a QTable
from their attributes and data::

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

