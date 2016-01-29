.. _GUIs:

****
GUIs
****

.. index:: GUIs

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Continuum fitting <Interactive_continuum_fitting_example>
   XSpecGui <xspecgui>

Overview
========

There are several GUIs included with linetools, primarily for
simple spectral inspection and analysis.

We caution that it is difficult (essentially impossible) to
generate unit tests for these GUIs.  As such, they are far
from bug free and may crash unexpectedly.  Buyer beware!

Continuum fitting
=================

This enables interactive fitting of the unabsorbed continuum for a
spectrum. A series of 'knots' positions are estimated, and these are
are joined with a spline to produce a continuum. Spline points can
then be interactively added, deleted or moved to improve the
continuum. See the notebook for an example.

XSpecGUI
========

This enables visual inspection of a spectrum.  Simple analysis
(e.g. EW measurements) may also be performed.  See the notebook for
details.



