.. _GUIs:

****
GUIs
****

.. index:: GUIs

Notebooks
=========

.. toctree::
   :maxdepth: 1

   XSpecGui <xspecgui>
   XAbsSysGui <xabssysgui>

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
spectrum. A series of 'knot' positions are estimated, and these are
then joined with a spline to produce a continuum. Spline points can be
interactively added, deleted or moved to improve the continuum. See
the notebook for an example or use the lt_continuumfit script directly.

XSpecGUI
========

This enables visual inspection of a spectrum.  Simple analysis
(e.g. equivalent width measurements) may also be performed.  See the
notebook for details.


XAbsSysGui
==========

This shows a velocity (stack) plot of the absorption lines from
an input absorption line system.  The user can then modify the
velocity limits that would be used for subsequent analysis, flag
bad lines, blends, set limits, etc.  The modified absoprtion system
is then written to the hard-drive as a JSON file.
