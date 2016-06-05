.. _RelAbund:

**************
RelAbund Class
**************

.. index:: RelAbund

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <RelAbund>

Overview
========

This class packages the relative abundances of an object, typically
an AbsSystem.

Instantation
============

To begin, make a new class instance::

   XY = RelAbund()
   Loading abundances from Asplund2009

Then load data into the data dict.  Here is an example::

   XY._data = {6: dict(flag=1, XH=-1., sigXH=0.2, sig=0.05),
            8: dict(flag=2, XH=-1.4, sigXH=0.25, sig=0.05),
            14: dict(flag=1, XH=-1.1, sigXH=0.25, sig=0.05),
            26: dict(flag=1, XH=-1.4, sigXH=0.25, sig=0.05),
            32: dict(flag=3, XH=-0.8, sigXH=0.25, sig=0.05),
            }

More frequent usage will be to instantate using an input
table of column density measurements, e.g.::

   dla.XY = RelAbund.from_ionclm_table((1,dla.NHI, dla.sig_NHI[0]), dla._ionN)

See pyigm DLA abund Notebook for more.

Usage
=====

You may grab the data for any element with item syntax::

   CH = XY[6]
   {'flag': 1, 'sig': 0.2, 'val': -1.0}
   CH = XY['C']

Element ratios can be accessed by providing a tuple of
atomic number or element name::

   SiFe = XY['Si', 'Fe']
   {'flag': 1, 'sig': 0.070710678118654766, 'val': 0.2999999999999998}

You may generate an astropy Table of the X/Y values::

   tbl = XY.table()  # For X/H
   tbl = XY.table('Fe')  # For X/Fe

