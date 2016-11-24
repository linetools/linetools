.. _AbsSightline:

******************
AbsSightline Class
******************

.. index:: AbsSightline

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <AbsSightline>

Overview
========

This Class is designed to organize the absorption systems
along a single sightline.  It may be most commonly used
for extragalactic sightlines, but it can be applied to the
ISM as well.

By definition, an AbsSightline is a unique collection of
absorption components.  The only quantities required
 to define the AbsSightline are its coordinates on the sky.


Instantiation
=============

The AbsSightline Class may be instantiated in a few ways.
The default sets the properties listed above::

	abssl = GenericAbsSightline((10.0*u.deg, 45*u.deg))

More commonly, one will instantiate with one
a set of components::

   lya = AbsLine('HI 1215', z=2.3)
    lya.limits.set([-300.,300.]*u.km/u.s)  # vlim
    lyb = AbsLine(1025.7222*u.AA, z=2.3)
    lyb.limits.set([-300.,300.]*u.km/u.s)  # vlim
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.coord = ltu.radec_to_coord((10.*u.deg, 45*u.deg))
    abssl = GenericAbsSightline.from_components([abscomp])



Inspecting
==========

Here are a few simple methods to explore/inspect the class.

Generate a Table
++++++++++++++++

If the class contains one or more AbsComponent obejcts, you may generate a
`~astropy.table.Table` from their attributes and data::

    comp_tbl = abssl.build_table()


I/O
===

One may generate a *dict* of the key properties of the AbsSystem
with the to_dict() method::

   asldict = abssl.to_dict()


This can then be written to disk with a JSON or yaml dump.
