.. highlight:: rest

******************
AbsComponent Class
******************

.. index:: AbsComponent

Overview
========

This Class is designed to organize and analyze a set of
absorption lines (:ref:`AbsLine Class`).

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
-------------

The AbsComponent Class may be instantiated in a few ways.
The default sets the properties listed above::

	abscomp = AbsComponent((10.0*u.deg, 45*u.deg), (14,2), 1.0, [-300,300]*u.km/u.s)

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

Inspecting
----------

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
