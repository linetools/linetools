.. _AbsSystem:

******************
AbsSystem Class
******************

.. index:: AbsSystem

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <AbsSystem_examples>

Overview
========

This Class is designed to organize and analyze an absorption system.
This is generally constructed of one or more :ref:`AbsComponent`.
The base class is abstract, i.e. one must instantiate one of its
flavors (e.g. HILyman, MgII, LLS, DLA).

By definition, an AbsSystem is a unique collection of
absorption components.  It is specified by:

=============== ========   ============== ============================================
Property        Variable   Type           Description
=============== ========   ============== ============================================
RA, DEC         radec      tuple or coord RA,DEC in deg or astropy.coordinate
Redshift        z          float          absorption redshift
Velocity limits vlim       Quantity array -/+ velocity limits of the system
=============== ========   ============== ============================================


Instantiation
=============

The AbsSystem Class may be instantiated in a few ways.
The default sets the properties listed above::

	gensys = GenericAbsSystem((15.23*u.deg,-23.22*u.deg), 1.244, [-500,500]*u.km/u.s, NHI=16.)

More commonly, one will instantiate with one or more AbsComponent objects::

    # HI Lya, Lyb
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    lya = AbsLine(1215.670*u.AA, z=2.92939)
    lya.limits.set([-300.,300.]*u.km/u.s)  # vlim
    lyb = AbsLine(1025.7222*u.AA, z=lya.attrib['z'])
    lyb.limits.set([-300.,300.]*u.km/u.s)  # vlim
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.coord = radec
    # Finish
    HIsys = LymanAbsSystem.from_components([abscomp])


One may also instantiate from a *dict*, usually read
from the hard-drive::

    abscomp = AbsSystem.from_dict(idict)

Attributes
==========

Sub Classes
===========

Generic
-------

A catch-all subclass for AbsSystem.
More options are provided in
`pyigm <https://github.com/pyigm/pyigm>`_.


Plots
=====

Methods
=======

AbsLines
--------

There are a few methods related to the AbsLine objects within
an AbsSystem.

List
++++

One can generate a list of all the AbsLine objects
with::

   lines = abssys.list_of_abslines()


get absline
+++++++++++

One can retrieve one or more AbsLine objects matching the name
or rest-wavelength of a transition, e.g. ::

   lyb = abssys.get_absline('HI 1025')
   # or
   lyb = abssys.get_absline(1025.72*u.AA)  # Nearest 0.01 A is required

measure_restew
++++++++++++++

Measure the rest-frame equivalent widths for all AbsLine objects
in the system.
In addition to the inputs described in :ref:`measure-ew`,
each line to be measured must have line.analy['do_anlaysis']=1.
Here is an example::

   abssys.measure_restew(spec=xspec_object)


fill transitions
++++++++++++++++

Generate an astropy.Table of information on the absorption lines
in the system.  This is stored in the ._trans attribute::

   abssys.fill_trans()
   print(abssys._trans)

Components
----------

`~linetools.igsm.abssystem.get_component`
grabs components matching an input where the input is either
a tuple of (Z, ion) or an AbsLine::

   SiII = gensys.get_component((14,2))

`~linetools.igsm.abssystem.update_component_colm` synthesizes
and updates the column densities for the components.::

   gensys.update_component_colm()


ionN
----

Fill the _ionN attribute with a QTable of column densities.
These are derived from the components ::

   abssys.fill_ionN()
   print(abssys._ionN)



Output
======

One may generate a *dict* of the key properties of the AbsSystem
with the to_dict() method::

   odict = HIsys.to_dict()

This dict is required to be JSON compatible.


