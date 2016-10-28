.. _EmLine:

************
EmLine Class
************

.. index:: EmLine

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <EmLine_examples>


Overview
========

This Class is a child of the abstract
:ref:`SpectralLine` Class.  See that
documentation for the base methods.

EmLine is designed to organize and analyze an emission line.
In addition to the attributes defaulted to SpectralLine,
this class has:

================ =================== ========= ========================================
Property         Variable            Type      Description
================ =================== ========= ========================================
Flux             attrib['flux']      Quantity  Line flux (erg/s)
Flux sigma       attrib['sig_flux']  Quantity  1 sigma uncertainty in flux
Flux flag        attrib['flag_flux'] int       Flux flag
================ =================== ========= ========================================


Instantiation
=============

The typical way to instantiate is a standard call with the
rest wavelength or name of the transition::

   emisslin = EmLine('Halpha')
   emisslin = EmLine(6564.613*u.AA)

By default the class searches the Galaxy LineList.

from_dict
---------

Instantiate from a dict.  The keys *ltype* ('Em')
and *trans* are required.

fill_data
---------

Attributes
==========

See the Table above.

Plots
=====

Methods
=======


Output
======
