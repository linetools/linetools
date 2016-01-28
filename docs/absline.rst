.. _AbsLine:

*************
AbsLine Class
*************

.. index:: AbsLine

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <AbsLine_examples>


Overview
========

This Class is a child of the abstract
:ref:`SpectralLine` Class.  See that
documentation for the base methods.

AbsLine is designed to organize and analyze an absorption line.
In addition to the attributes defaulted to SpectralLine,
this class has:

================ ================= ========= ========================================
Property         Variable          Type      Description
================ ================= ========= ========================================
Doppler param.   attrib['b']       Quantity  Doppler parameter
b sigma          attrib['sig_b']   Quantity  1 sigma uncertainty in b
Column density   attrib['N']       Quantity  Column density
N sigma          attrib['sig_N']   Quantity  1 sigma uncertainty in N
N flag           attrib['flag_N']  int       Column density flag
================ ================= ========= ========================================


Instantiation
=============

from_dict
---------

Instantiate from a dict.  The keys *ltype* ('Abs')
and *trans* are required.

fill_data
---------

Attributes
==========

See the Table above.  logN and sig_logN are commonly used
as well.

Plots
=====

Methods
=======

generate_voigt
--------------

measure_aodm
------------

Output
======
