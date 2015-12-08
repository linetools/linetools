.. highlight:: rest

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

This Class is designed to organize and analyze an spectral line,
either emission or absorption. This is an abstract base class,
i.e. one must instantiate one of its subclasses (e.g. AbsLine).

By definition, a AbsLine is specified by:

================ ================= ============== ========================================
Property         Variable          Type           Description
================ ================= ============== ========================================
RA, Dec          attrib['coord']   coord          astropy.coordinate
Redshift         attrib['z']       float          Reference redshift
Velocity         attrib['v']       Quantity       line velocity relative to its redshift
Velocity sigma   attrib['sig_v']   Quantity       1 sigma uncertainty in the velocity
Equivalent Width attrib['EW']      Quantity       Equivalent width
EW sigma         attrib['sig_EW']  Quantity       1 sigma uncertainty in EW
EW flag          attrib['flag_EW'] int            Equivalent width flag
Doppler param.   attrib['b']       Quantity       doppler parameter
b sigma          attrib['sig_b']   Quantity       1 sigma uncertainty in b
Column density   attrib['N']       Quantity       Column density
N sigma          attrib['sig_N']   Quantity       1 sigma uncertainty in N
N flag           attrib['flag_N']  int            Column density flag
================ ================= ============== ========================================


Instantiation
=============

Attributes
==========

Sub Classes
===========

Generic
-------

Plots
=====

Methods
=======

Output
======
