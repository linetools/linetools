.. _SolarAbund:

******************
SolarAbund Class
******************

.. index:: SolarAbund

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <SolarAbund>

Overview
========

This class provides access to element abundance ratios measured in (or
nearby) the Sun.

To access the abundances, make a new class instance::

  >>> from linetools.abund import solar as labsol
  >>> sol = labsol.SolarAbund()
  Loading abundances from Asplund2009
  Abundances are relative by number on a logarithmic scale with H=12

Then select the element you want by either its name or atomic number::

  >>> print(sol['C'])
  8.43
  >>> print(sol[6])
  8.43

Currently the abundances from Asplund et al. 2009 are available, and
in future more references will be included.

Multiple elements can also be selected::

  >>> print(sol['C', 'O']) # doctest: +SKIP
  [8.43 8.69]

Element ratios can be accessed using the ``get_ratio`` method::

  >>> print(sol.get_ratio('C/Fe')) # doctest: +SKIP
  0.98
