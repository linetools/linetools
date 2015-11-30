.. highlight:: rest

******************
AbsSurvey Class
******************

.. index:: AbsSystem

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <AbsSurvey_examples>

Overview
========

This Class is designed to organize and analyze a survey of
absorption systems (defined as AbsSystem objects).
The base class is abstract, i.e. one must instantiate
one of its flavors (e.g. Generic).

By definition, an AbsSurvey is a unique collection of
AbsSystem objects.  It is specified by the number of
systems and the references.


Instantiation
-------------

The AbsSystem Class may be instantiated in a few ways.
The default sets the properties listed above::

	gensurvey = GenericAbsSurvey()

More commonly, one will instantiate with one or more AbsSystem objects::

    coord = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gensys = GenericAbsSystem(coord, 1.244, [-300,300.]*u.km/u.s, NHI=16.)
    gensys.name = 'Sys1'
    #
    coord2 = SkyCoord(ra=223.1143*u.deg, dec=42.4321*u.deg)
    gensys2 = GenericAbsSystem(coord2, 1.744, [-300,300.]*u.km/u.s, NHI=17.)
    gensys2.name = 'Sys2'

::::

Attributes
----------

Plots
-----

Methods
-------

Output
------
