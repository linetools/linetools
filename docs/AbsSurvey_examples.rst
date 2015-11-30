
Examples for AbsSurvey (v1.1)
=============================

.. code:: python

    # imports
    from astropy.coordinates import SkyCoord
    
    from linetools.isgm.abssystem import GenericAbsSystem
    from linetools.isgm.abssurvey import GenericAbsSurvey

Simple instantiation
--------------------

.. code:: python

    gensurvey = GenericAbsSurvey()
    gensurvey




.. parsed-literal::

    [AbslineSurvey: nsys=0, type=Generic, ref=]



.. code:: python

    coord = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gensys = GenericAbsSystem(coord, 1.244, [-300,300.]*u.km/u.s, NHI=16.)
    gensys.name = 'Sys1'
    #
    coord2 = SkyCoord(ra=223.1143*u.deg, dec=42.4321*u.deg)
    gensys2 = GenericAbsSystem(coord2, 1.744, [-300,300.]*u.km/u.s, NHI=17.)
    gensys2.name = 'Sys2'

.. code:: python

    gensurvey.nsys = 2
    gensurvey._abs_sys.append(gensys)
    gensurvey._abs_sys.append(gensys2)

.. code:: python

    gensurvey._abs_sys




.. parsed-literal::

    [[GenericAbsSystem: name=Sys1 type=Generic, 08:12:27.432 -12:25:55.56, z=1.244, NHI=16],
     [GenericAbsSystem: name=Sys2 type=Generic, 14:52:27.432 +42:25:55.56, z=1.744, NHI=17]]



Parsing
-------

.. code:: python

    gensurvey.NHI




.. parsed-literal::

    array([ 16.,  17.])



.. code:: python

    gensurvey.zabs




.. parsed-literal::

    array([ 1.244,  1.744])



