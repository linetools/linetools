
Examples for AbsSystem Class (v1.0)
===================================

.. code:: python

    # imports
    import imp
    from astropy.coordinates import SkyCoord
    
    from linetools.isgm import abssystem as lt_absys
    from linetools.spectralline import AbsLine
    from linetools.isgm.abscomponent import AbsComponent

Simple instantiation
--------------------

Standard init
~~~~~~~~~~~~~

.. code:: python

    reload(lt_absys)
    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gensys = lt_absys.GenericAbsSystem(radec, 1.244, [-500,500]*u.km/u.s, NHI=16.)
    gensys




.. parsed-literal::

    [GenericAbsSystem: name=Foo type=Generic, 08:12:27.432 -12:25:55.56, z=1.244, NHI=16]



From components
~~~~~~~~~~~~~~~

One component
^^^^^^^^^^^^^

.. code:: python

    # HI Lya, Lyb
    lya = AbsLine(1215.670*u.AA)
    lya.analy['vlim'] = [-300.,300.]*u.km/u.s
    lya.attrib['z'] = 2.92939
    lyb = AbsLine(1025.7222*u.AA)
    lyb.analy['vlim'] = [-300.,300.]*u.km/u.s
    lyb.attrib['z'] = lya.attrib['z']
    abscomp = AbsComponent.from_abslines([lya,lyb])
    abscomp.coord = radec

.. code:: python

    # HILyman system
    reload(lt_absys)
    HIsys = lt_absys.LymanAbsSystem.from_components([abscomp])
    print(HIsys)
    print(HIsys._components)


.. parsed-literal::

    [LymanAbsSystem: name= type=HILyman, 08:12:27.432 -12:25:55.56, z=2.92939, NHI=0]
    [[AbsComponent: 08:12:27.432 -12:25:55.56, Zion=(1,1), z=2.92939]]


Multiple
^^^^^^^^

.. code:: python

    # SiII
    SiIItrans = ['SiII 1260', 'SiII 1304', 'SiII 1526', 'SiII 1808']
    abslines = []
    for trans in SiIItrans:
        iline = AbsLine(trans)
        iline.attrib['z'] = 2.92939
        iline.analy['vlim'] = [-250.,80.]*u.km/u.s
        abslines.append(iline)
    #
    SiII_comp = AbsComponent.from_abslines(abslines)
    SiII_comp.coord = radec

.. code:: python

    # LLS (coming)
    reload(lt_absys)
    LLSsys = lt_absys.GenericAbsSystem.from_components([abscomp,SiII_comp])
    print(LLSsys)
    print(LLSsys._components)


.. parsed-literal::

    [GenericAbsSystem: name=Foo type=Generic, 08:12:27.432 -12:25:55.56, z=2.92939, NHI=0]
    [[AbsComponent: 08:12:27.432 -12:25:55.56, Zion=(1,1), z=2.92939], [AbsComponent: 08:12:27.432 -12:25:55.56, Zion=(14,2), z=2.92939]]


