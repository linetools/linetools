
Examples for AbsSystem Class (v1.1)
===================================

:download:`Download <examples/AbsSystem_examples.ipynb>` this notebook.

.. code:: python

    # suppress warnings for these examples
    import warnings
    warnings.filterwarnings('ignore')
    
    # imports
    import imp
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    
    from linetools.isgm import abssystem as lt_absys
    from linetools.spectralline import AbsLine
    from linetools.isgm.abscomponent import AbsComponent

Simple instantiation
--------------------

Standard init
~~~~~~~~~~~~~

.. code:: python

    radec = SkyCoord(ra=123.1143*u.deg, dec=-12.4321*u.deg)
    gensys = lt_absys.GenericAbsSystem(radec, 1.244, [-500,500]*u.km/u.s, NHI=16.)
    gensys




.. parsed-literal::

    <GenericAbsSystem: name=Foo type=Generic, 08:12:27.432 -12:25:55.56, z=1.244, NHI=16>



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


.. parsed-literal::

    linetools.lists.parse: Reading linelist --- 
       /Users/ncrighton/Code/Repo/linetools/build/lib.macosx-10.5-x86_64-3.4/linetools/data/lines/morton03_table2.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/ncrighton/Code/Repo/linetools/build/lib.macosx-10.5-x86_64-3.4/linetools/data/lines/morton00_table2.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/ncrighton/Code/Repo/linetools/build/lib.macosx-10.5-x86_64-3.4/linetools/data/lines/verner96_tab1.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/ncrighton/Code/Repo/linetools/build/lib.macosx-10.5-x86_64-3.4/linetools/data/lines/verner94_tab6.fits
    linetools.lists.parse: Reading linelist --- 
       /Users/ncrighton/Code/Repo/linetools/build/lib.macosx-10.5-x86_64-3.4/linetools/data/lines/EUV_lines.ascii
    read_sets: Using set file -- 
      /Users/ncrighton/Code/Repo/linetools/build/lib.macosx-10.5-x86_64-3.4/linetools/lists/sets/llist_v1.0.ascii


.. code:: python

    # HILyman system
    HIsys = lt_absys.LymanAbsSystem.from_components([abscomp])
    print(HIsys)
    print(HIsys._components)


.. parsed-literal::

    <LymanAbsSystem: name=J081227.432-122555.56_z2.929 type=HILyman, 08:12:27.432 -12:25:55.56, z=2.92939, NHI=0>
    [<AbsComponent: 08:12:27.432 -12:25:55.56, Name=HI_z2.92939, Zion=(1,1), Ej=0 1 / cm, z=2.92939, vlim=-300 km / s,300 km / s>]


Multiple components
^^^^^^^^^^^^^^^^^^^

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

    # Generic 
    imp.reload(lt_absys)
    LLSsys = lt_absys.GenericAbsSystem.from_components([abscomp,SiII_comp])
    print(LLSsys)
    print(LLSsys._components)


.. parsed-literal::

    <GenericAbsSystem: name=Foo type=Generic, 08:12:27.432 -12:25:55.56, z=2.92939, NHI=0>
    [<AbsComponent: 08:12:27.432 -12:25:55.56, Name=HI_z2.92939, Zion=(1,1), Ej=0 1 / cm, z=2.92939, vlim=-300 km / s,300 km / s>, <AbsComponent: 08:12:27.432 -12:25:55.56, Name=SiII_z2.92939, Zion=(14,2), Ej=0 1 / cm, z=2.92939, vlim=-250 km / s,80 km / s>]


Methods
-------

List of AbsLines
~~~~~~~~~~~~~~~~

.. code:: python

    lines = LLSsys.list_of_abslines()
    lines




.. parsed-literal::

    [<AbsLine: HI 1215, wrest=1215.6700 Angstrom>,
     <AbsLine: HI 1025, wrest=1025.7222 Angstrom>,
     <AbsLine: SiII 1260, wrest=1260.4221 Angstrom>,
     <AbsLine: SiII 1304, wrest=1304.3702 Angstrom>,
     <AbsLine: SiII 1526, wrest=1526.7070 Angstrom>,
     <AbsLine: SiII 1808, wrest=1808.0129 Angstrom>]



Single Line
~~~~~~~~~~~

.. code:: python

    lyb = LLSsys.get_absline('HI 1025')
    lyb




.. parsed-literal::

    <AbsLine: HI 1025, wrest=1025.7222 Angstrom>



.. code:: python

    lyb = LLSsys.get_absline(1025.72*u.AA)
    lyb




.. parsed-literal::

    <AbsLine: HI 1025, wrest=1025.7222 Angstrom>



.. code:: python

    lyb.wrest




.. math::

    1025.7222 \; \mathrm{\mathring{A}}



