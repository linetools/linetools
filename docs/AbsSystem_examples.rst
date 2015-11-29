
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


.. parsed-literal::

    WARNING: UnitsWarning: The unit 'Angstrom' has been deprecated in the FITS standard. Suggested: 10**-1 nm. [astropy.units.format.utils]
    WARNING:astropy:UnitsWarning: The unit 'Angstrom' has been deprecated in the FITS standard. Suggested: 10**-1 nm.


.. parsed-literal::

    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/morton03_table2.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/morton00_table2.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/verner94_tab6.fits

.. parsed-literal::

    WARNING: UnitsWarning: '0.1nm' did not parse as fits unit: Numeric factor not supported by FITS [astropy.units.core]
    WARNING:astropy:UnitsWarning: '0.1nm' did not parse as fits unit: Numeric factor not supported by FITS


.. parsed-literal::

    
    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/EUV_lines.ascii
    read_sets: Using set file -- 
      /Users/xavier/local/Python/linetools/linetools/lists/sets/llist_v0.3.ascii


.. code:: python

    # HILyman system
    reload(lt_absys)
    HIsys = lt_absys.LymanAbsSystem.from_components([abscomp])
    print(HIsys)
    print(HIsys._components)


.. parsed-literal::

    [LymanAbsSystem: name= type=HILyman, 08:12:27.432 -12:25:55.56, z=2.92939, NHI=0]
    [[AbsComponent: 08:12:27.432 -12:25:55.56, Zion=(1,1), z=2.92939, vlim=-300 km / s,300 km / s]]


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
    [[AbsComponent: 08:12:27.432 -12:25:55.56, Zion=(1,1), z=2.92939, vlim=-300 km / s,300 km / s], [AbsComponent: 08:12:27.432 -12:25:55.56, Zion=(14,2), z=2.92939, vlim=-250 km / s,80 km / s]]


.. code:: python

    lya.data




.. parsed-literal::

    {'A': <Quantity 626500000.0 1 / s>,
     'Am': 0,
     'Ej': <Quantity 0.0 1 / cm>,
     'Ek': <Quantity 2259.163 1 / cm>,
     'Ex': <Quantity 0.0 1 / cm>,
     'Jj': 0.0,
     'Jk': 0.0,
     'Ref': 'Morton2003',
     'Z': 1,
     'col0': masked,
     'col6': masked,
     'el': 0,
     'f': 0.41639999999999999,
     'gamma': <Quantity 626500000.0 1 / s>,
     'gj': 2,
     'gk': 6,
     'group': 1,
     'ion': 1,
     'mol': '',
     'name': 'HI 1215',
     'nj': 0,
     'nk': 0,
     'wrest': <Quantity 1215.67 Angstrom>}



