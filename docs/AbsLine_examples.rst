
Examples for AbsLine class (v1.2)
=================================

:download:`Download <examples/AbsLine_examples.ipynb>` this notebook.

.. code:: python

    # suppress warnings for these examples
    import warnings
    warnings.filterwarnings('ignore')
    
    # import
    import astropy.units as u
    from linetools.spectralline import AbsLine, SpectralLine
    from linetools import spectralline as ltsp
    from linetools.spectra.xspectrum1d import XSpectrum1D

Generate a line
---------------

.. code:: python

    abslin = AbsLine(1548.195*u.AA)
    abslin


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




.. parsed-literal::

    <AbsLine: CIV 1548, wrest=1548.1950 Angstrom>



Data
~~~~

.. code:: python

    abslin.data




.. parsed-literal::

    {'A': <Quantity 0.0 1 / s>,
     'Am': 0,
     'Ej': <Quantity 0.0 1 / cm>,
     'Ek': <Quantity 0.0 1 / cm>,
     'Ex': <Quantity 0.0 1 / cm>,
     'Jj': 0.0,
     'Jk': 0.0,
     'Ref': 'Verner1994',
     'Z': 6,
     'col0': masked,
     'col7': masked,
     'el': 0,
     'f': 0.18999999761581421,
     'gamma': <Quantity 0.0 1 / s>,
     'gj': 2,
     'gk': 4,
     'group': 1,
     'ion': 4,
     'mol': '',
     'name': 'CIV 1548',
     'nj': 0,
     'nk': 0,
     'wrest': <Quantity 1548.195 Angstrom>}



As dict
~~~~~~~

.. code:: python

    abslin = AbsLine(1548.195*u.AA)
    tmp = abslin.to_dict()
    tmp




.. parsed-literal::

    {'analy': {'datafile': '',
      'do_analysis': 1,
      'flag_kin': 0,
      'flg_eye': 0,
      'flg_limit': 0,
      'name': 'CIV 1548',
      'vlim': {'unit': 'km / s', 'value': [0.0, 0.0]},
      'wvlim': {'unit': 'Angstrom', 'value': [0.0, 0.0]}},
     'attrib': {'DEC': 0.0,
      'EW': {'unit': 'Angstrom', 'value': 0.0},
      'N': {'unit': '1 / cm2', 'value': 0.0},
      'RA': 0.0,
      'b': {'unit': 'km / s', 'value': 0.0},
      'flag_EW': 0,
      'flag_N': 0,
      'sig_EW': {'unit': 'Angstrom', 'value': 0.0},
      'sig_N': {'unit': '1 / cm2', 'value': 0.0},
      'sig_b': {'unit': 'km / s', 'value': 0.0},
      'sig_v': {'unit': 'km / s', 'value': 0.0},
      'sig_z': 0.0,
      'v': {'unit': 'km / s', 'value': 0.0},
      'z': 0.0},
     'data': {'A': {'unit': '1 / s', 'value': 0.0},
      'Am': 0,
      'Ej': {'unit': '1 / cm', 'value': 0.0},
      'Ek': {'unit': '1 / cm', 'value': 0.0},
      'Ex': {'unit': '1 / cm', 'value': 0.0},
      'Jj': 0.0,
      'Jk': 0.0,
      'Ref': 'Verner1994',
      'Z': 6,
      'el': 0,
      'f': 0.18999999761581421,
      'gamma': {'unit': '1 / s', 'value': 0.0},
      'gj': 2,
      'gk': 4,
      'group': 1,
      'ion': 4,
      'mol': '',
      'name': 'CIV 1548',
      'nj': 0,
      'nk': 0,
      'wrest': {'unit': 'Angstrom', 'value': 1548.195}},
     'ltype': 'Abs',
     'name': 'CIV 1548',
     'wrest': {'unit': 'Angstrom', 'value': 1548.195}}



From dict
~~~~~~~~~

.. code:: python

    tmp2 = SpectralLine.from_dict(tmp)
    tmp2




.. parsed-literal::

    <AbsLine: CIV 1548, wrest=1548.1950 Angstrom>



Measure an EW
-------------

.. code:: python

    # Set spectrum
    abslin.analy['spec'] = XSpectrum1D.from_file('../../linetools/spectra/tests/files/UM184_nF.fits')

.. code:: python

    # Set analysis range
    abslin.analy['wvlim'] = [6080.78, 6087.82]*u.AA

.. code:: python

    # Measure
    abslin.measure_ew() # Observer frame
    print('EW = {:g} with error {:g}'.format(abslin.attrib['EW'],abslin.attrib['sig_EW']))


.. parsed-literal::

    EW = 0.993502 Angstrom with error 0.0527114 Angstrom


Measure AODM
------------

.. code:: python

    abslin.analy['wvlim'] = [0.,0.]*u.AA # Zero out for test
    #
    abslin.analy['spec'] = lsio.readspec('../../linetools/spectra/tests/files/UM184_nF.fits')
    abslin.analy['vlim'] = (-150., 150.)*u.km/u.s
    abslin.attrib['z'] = 2.92929

.. code:: python

    abslin.measure_aodm()
    N, sigN, flgN = [abslin.attrib[key] for key in ['N','sig_N','flag_N']] 
    print('logN = {:g}, siglogN = {:g}'.format(abslin.attrib['logN'], abslin.attrib['sig_logN']))


.. parsed-literal::

    logN = 13.9051, siglogN = 0.0207026


