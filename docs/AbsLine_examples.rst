
Examples for AbsLine class (v1.2)
=================================

Download :download:`examples/AbsLine_examples.ipynb` this notebook.

.. code:: python

    # import
    from linetools.spectralline import AbsLine, SpectralLine
    from linetools import spectralline as ltsp
    
    from linetools.spectra import io as lsio

Generate a line
---------------

.. code:: python

    abslin = AbsLine(1548.195*u.AA)
    abslin


.. parsed-literal::

    WARNING: UnitsWarning: The unit 'Angstrom' has been deprecated in the FITS standard. Suggested: 10**-1 nm. [astropy.units.format.utils]
    WARNING:astropy:UnitsWarning: The unit 'Angstrom' has been deprecated in the FITS standard. Suggested: 10**-1 nm.


.. parsed-literal::

    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/morton03_table2.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/morton00_table2.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/verner96_tab1.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/verner94_tab6.fits

.. parsed-literal::

    WARNING: UnitsWarning: '0.1nm' did not parse as fits unit: Numeric factor not supported by FITS [astropy.units.core]
    WARNING:astropy:UnitsWarning: '0.1nm' did not parse as fits unit: Numeric factor not supported by FITS


.. parsed-literal::

    
    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/EUV_lines.ascii
    read_sets: Using set file -- 
      /Users/xavier/local/Python/linetools/linetools/lists/sets/llist_v1.0.ascii




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

    {'analy': {u'datafile': u'',
      u'do_analysis': 1,
      u'flg_eye': 0,
      u'flg_limit': 0,
      u'name': 'CIV 1548',
      u'vlim': {'unit': u'km / s', 'value': [0.0, 0.0]},
      u'wvlim': {'unit': u'Angstrom', 'value': [0.0, 0.0]}},
     'attrib': {u'DEC': 0.0,
      u'EW': {'unit': u'Angstrom', 'value': 0.0},
      u'N': {'unit': u'1 / cm2', 'value': 0.0},
      u'RA': 0.0,
      u'b': {'unit': u'km / s', 'value': 0.0},
      u'flag_EW': 0,
      u'flag_N': 0,
      u'sig_EW': {'unit': u'Angstrom', 'value': 0.0},
      u'sig_N': {'unit': u'1 / cm2', 'value': 0.0},
      u'sig_b': {'unit': u'km / s', 'value': 0.0},
      u'sig_v': {'unit': u'km / s', 'value': 0.0},
      u'sig_z': 0.0,
      u'v': {'unit': u'km / s', 'value': 0.0},
      u'z': 0.0},
     'data': {'A': {'unit': u'1 / s', 'value': 0.0},
      'Am': 0,
      'Ej': {'unit': u'1 / cm', 'value': 0.0},
      'Ek': {'unit': u'1 / cm', 'value': 0.0},
      'Ex': {'unit': u'1 / cm', 'value': 0.0},
      'Jj': 0.0,
      'Jk': 0.0,
      'Ref': 'Verner1994',
      'Z': 6,
      'el': 0,
      'f': 0.18999999761581421,
      'gamma': {'unit': u'1 / s', 'value': 0.0},
      'gj': 2,
      'gk': 4,
      'group': 1,
      'ion': 4,
      'mol': '',
      'name': 'CIV 1548',
      'nj': 0,
      'nk': 0,
      'wrest': {'unit': u'Angstrom', 'value': 1548.195}},
     'ltype': u'Abs',
     'trans': 'CIV 1548',
     'wrest': {'unit': u'Angstrom', 'value': 1548.195}}



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
    abslin.analy['spec'] = lsio.readspec('../../linetools/spectra/tests/files/UM184_nF.fits')

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

    logN = 13.9051, siglogN = 0.0207027


