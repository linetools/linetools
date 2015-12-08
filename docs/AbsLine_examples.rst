
Examples for AbsLine class (v1.1)
=================================

.. code:: python

    # import
    from linetools.spectralline import AbsLine
    from linetools.spectra import io as lsio

Generate a line
---------------

.. code:: python

    abslin = AbsLine(1548.195*u.AA)
    abslin


.. parsed-literal::

    WARNING: UnitsWarning: The unit 'Angstrom' has been deprecated in the FITS standard. Suggested: nm (with data multiplied by 0.1). [astropy.units.format.utils]
    WARNING:astropy:UnitsWarning: The unit 'Angstrom' has been deprecated in the FITS standard. Suggested: nm (with data multiplied by 0.1).


.. parsed-literal::

    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/morton03_table2.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/morton00_table2.fits.gz
    linetools.lists.parse: Reading linelist --- 
       /Users/xavier/local/Python/linetools/linetools/data/lines/verner94_tab6.fits




.. parsed-literal::

    [AbsLine: CIV 1548, wrest=1548.1950 Angstrom]



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



Measure an EW
-------------

.. code:: python

    # Set spectrum
    abslin.analy['spec'] = lsio.readspec('../../linetools/spectra/tests/files/UM184_nF.fits')


.. parsed-literal::

    /Users/xavier/local/Python/linetools/linetools/spectra/io.py:273: UserWarning: WARNING: CDELT1 < 1e-4, Assuming log wavelength scale
      warnings.warn('WARNING: CDELT1 < 1e-4, Assuming log wavelength scale')


.. code:: python

    # Set analysis range
    abslin.analy['wvlim'] = [6080.78, 6087.82]*u.AA

.. code:: python

    # Measure
    abslin.measure_ew() # Observer frame
    print('EW = {:g} with error {:g}'.format(abslin.attrib['EW'],abslin.attrib['sigEW']))


.. parsed-literal::

    EW = 0.990466 Angstrom with error 0.053301 Angstrom


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
    N, sigN, flgN = [abslin.attrib[key] for key in ['N','sigN','flagN']] 
    print('logN = {:g}, siglogN = {:g}'.format(abslin.attrib['logN'], abslin.attrib['sig_logN']))


.. parsed-literal::

    logN = 13.9037, siglogN = 0.0211602

