
Examples with the SolarAbund Class (v1.1)
=========================================

:download:`Download <examples/SolarAbund.ipynb>` this notebook.

.. code:: python

    # import
    from linetools.abund import solar as labsol

Init
----

.. code:: python

    sol = labsol.SolarAbund()


.. parsed-literal::

    Loading abundances from Asplund2009
    Abundances are relative by number on a logarithmic scale with H=12


.. code:: python

    sol




.. parsed-literal::

    <SolarAbund: Asplund2009>



Usage
-----

.. code:: python

    # Simple calls
    print(sol['C'])
    print(sol[6])


.. parsed-literal::

    8.43
    8.43


.. code:: python

    # Ratio
    print(sol.get_ratio('C/Fe'))


.. parsed-literal::

    0.98


.. code:: python

    # Multiple calls
    print(sol[6,12,'S'])


.. parsed-literal::

    [ 8.43  7.53  7.15]


Bits and pieces
---------------

.. code:: python

    from linetools.abund import ions as laions

.. code:: python

    # Ion name
    laions.ion_name((6,2))




.. parsed-literal::

    'CII'



.. code:: python

    # Name to ion
    laions.name_ion('CII')




.. parsed-literal::

    (6, 2)



.. code:: python

    from linetools.abund.elements import ELEMENTS


.. code:: python

    ele = ELEMENTS['C']

.. code:: python

    ele.eleconfig_dict




.. parsed-literal::

    {(1, 's'): 2, (2, 'p'): 2, (2, 's'): 2}



