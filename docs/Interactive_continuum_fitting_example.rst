
Interactive continuum fitting
=============================

:download:`Download <examples/Interactive_continuum_fitting_example.ipynb>` this notebook.

.. code:: python

    # suppress warnings for these examples
    import warnings
    warnings.filterwarnings('ignore')
    
    import imp
    prefix = imp.find_module('linetools')[1] + '/spectra/tests/files/'
    import linetools.spectra.xspectrum1d as lsx
    spec = lsx.XSpectrum1D.from_file(prefix + 'q0002m422.txt.gz')
    # keep the old continuum to compare later on
    co_old = spec.co.copy()

.. code:: python

    %pylab


.. parsed-literal::

    Using matplotlib backend: TkAgg
    Populating the interactive namespace from numpy and matplotlib


.. code:: python

    # now fit the continuum interactively. We say we're fitting a QSO, 
    # so it can make intelligent guesses about where to put the spline
    # points that define the continuum.
    spec.fit_continuum(kind='QSO', redshift=2.76)
    
    # now you can interactively tweak these spline points, adding or
    # removing them as necessary. Once you're finished, press 'q' to
    # close the window.


.. parsed-literal::

    knots file exists, use this? (y) n
    
    i,o          Zoom in/out x limits
    y            Zoom out y limits
    Y            Guess y limits
    t,b          Set y top/bottom limit
    l,r          Set left/right x limit
    [,]          Pan left/right
    w            Plot the whole spectrum
    
    S,U          Smooth/unsmooth spectrum
    
    
    i,o          Zoom in/out x limits
    y            Zoom out y limits
    Y            Guess y limits
    t,b          Set y top/bottom limit
    l,r          Set left/right x limit
    [,]          Pan left/right
    w            Plot the whole spectrum
    
    S,U          Smooth/unsmooth spectrum
    
    a        : add a new spline knot
    A        : add a new spline knot, and use a flux median to guess y position
    +        : double the number of spline knots
    _        : halve the number of spline knots
    d        : delete the nearest knot
    m        : move the nearest knot
    M        : move the nearest knot, and use a flux median to guess y position
    
    q        : quit
    
    Updating continuum


.. code:: python

    # the New continuum is now saved in spec.co, and the spline knots are in
    # spec.meta['contpoints']
    #
    # Let's compare the old and new continuum
    plt.figure()
    wa = spec.dispersion.value
    plt.plot(wa, co_old)
    plt.plot(wa, spec.co)




.. parsed-literal::

    [<matplotlib.lines.Line2D at 0x10c8df978>]



.. code:: python

    co_old2 = spec.co.copy()
    
    # we can also tweak a small section of the continuum without affecting the whole spectrum.
    spec.fit_continuum(wlim=(5000, 5100))


.. parsed-literal::

    knots file exists, use this? (y) n
    
    i,o          Zoom in/out x limits
    y            Zoom out y limits
    Y            Guess y limits
    t,b          Set y top/bottom limit
    l,r          Set left/right x limit
    [,]          Pan left/right
    w            Plot the whole spectrum
    
    S,U          Smooth/unsmooth spectrum
    
    
    i,o          Zoom in/out x limits
    y            Zoom out y limits
    Y            Guess y limits
    t,b          Set y top/bottom limit
    l,r          Set left/right x limit
    [,]          Pan left/right
    w            Plot the whole spectrum
    
    S,U          Smooth/unsmooth spectrum
    
    a        : add a new spline knot
    A        : add a new spline knot, and use a flux median to guess y position
    +        : double the number of spline knots
    _        : halve the number of spline knots
    d        : delete the nearest knot
    m        : move the nearest knot
    M        : move the nearest knot, and use a flux median to guess y position
    
    q        : quit
    
    Updating continuum


.. code:: python

    # check it works without a predefined continuum
    spec = lsx.XSpectrum1D.from_file(prefix + 'q0002m422.txt.gz')
    spec.co = None
    spec.fit_continuum(kind='QSO', redshift=2.76)


.. parsed-literal::

    knots file exists, use this? (y) n
    
    i,o          Zoom in/out x limits
    y            Zoom out y limits
    Y            Guess y limits
    t,b          Set y top/bottom limit
    l,r          Set left/right x limit
    [,]          Pan left/right
    w            Plot the whole spectrum
    
    S,U          Smooth/unsmooth spectrum
    
    
    i,o          Zoom in/out x limits
    y            Zoom out y limits
    Y            Guess y limits
    t,b          Set y top/bottom limit
    l,r          Set left/right x limit
    [,]          Pan left/right
    w            Plot the whole spectrum
    
    S,U          Smooth/unsmooth spectrum
    
    a        : add a new spline knot
    A        : add a new spline knot, and use a flux median to guess y position
    +        : double the number of spline knots
    _        : halve the number of spline knots
    d        : delete the nearest knot
    m        : move the nearest knot
    M        : move the nearest knot, and use a flux median to guess y position
    
    q        : quit
    
    Updating continuum


