.. highlight:: rest

************
Installation
************

Dependencies
============

Linetools depends on these packages:

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later.

* `numpy <http://www.numpy.org/>`_ version 1.9 or later.

* `astropy <http://www.astropy.org>`_ version 1.0 or later.

* `scipy <http://www.scipy.org/>`_ version 0.16 or later.

* `matplotlib <http://matplotlib.org/>`_  version 2 or later.

* `pyyaml <http://pyyaml.org/wiki/PyYAML>`_ version 3.11 or later.

* `specutils <https://github.com/astropy/specutils>`_ version 0.2 or later.

We strongly recommend that you use `Anaconda
<https://www.continuum.io/downloads>`_ to install them. With Anaconda
you can check for the presence and versions of the dependencies with::

  conda list "^python$|numpy|astropy$|scipy$|matplotlib|pyyaml|specutils"

If you're missing any, install them with (for example)::

  conda install astropy pyyaml matplotlib

If their versions are too old, update them with::

  conda update astropy

Specutils can't be installed with conda; instead it needs to be
installed using `pip <https://pip.pypa.io/en/latest/>`_::
  
  pip install --no-deps specutils

If you aren't using Anaconda then all of the dependencies can be
installed with pip.


Installing Linetools
====================

If you plan to play around with the code and possibly contribute
changes, then follow the :ref:`installsource` instructions
below. Otherwise simply use::

    pip install --no-deps linetools

and you're done!


Testing Linetools
=================

You can test whether your linetools installation is ok by using:

.. doctest-skip::

    >>> import linetools
    >>> linetools.test()

The tests should run and print out any failures, which we'd love you
to report at the `linetools issue tracker
<http://github.com/linetools/linetools/issues>`_.


.. _installsource:

Installing Linetools from Source
================================

*I just want to play with the code:*
------------------------------------

Install the development version like this::

    git clone https://github.com/linetools/linetools.git
    cd linetools
    python setup.py develop

Now you can easily make tweaks to the code, which are immediately
applied to your installed version (you'll have to reload the relevant
modules to see those changes in an existing Python session, though).

*I want to make a code contribution to linetools:*
--------------------------------------------------

Fantastic! In that case, follow the `Astropy developer guidelines
<http://docs.astropy.org/en/stable/development/workflow/development_workflow.html>`_,
, replacing every instance of `astropy` in those instructions with
`linetools`. This will install a 'fork' of linetools that you can use
to submit pull requests to the main repository.
