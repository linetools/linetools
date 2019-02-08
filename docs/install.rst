************
Installation
************

Dependencies
============

Linetools depends on these packages:

* `python <http://www.python.org/>`_ version 3.6 or later
* `numpy <http://www.numpy.org/>`_ version 1.12 or later
* `astropy`_ version 3.0 or later
* `scipy <http://www.scipy.org/>`_ version 0.19 or later
* `matplotlib <http://matplotlib.org/>`_  version 2.0 or later
* `PyQT5 <https://wiki.python.org/moin/PyQt>`_ version 5 (for GUIs)
* `h5py <https://www.h5py.org/>`_ version 2.6 (for data I/O)

We strongly recommend that you use `Anaconda
<https://www.continuum.io/downloads>`_ to install them. With Anaconda
you can check for the presence and versions of the dependencies with::

  conda list "^python|numpy|astropy|scipy|matplotlib|PyQT|h5py"

If you're missing any, install them with (for example)::

  conda install astropy scipy matplotlib h5py

For PyQT, the current conda version is PyQT5.  If you are still using
PyQT4, then you will need to update to use linetools GUIs:

    conda install pyqt=5

If their versions are too old, update them with (for example)::

  conda update astropy

If you aren't using Anaconda then all of the dependencies can also be
installed with pip.

Also note, if you wish to have full functionality of the GUIs and are
using MacOSX, then you probably need to change
your *backend* from macosx to TkAgg in the matplotlibrc file.

Installing Linetools
====================

pip
---

Do *not* use the version of linetools on PyPi (i.e. pip install).
For now, see :ref:`installsource` to install linetools from source.

.. _installsource:

Installing Linetools from Source
--------------------------------

There are two modes to consider and choose from:

*I just want to play with the code*
+++++++++++++++++++++++++++++++++++

Install the development version like this::

    git clone https://github.com/linetools/linetools.git
    cd linetools
    python setup.py develop

Now you can easily make tweaks to the code, which are immediately
applied to your installed version (you'll have to reload the relevant
modules to see those changes in an existing Python session, though).

*I want to make a code contribution to linetools*
+++++++++++++++++++++++++++++++++++++++++++++++++

Fantastic! In that case, follow the `Astropy developer guidelines
<http://docs.astropy.org/en/stable/development/workflow/development_workflow.html>`_,
replacing every instance of **astropy** in those instructions with
**linetools**. This will install a 'fork' of linetools that you can
use to submit code changes to the main repository.


Running Tests
=============

To test your installation, run::

    python -c 'import linetools; linetools.test()'

The tests take a couple of minutes to finish. If you notice any
failures, we'd love you to report them on the `linetools issue tracker
<http://github.com/linetools/linetools/issues>`_.

Before Launching GUIs
=====================

If you are a Mac user, we **highly** recommend that you set your
matplotlib backend from MacOSX to TkAgg (or another option, see
`backends <http://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`__).


Building Documentation
======================

Only do this if you're a developer! If you want build the
documentation, you also need to install Sphinx (version 1.3+)::

  conda install sphinx

If you'd like to generate inheritance diagrams in the docs then you
also need to install graphviz (`MacOSX
<http://www.graphviz.org/Download_macos.php>`_, `Ubuntu
<http://www.graphviz.org/Download_linux_ubuntu.php>`_), but this isn't
required. Once sphinx is installed, change to the `/docs` directory
under the source directory and run::

  make html

The documentation should now be in _build/html.
