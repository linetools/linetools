************
Installation
************

Dependencies
============

Linetools depends on these packages:

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later
* `numpy <http://www.numpy.org/>`_ version 1.9 or later
* `astropy`_ version 1.0 or later
* `scipy <http://www.scipy.org/>`_ version 0.16 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.4 or later
* `specutils <https://github.com/astropy/specutils>`_ version 0.2 or later
* `PyQT4 <https://wiki.python.org/moin/PyQt>`_ version 4 (for GUIs)

We strongly recommend that you use `Anaconda
<https://www.continuum.io/downloads>`_ to install them. With Anaconda
you can check for the presence and versions of the dependencies with::

  conda list "^python$|numpy|astropy$|scipy$|matplotlib|specutils|PyQT"

If you're missing any, install them with (for example)::

  conda install astropy scipy matplotlib PyQT

If their versions are too old, update them with (for example)::

  conda update astropy

Specutils can't be installed with conda; use `pip
<https://pip.pypa.io/en/latest/>`_ instead::
  
  pip install --no-deps specutils

If you aren't using Anaconda then all of the dependencies can also be
installed with pip.


Installing Linetools
====================

If you plan to play around with the code and possibly contribute
changes, then follow the instructions in the section below,
:ref:`installsource`. Otherwise simply use::

    pip install linetools

and you're done!

If you wish to have full functionality of the GUIs and are
using MacOSX, then you probably need to change
your *backend* from macosx to TkAgg in the matplotlibrc file.

.. _installsource:

Installing Linetools from Source
================================

*I just want to play with the code*
-----------------------------------

Install the development version like this::

    git clone https://github.com/linetools/linetools.git
    cd linetools
    python setup.py develop

Now you can easily make tweaks to the code, which are immediately
applied to your installed version (you'll have to reload the relevant
modules to see those changes in an existing Python session, though).

*I want to make a code contribution to linetools*
-------------------------------------------------

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
