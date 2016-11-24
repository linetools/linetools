.. PYPIT documentation master file, created by
   sphinx-quickstart on Fri Nov 13 13:39:35 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Linetools
=========

linetools is an in-development package for the analysis of 1-D
spectra, with the aim to become an `Astropy`_ `affiliated package
<http://www.astropy.org/affiliated/index.html>`_. Its core developers
work primarily on UV/optical/IR absorption line research, so most of the
functionality is aimed at the identification and analysis of
absorption lines. The eventual goal is to provide a set of tools
useful for both absorption and emission lines.


.. note::

    linetools is still under active development. While the developers
    strive to maintain compatibility in new releases, there may be
    backwards-incompatible changes in future versions.


Getting Started
---------------

.. toctree::
   :maxdepth: 1

   install
   changelog

Core classes
------------

.. toctree::
   :maxdepth: 1

   AbsComponent <abscomp>
   AbsSystem <abssys>
   AbsSightline <abssightline>
   RelAbund <relabund>
   SolarAbund <solar>
   LineList <linelist>
   SpectralLine <specline>
   XSpectrum1D <xspectrum1d>

Graphical User Interfaces (GUIs)
--------------------------------

.. toctree::
   :maxdepth: 2

   guis.rst

Command line tools
------------------

.. toctree::
   :maxdepth: 2

   Scripts <scripts>


Reference & API
---------------

.. toctree::
   :maxdepth: 1

   api


..
   **Project details**

   .. toctree::
      :maxdepth: 1

      whatsnew
      credits
      license


Indices and tables
++++++++++++++++++

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
