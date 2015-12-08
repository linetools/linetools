.. highlight:: rest

.. _LineList:

******************
LineList Class
******************

.. index:: LineList

Overview
========

This Class is designed to organize and handle information about atomic
and/or molecular transition lines (e.g. HI Lya, CIV 1548, Hydrogen
Balmer series, etc.) observed in a variety of astrophysical
environments. It is currently implemented for absorption lines, but we
expect to also include common emission lines in the near future.

..
   (:ref:`AbsLine Class`).  add this back in when written

By definition, a LineList is a unique collection of transitions lines
specified by:


============ ========= =========== ====================================================
Property     Variable  Type        Description
============ ========= =========== ====================================================
LineList key llst_keys str or list A key to define a subsample of transitions to load:
                                     'ISM' - "All" ISM lines (can be overwhelming!)
                                     'Strong' - Strong ISM lines
                                     'HI' - HI Lyman series
                                     'H2' - H2 (Lyman-Werner)
                                     'CO' - CO UV band-heads
                                     'EUV' - Extreme UV lines
============ ========= =========== =====================================================


Instantiation
=============

The LineList Class may be instantiated using a single key::

	hi = LineList('HI')

  linetools.lists.parse: Reading linelist --- 
   /home/ntejos/python/linetools/linetools/data/lines/morton03_table2.fits.gz
  WARNING: UnitsWarning: The unit 'Angstrom' has been deprecated in the FITS standard. Suggested: nm (with data multiplied by 0.1). [astropy.units.format.utils]
  read_sets: Using set file -- 
   /home/ntejos/python/linetools/linetools/lists/sets/llist_v0.4.ascii
  
or a list of keys:

  euv = LineList(['HI','EUV'])

  linetools.lists.parse: Reading linelist --- 
   /home/ntejos/python/linetools/linetools/data/lines/morton03_table2.fits.gz
  linetools.lists.parse: Reading linelist --- 
   /home/ntejos/python/linetools/linetools/data/lines/morton00_table2.fits.gz
  linetools.lists.parse: Reading linelist --- 
   /home/ntejos/python/linetools/linetools/data/lines/verner94_tab6.fits
  linetools.lists.parse: Reading linelist --- 
   /home/ntejos/python/linetools/linetools/data/lines/EUV_lines.ascii
  read_sets: Using set file -- 
  /home/ntejos/python/linetools/linetools/lists/sets/llist_v0.4.ascii

In these examples, the object ``hi`` has purely HI Lyman series
transitions (e.g. HI Lya) and ``euv`` has HI Lyman series and Extreme
UV transitions (e.g. Ne VIII, MgX). The available keys are listed in
Table XX.

We can now easily access atomic information regarding individual
transitions either by the rest-frame wavelength::

  wrest = 1215.67 * u.AA  # HI Lya
  hi[wrest]

or by the name convention within linetools, which in the case of HI
Lya is ``HI 1215``::

  name = 'HI 1215'
  hi[name]

both cases will provide the following dictionary::

  {'A': <Quantity 626500000.0 1 / s>,   # Einstein coefficient
  'Am': 0,                              # Mass number (often written as "A"; only used for D) 
  'Ej': <Quantity 0.0 1 / cm>,          # Energy of lower level (relative to ground state)
  'Ek': <Quantity 2259.163 1 / cm>,     # Energy of upper level (relative to ground state)
  'Ex': <Quantity 0.0 1 / cm>,          # Excitation energy (cm^-1)
  'Jj': 0.0,                            # Tot ang mom (z projection) of lower state (or rotation level)
  'Jk': 0.0,                            # Tot ang mom (z projection) of upper state (or rotation level)
  'Ref': 'Morton2003',                  # References
  'Z': 1,                               # Atomic number (for atoms)       
  'col0': masked,                       # (Reserved)
  'col6': masked,                       # (Reserved)
  'el': 0,                              # Electronic transition (2=Lyman (B-X), 3=Werner (C-X)) 
  'f': 0.41639999999999999,             # Oscillator strength
  'gamma': <Quantity 626500000.0 1 / s>,# Sum of A 
  'gj': 2,                              # Lower statistical weight (2J+1)
  'gk': 6,                              # Upper statistical weight (2J+1)
  'group': 1,                           # Flag for grouping
  'ion': 1,                             # Ionic state (1=Neutral)
  'mol': '',                            # Molecular name (H2, HD, CO, C13O)
  'name': 'HI 1215',                    # Name
  'nj': 0,                              # Orbital level of lower state (or vibrational level)
  'nk': 0,                              # Orbital level of upper state (or vibrational level)
  'wrest': <Quantity 1215.67 Angstrom>} # Rest Wavelength (Quantity)  

which summarizes the most important atomic information of HI Lya
transition, including the reference where these values come from
(i.e., ``Morton2003``). One can therefore access any of these by
calling its dictionary keywords::

  hi['HI 1215']['wrest']
  <Quantity 1215.67 Angstrom>

or::

  euv['NeVIII 780']['f']
  0.050500001758337021

