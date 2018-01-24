.. highlight:: rest

.. _LineList:

**************
LineList Class
**************

.. index:: LineList

Overview
========

This class organizes information about atomic and molecular transition
lines (e.g. `HI Lya`, `CIV 1548`, Hydrogen Balmer series) observed
in astrophysical environments.

..
   (:ref:`AbsLine Class`).  add this back in when written

The following lists are currently avaliable:

* 'ISM' : "All" ISM lines (can be overwhelming!)
* 'Strong' : Strong ISM lines (most common absorption line transitions observed)
* 'HI' : HI Lyman series
* 'H2' : H2 (Lyman-Werner)
* 'CO' : CO UV band-heads
* 'EUV' :  Extreme UV lines
* 'Galaxy' :  Lines typically identified in galaxy spectra
* 'AGN': Lines tipically identified in AGN spectra


Instantiation
=============

The LineList Class may be instantiated using one of the keys in the
list above::

    from linetools.lists.linelist import LineList
    hi = LineList('HI')
    #
    euv = LineList('EUV')

  
``hi``, for example, contains only HI Lyman series transitions
(e.g. `HI Lya`), and ``euv`` contains both HI Lyman series and extreme
UV metal transitions (e.g. `HI Lyb`, `NeVIII`, `MgX`).


Accessing single transitions
++++++++++++++++++++++++++++

We can now easily access atomic information regarding individual
transitions either by the rest-frame wavelength::

    wrest = 1215.67 * u.AA  # HI Lya
    hi[wrest]

or by the name convention within linetools, which in the case of `HI
Lya` is ``HI 1215``::

    name = 'HI 1215'   #  We adopt the convention of *not* rounding in the name
    hi[name]

Both cases will provide the following dictionary::

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

which summarizes the most important atomic information of `HI Lya`
transition, including the reference where these values come from
(i.e., `Morton2003`). One can therefore access any of these by
calling its dictionary keywords::

    hi['HI 1215']['wrest']
    <Quantity 1215.67 Angstrom>

is the rest-frame wavelength of the `HI Lya` transition. Similarly,::

    euv['NeVIII 780']['f']
    0.050500001758337021

is the oscillator strength of the `NeVIII 780` transition.


::::

Methods
=======

subset_lines()
++++++++++++++

This method provides a way to define a subset of lines drawn from the
original` LineList` object. Consider that for some reason you may want
only `HI Lya` and `Lyb` in your `LineList`, then you can achieve this by::

  hi = LineList('HI')
  hi = hi.subset_lines(['HI 1215', 'HI 1025'])

Which has only those two transitions loaded.

You may also want to use rest-frame wavelength to define a subset, for
instance::

    ism = LineList('ISM')
    lines = [2796.3543, 2803.5315, 1548.195, 1550.77] * u.AA
    ism = ism.subset_lines(lines)
    print(ism)
    <LineList: ISM; 4 transitions>

selects only those four transitions of `MgII` and `CIV`. In order to
avoid loading the ``LineList('ISM')`` again, you can use the keyword
`reset_data` in `subset_lines()` to make another arbitrarily different
subset of lines from the original `LineList`::

    lines = ['HI 1215', 'HI 1025']
    ism = ism.subset_lines(lines, reset_data=True)
    print(ism)
    <LineList: ISM; 2 transitions>

which now has only `HI Lya` and `Lyb`.

Finally, if you want the transitions to be sorted by rest-frame
wavelength you can use the optional keyword `sort`::

    lines = [2796.3543, 2803.5315, 1548.195, 1550.77] * u.AA
    ism = ism.subset_lines(lines, reset_data=True, sort=True)
    ism._data['wrest']
    <Quantity [ 1548.195 , 1550.77  , 2796.3543, 2803.5315] Angstrom>


set_lines()
+++++++++++

Another way to reset the LineList to its original form is by using
`set_lines()`. Following the previous example, we have a ism Linelist
with only 4 transitions::

    print(ism._data['name'])
       name
    ---------
    CIV 1548
    CIV 1550
    MgII 2796
    MgII 2803

    print(ism)
    <LineList: ISM; 4 transitions>

    ism.set_lines()
    print(ism)
    <LineList: ISM; 412 transitions>

Give us the original ism `LineList` with 412 unique transitions.

You may also want to use rest-frame wavelength to define a subset, for
instance::

    ism = LineList('ISM')
    sub_lines = [2796.3543, 2803.5315, 1548.195, 1550.77] * u.AA
    civ_mgii = ism.subset(sub_lines)

all_transitions()
+++++++++++++++++

Sometimes it may be useful to know all the transitions associated
to a given ion species. This can be achieved by the
`all_transitions()` method::

    ism = LineList('ISM')
    mgii = ism.all_transitions('MgII')

Which gives us the information of all the 6 transitions of `MgII`::

    print(mgii)
         A       el  nj  nk group    name       Ek    ...  Jk  Z   gk  gj    gamma    col0 col6
        1 / s                                  1 / cm  ...                    1 / s
    ----------- --- --- --- ----- --------- --------- ... --- --- --- --- ----------- ---- ----
      2350000.0   0   0   0     1 MgII 1025  97468.92 ... 0.0  12   4   2   2350000.0   --   --
      2480000.0   0   0   0     1 MgII 1026  97455.12 ... 0.0  12   2   2   2480000.0   --   --
      1370000.0   0   0   0     1 MgII 1239  80650.02 ... 0.0  12   4   2   1370000.0   --   --
      1540000.0   0   0   0     1 MgII 1240   80619.5 ... 0.0  12   2   2   1540000.0   --   --
    262500000.0   0   0   0     1 MgII 2796 35760.848 ... 0.0  12   4   2 262500000.0   --   --
    259500000.0   0   0   0     1 MgII 2803 35669.298 ... 0.0  12   2   2 259500000.0   --   --

In this case ``mgii`` is a Table because
more than 1 transition was found. In cases were only 1 transition
exists, the output of `all_transitions()` is a dictionary
with the same keywords as the columns of ``ism._data`` QTable::

    ciii = ism.all_transitions('CIII')
    type(ciii)
    dict
    print(ciii)
    {'A': <Quantity 1760000000.0 1 / s>,
    'Am': 0,
    'Ej': <Quantity 0.0 1 / cm>,
    'Ek': <Quantity 2352.04 1 / cm>,
    'Ex': <Quantity 0.0 1 / cm>,
    'Jj': 0.0,
    'Jk': 0.0,
    'Ref': 'Morton2003',
    'Z': 6,
    'col0': masked,
    'col6': masked,
    'el': 0,
    'f': 0.75700000000000001,
    'gamma': <Quantity 1760000000.0 1 / s>,
    'gj': 1,
    'gk': 3,
    'group': 1,
    'ion': 3,
    'mol': '',
    'name': 'CIII 977',
    'nj': 0,
    'nk': 0,
    'wrest': <Quantity 977.0201 Angstrom>}
Alternatively, we could call this method via tuple, e.g., (8,6),
where the first entry is the atomic number (of O) and the second is
the ionization state (VI) ::

    ovi = ism.all_transitions((8,6))
    print(ovi['name', 'wrest', 'f'])
      name     wrest     f
              Angstrom
    -------- --------- ------
    OVI 1031 1031.9261 0.1325
    OVI 1037 1037.6167 0.0658
You can also use a rest-frame wavelength to identify the ion species
of interest::

    wrest = 1260.4221 * u.AA
    si2 = ism.all_transitions(wrest)
    print(si2['name', 'wrest', 'f'])
       name     wrest          f
               Angstrom
    --------- --------- ---------------
    SiII 889  889.7228 0.0434000007808
    SiII 989  989.8731           0.171
    SiII 1020 1020.6989          0.0168
    SiII 1190 1190.4158           0.292
    SiII 1193 1193.2897           0.582
    SiII 1260 1260.4221            1.18
    SiII 1304 1304.3702          0.0863
    SiII 1526  1526.707           0.127
    SiII 1808 1808.0129         0.00208
    SiII 2335  2335.123        4.25e-06

For the purposes of `all_transitions`, it does not matter which
transition of a given ion species you choose, it will still retrieve
the same answer, e.g.::

    hi = ism.all_transitions('HI 1215')
    hi = ism.all_transitions('HI 1025')
    hi = ism.all_transitions(972.5367 * u.AA)
    hi = ism.all_transitions('HI')

are all equivalent. Note that in the last example we only used the
root name of the transition (i.e. the string before the blank space,
``'HI'``), so no prior knowledge of the `linetools` naming convention is
needed.


strongest_transitions()
+++++++++++++++++++++++

Sometimes it is useful to know the strongest transition for an ion in
the `LineList` within some wavelength range. `strongest_transitions()`
gives the strongest `n_max` transitions of a given ion
between a wavelength range, sorted by relative strength (defined as
the product of its rest-frame wavelength `wrest` and oscillator
strength `f`)::

    wvlims = [1000, 3000] * u.AA
    line = 'SiII'
    si2_strong = ism.strongest_transitions(line, wvlims, n_max=4)
    print(si2_strong['name'])
       name
    ---------
    SiII 1260
    SiII 1193
    SiII 1190
    SiII 1526

The syntax is the same as for `all_transitions()`. Note that you will
get the same result if you use ``line='SiII'``, ``line='SiII 1190'``,
``line='SiII 889'``, ``line=889.7228*u.AA``, or ``line=(8,6)``. By default `n_max=3`.
Depending on the wavelength range, however, the output may vary::

    wvlims = [500, 1100] * u.AA
    line = 'SiII 1260'
    si2_strong = ism.strongest_transitions(line, wvlims, n_max=4)
    print(si2_strong['name'])
       name
    ---------
    SiII 989
    SiII 889
    SiII 1020

Note that despite `n_max=4` we have only retrieved the 3 transitions
satisfying the criteria of belonging to ``wvlims = [500, 1100] * u.AA``.
Again, note that even though `SiII 1260` is out of ``wvlims`` range, it
can still be used to identify that you are interested in the `SiII` ion
species.

If you would like to retrieve all the transitions in a given ``wvlims``
regardless of its relative strength, you can set `n_max=None`.

Following the convention within `LineList`, if only 1 transition is
retrieved, the output of `strongest_transitions()` is a dictionary; if
more than 1 transition are retrieved the output is a `QTable`. If no
transition exist the output is `None`.


available_transitions()
+++++++++++++++++++++++

Sometimes it may be useful to know what are the available
transition in a given wavelength range found in the LineList
regardless of the ion species. This is particularly the case when
someone is trying to identify unknown emission/absorption lines
in a spectrum. Let us then illustrate the use of this method
with an example. Imagine that you have an observed spectrum
covering the following wavelength range::

    wvlims = [3500,5000] * u.AA

Let us now imagine that we are interested in a particular redshift, say
``z=0.67``. Then, we can do::

    z = 0.67
    transitions = ism.available_transitions(wvlims/(1+z), n_max_tuple=None, min_strength=0.)
    print(len(transitions))
    33

Will give the 33 transitions available that could correspond to having
``z=0.67`` in the form of a `QTable`. The output is sorted by strength of
the strongest available transition per ion species, and strength is defined
as `log10(wrest * fosc * abundance)`, where `abundance` is that of the solar
composition given by Asplund2009. As optional keyword parameters one can
specify a minimum strength as `min_strength`, so transitions below this
value are omitted, e.g.::

    transitions = ism.available_transitions(wvlims/(1+z), n_max_tuple=None, min_strength=10.5)
    print(len(transitions))
    3

Which correspond to `MgI 2852`, `MgII 2796` and `MgII 2803`. Note than this
method does not correct for ionization state. Similarly, one can also specify
the maximum number of transitions per ion species
tuple using the optional keyword parameter `n_max_tuple`, e.g.::

    transitions = ism.available_transitions(wvlims/(1+z), n_max_tuple=1, min_strength=0.)
    print(transitions['name'])
        name
    -----------
    MgI 2852
    MgII 2796
    FeII 2382
    FeII* 2396b
    MnII 2576
    VII  2683
    ...

Which for the case of `MgII` only retrieves ``'MgII 2796'``. Again, following the convention within
`LineList`, if only 1 transition is retrieved, the output of `available_transitions()`
is a dictionary; if more than 1 transition are retrieved the output is a `QTable`. If no
transition exist satisfying the criteria the output is `None`.
