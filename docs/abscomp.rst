.. _AbsComponent:

******************
AbsComponent Class
******************

.. index:: AbsComponent

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <AbsComponent_examples>
   Column Densities <AbsComponent_ColumnDensities>

Overview
========

This Class is designed to organize and analyze a set of
absorption lines.

By definition, an AbsComponent is a unique collection of
absorption lines specified by:

=============== ========   ============== ============================================
Property        Variable   Type           Description
=============== ========   ============== ============================================
RA, DEC         radec      tuple or coord RA,DEC in deg or astropy.coordinate
Z, ion          Zion       tuple          Atomic Number (Z) and ionization state (ion)
Redshift        zcomp      float          Absorption redshift of the component
Velocity limits vlim       Quantity array -/+ velocity limits of the component
Energy level    Ej         Quantity       Energy of the level relative to ground
Isotope         A          int            Atomic Mass number (optional)
=============== ========   ============== ============================================


Instantiation
=============

The AbsComponent Class may be instantiated in a few ways.
The default sets the properties listed above::

	abscomp = AbsComponent((10.0*u.deg, 45*u.deg), (14,2), 1.0, [-300,300]*u.km/u.s)

More commonly, one will instantiate with one
`~linetools.spectralline.AbsLine` object::

    lya = AbsLine(1215.670*u.AA, z=2.92939)
    lya.limits.set([-300.,300.]*u.km/u.s)  # vlim
    abscomp1 = AbsComponent.from_abslines([lya])

or multiple::

    lyb = AbsLine(1025.7222*u.AA, z=lya.z)
    lyb.limits.set([-300.,300.]*u.km/u.s)  # vlim
    abscomp = AbsComponent.from_abslines([lya,lyb])

One may also instantiate from a *dict*, usually read
from the hard-drive (e.g. a JSON file)::

    abscomp = AbsComponent.from_dict(idict)

::::

One may also generate a set of components from a larger
list of AbsLines::

   complist = ltiu.build_components_from_abslines([lya,lyb,SiIIlines[0],SiIIlines[1]])

Attributes
==========

There is a set of default attributes that are initialized at Instantiation
and kepy in an attrib *dict*.  These are::

    init_attrib = {'N': 0./u.cm**2, 'sig_N': [0.,0.]/u.cm**2, 'flag_N': 0, # Column
              'logN': 0., 'sig_logN': np.array([0.,0.]),
              'b': 0.*u.km/u.s, 'sig_b': 0.*u.km/u.s,  # Doppler
              'vel': 0*u.km/u.s, 'sig_vel': 0*u.km/u.s
              }


One can access these attributes with standard . syntax, e.g.::

   logN = abscomp.logN   # This accesses the logN value from the internal attrib dict

Inspecting
==========

Here are a few simple methods to explore/inspect the class.

Generate a Table
++++++++++++++++

If the class contains one or more AbsLines, you may generate a
`~astropy.table.Table` from their attributes and data::

    comp_tbl = abscomp.build_table()

Show a Stack Plot
+++++++++++++++++

If the AbsLine have spectra attached to them (in attrib['spec']),
a stack plot (aka velocity plot) is generated with::

    abscomp.stack_plot()

Apparent Column Densitities
+++++++++++++++++++++++++++

Show a plot of the apparent column density profiles, :math:`N_a`::

    abscomp.plot_Na()

::::

Analysis
========

Here are some methods related to analysis.

Synthesize Columns
++++++++++++++++++


If one inputs a set of AbsLine(s) with column density measurements,
the column densities are synthesized at instantiation unless one
sets skip_synth=True.

There are two approaches to synthesis:

Standard
--------

The default uses the synthesize_colm() method.
Positive, unsaturated detections
are combined in a weighted mean whereas limits are compared
and the most sensitive one is adopted.::

Here is the set of rules:

1.  If all measurements are upper limits, take the lowest value and flag as an upper limit (*flag_N=3*).
2.  If all measurements are a mix of upper and lower limits, take the highest lower limit and flag as a lower limit (*flag_N=2*).
3.  If one or more measurements are a proper detection, take the weighted mean of these and flag as a detection (*flag_N=1*).

A future implementation will introduce a flag for
bracketing an upper and lower limit value.

Median
------

Another option for synthesizing the column densities and other attributes
is to take the median values.  This is performed when one
passes adopt_median=True to the from_abslines() call, e.g.::

    abscomp = AbsComponent.from_abslines([lya,lyb], adopt_median=True)


The current code computes the median on all input
values and does not consider the flag_N values for
the input AbsLine objects.

Curve of Growth
+++++++++++++++

A standard, single-component curve-of-growth (COG) analysis may be
performed on the set of AbsLines::

    COG_dict = abscomp.cog(show_plot=True)

The output dict includes:

========== ============== =====================================
Key        Type           Description
========== ============== =====================================
EW         Quantity array Input equivalent widths
sigEW      Quantity array Input error in equivalent widths
f          ndarray        Input f-values
wrest      Quantity array Input rest wavelengths
logN       float          Output fitted column density (log10)
sig_logN   float          Output error in fitted logN
b          Quantity       Output b-value (km/s)
sig_b      Quantity       Output error in b-value (km/s)
========== ============== =====================================

Misc
====

I/O
+++

One may generate a *dict* of the key properties of the AbsComponent
with the to_dict() method::

   cdict = component.to_dict()

One may also wish to Voigt profile fit components with one
of a number of software packages (e.g. ALIS, JoeBVP).  To
generate an input file for JoeBVP use::

    from linetools.isgm.io import write_joebvp_from_components
    write_joebvp_from_components(component_list, specfile, 'output_file.ascii')

Similarly, one can generate a list of components from an outputted
JoeBVP file::

    from linetools.isgm.io import read_joebvp_to_components
    comp_list = read_joebvp_to_components('joebvp_file.out')


Synthesize Components
+++++++++++++++++++++

This method combines a list of two or more components into a new one.
It checks first for consistent RA/DEC, Zion, and Ej.  It does
not place any constraints on z and vlim.  The column density of
the new component is the sum of the input ones (with rules for
limits).  And the redshift and vlim are set to encompass the
velocities of the input components.::

   from linetools.isgm import utils as ltiu
   synth_SiII = ltiu.synthesize_components([SiIIcomp1,SiIIcomp2])

See the :doc:`AbsComponent_examples` notebook for a complete example.


Tables
++++++

You can also create a list of components using an input `astropy.table.Table` object
(mandatory column names are
`['RA', 'DEC', 'ion_name', 'z_comp', 'vmin', 'vmax', 'Z', 'ion', 'Ej']`)::

   from astropy.table import Table
   tab = Table()
   tab['ion_name'] = ['HI', 'HI', 'CIV', 'SiII', 'OVI']
   tab['Z'] = [1,1,4,14,8]
   tab['ion'] = [1,1,4,2,6]
   tab['Ej'] = [0.,0.,0.,0.,0.]/u.cm
   tab['z_comp'] = [0.05, 0.0999, 0.1, 0.1001, 0.3]
   tab['logN'] = [12.0, 13.0., 13.0, 14.0, 13.5]
   tab['sig_logN'] = [0.1, 0.12., 0.15, 0.2, 0.1]
   tab['flag_logN'] = [1, 1, 2, 0, 1]
   tab['RA'] = 100.0 * u.deg * np.ones(len(tab))
   tab['DEC'] = -0.8 * u.deg * np.ones(len(tab))
   tab['vmin'] = -50 * u.km/u.s * np.ones(len(tab))
   tab['vmax'] =  50 * u.km/u.s * np.ones(len(tab))

   complist = ltiu.complist_from_table(tab)

These components will not have AbsLines defined by default.
However, it is very easy to append
AbsLines to these components for a given observed wavelength
and minimum rest-frame equivalent width::

   # add abslines
   wvlim = [1100, 1500]*u.AA
   for comp in complist:
      comp.add_abslines_from_linelist(llist='ISM', wvlim=wvlim, min_Wr=0.01*u.AA, chk_sep=False)

This will look for transitions in the `LineList('ISM')` object and append those
expected to be within `wvlim`
to each corresponding AbsComponent in the `complist`. In this example, only
those AbsLines expected to have rest-frame
equivalent widths larger than `0.01*u.AA` will be appended, but if you wish to
include all of the available AbsLines
you can set `min_Wr=None`. Here, we have set `chk_sep=False` to avoid checking
for coordinates because by construction the
individual AbsLines have the same coordinates as the corresponding AbsComponent.

One can, of course, go the other way and generate a Table
from a list of components::

    tab2 = ltiu.table_from_complist(complist)

If you wish to write to disk, we recommend that you use
the *astropy* format: ascii.ecsv


Use Components to create a spectrum model
+++++++++++++++++++++++++++++++++++++++++

You can create a spectrum model using a list of AbsComponents, e.g.::

   from linetools.analysis import voigt as lav
   wv_array = np.array(1100,1500, 0.1)*u.AA
   model = lav.voigt_from_components(wv_array, complist)

In this manner, model is a XSpectrum1D object with the the AbsComponents contained in the `complist` list.
You may also add a convolution with a given kernel to compare with observations from different spectrographs,
in which case you can use::

   model = lav.voigt_from_components(wv_array, complist, fwhm=3)

This is same as before but the model is convolved with a Gaussian kernel of FWHM of 3 pixels.

