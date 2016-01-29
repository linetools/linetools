****************
linetool Scripts
****************

lt_xspec
--------

.. toctree::
   :maxdepth: 1

   Example Notebook <xspecgui>

Launch a QT Gui to examine an input spectrum.

For example::

   lt_xspec filename

You can explore the data, perform simple analysis (e.g. EW
measurements) and more. See the Notebook for more.

.. warning::
 
  The xspec gui is still in an experimental stage and some aspects
  might not work as expected.


lt_plot
-------

Plot one or more spectra. This has fewer features than lt_xspec above,
but is faster.

For example::

	lt_plot filename1 filename2

A plot of the spectrum in ``filename1`` appears, and you can navigate
around it using the same commands as in lt_xspec. Move to the next or
previous spectrum with the right or left arrow key.

To list all the command line options available, use::

        lt_plot -h

lt_absline
----------

This plots a single absorption line, given a transition rest
wavelength (Angstroms), log10 column density, and Doppler parameter
(km/s).

For example::

	lt_absline 1215.6701 14.0 30

plots a Hydrogen Ly-a line with column density of 10\ :sup:`14` cm\
:sup:`-2` and b=30 km/s. A plot will appear and the line info and EW
as well, i.e. ::

	[AbsLine: HI 1215, wrest=1215.6700 Angstrom]
	EW = 0.268851 Angstrom

Try:: 

	lt_absline -h

for the full set of options.
