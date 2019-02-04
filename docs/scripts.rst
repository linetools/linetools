****************
linetool Scripts
****************

There are a number of scripts, many of which are GUIs,
provided with linetools.  As regards the GUIs we warn
again that Mac users will need to set their matplotlib to
something other than MacOSX. See
`backends <http://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`__.

lt_xspec
--------

.. toctree::
   :maxdepth: 1

   Example Notebook <xspecgui>

Launch a QT Gui to examine an input spectrum.

For example::

   lt_xspec filename.fits
   lt_xspec filename.fits#1#   -- Specifies extension 1 for a multi-extension FITS file of binary tables

You can explore the data, perform simple analysis (e.g. EW
measurements) and more. See the Notebook for more.

Here is the full usage::

    usage: lt_xspec [-h] [-z ZSYS] [--norm] [--air] [--exten EXTEN]
                    [--splice SPLICE] [--scale SCALE]
                    file

    Parser for lt_xspec v1.2; Note: Extra arguments are passed to read_spec (e.g.
    --flux_tag=FX)

    positional arguments:
      file                  Spectral file; specify extension by appending #exten#

    optional arguments:
      -h, --help            show this help message and exit
      -z ZSYS, --zsys ZSYS  System Redshift
      --norm                Show spectrum continuum normalized (if continuum is
                            provided)
      --air                 Convert input spectrum wavelengths from air to vacuum
      --exten EXTEN         FITS extension
      --splice SPLICE       Splice with the input file; extension convention
                            applies
      --scale SCALE         Scale factor for GUI size [1. is default]



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
wavelength (Angstroms) or name (e.g. CIV1548),
log10 column density, and Doppler parameter
(km/s).

For example::

	lt_absline 1215.6701 14.0 30
	lt_absline HI1215 14.0 30

plots a Hydrogen Ly-a line with column density of 10\ :sup:`14` cm\
:sup:`-2` and b=30 km/s. A plot will appear and the line info and EW
as well, i.e. ::

	[AbsLine: HI 1215, wrest=1215.6700 Angstrom]
	EW = 0.268851 Angstrom

Try:: 

	lt_absline -h

for the full set of options.

lt_radec
--------

Input coordinates in one format and print them out
in several formats::

   lt_radec 152.25900,7.22885
   lt_radec J100902.16+071343.86
   lt_radec 10:09:02.16,+07:13:43.8


lt_line
-------

Print the atomic data for an input ion, transition or for an
entire linelist.::

    lt_line -h
    lt_line HI
    lt_line HI -z 3.0
    lt_line HI1215
    lt_line 1215
    lt_line --all

Here is the usage::

    beast> lt_line -h
    usage: lt_line [-h] [-a] [--llist LLIST] [-t TOLER] [-z REDSHIFT] [inp]

    Print spectral line data of a line or lines.

    positional arguments:
      inp                   Ion, transition name, or Rest wavelength (e.g. HI,
                            HI1215 or 1215)

    optional arguments:
      -h, --help            show this help message and exit
      -a, --all             Print all lines
      --llist LLIST         Name of LineList: ISM, HI, H2, CO, etc.
      -t TOLER, --toler TOLER
                            Matching tolerance (in Ang) on input wavelength
      -z REDSHIFT, --redshift REDSHIFT
                            Matching tolerance (in Ang) on input wavelength


lt_solabnd
----------

Print the solar abundances by number relative to Hydrogen on the
traditional log scale of 12, e.g. e(Fe) = 7.55::

   lt_solabnd Fe
   lt_solabnd -a
   lt_solabnd -a --sortZ

Here is the usage::

    beast> lt_solabnd -h
    usage: lt_solabnd [-h] [-a] [--sortZ] [inp]

    Print Solar abundance data for an element or all elements.

    positional arguments:
      inp         Elm (e.g. H, Fe)

    optional arguments:
      -h, --help  show this help message and exit
      -a, --all   Print all values
      --sortZ     Sort on Atomic Number


lt_continuumfit
---------------

Launch the GUI to continuum fit a spectrum.
If a redshift is supplied by zsys, then the
script assumes this is a QSO.  Enables spectra
to be read in via a `specdb <https://specdb.readthedocs.io/en/latest/>`__ file
which then requires specdb be installed.::

   lt_continuumfit input_file output_filename --redshift 0.867

Here is the current usage message::

    usage: lt_continuumfit [-h] [--redshift REDSHIFT] [--wchunk WCHUNK] [--native]
                           [--specdb SPECDB]
                           file outfil

    GUI to fit a continuum to a spectrum

    positional arguments:
      file                 Input spectral file (FITS, ASCII, etc.)
      outfil               Output, normalized spectrum filename; FITS [can be the
                           same]

    optional arguments:
      -h, --help           show this help message and exit
      --redshift REDSHIFT  Redshift of the Source
      --wchunk WCHUNK      Width of a 'chunk' (Ang)
      --native             Do not mask input spectrum
      --specdb SPECDB      Input file is specdb. Input (ra,dec,group) in this
                           order without spaces

