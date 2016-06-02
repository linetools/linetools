
xabssys Documentation
=====================

.. :download:`Download <examples/xspecgui.ipynb>` this notebook.

This ipython Notebook is intended to provide documentation for the
linetools GUI named XAbsSysGui.

Enjoy and feel free to suggest edits/additions, etc.


Before Launching the GUI
------------------------

If you are a Mac user, we **highly** recommend that you set your
matplotlib backend from MacOSX to TkAgg (or another option, see
`backends <http://matplotlib.org/faq/usage_faq.html#what-is-a-backend>`__).

Launching the GUI
-----------------

From the command line (recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We recommend you use the script provided with linetools.

Then it is as simple as:

::

    > lt_xabssys spec_file abssys_file

The abssys_file is expected to be a JSON file that contains
an AbsSystem (likely written with the write_json method).

Here are the current command-line options:

::

    > lt_xspec -h
    usage: lt_xabssys [-h] [-outfile OUTFILE] [-llist LLIST] [--un_norm]
                  spec_file abssys_file

    Parse for XAbsSys

    positional arguments:
      spec_file         Spectral file
      abssys_file       AbsSys file (JSON)

    optional arguments:
      -h, --help        show this help message and exit
      -outfile OUTFILE  Output filename
      -llist LLIST      Name of LineList
      --un_norm         Spectrum is NOT normalized


From within ipython or equivalent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Not yet implemented

--------------

Navigating in the Main Window- These key strokes help you explore the spectrum (be sure to click in the spectrum panel first!)
------------------------------------------------------------------------------------------------------------------------------

Setting the window edges (mouse+keystroke)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  l -- Set left edge of window
-  r -- Set right edge of window
-  t -- Set top edge of window
-  b -- Set bottom edge of window

Zoom in/out Wavelength
~~~~~~~~~~~~~~~~~~~~~~

-  i -- Zoom in on cursor
-  I -- Zoom in extra fast
-  o -- Zoom out
-  O -- Zoom out extra fast

Zoom out Flux
~~~~~~~~~~~~~

-  Y -- Zoom out

Pan
~~~

-  [ -- Pan left
-  { -- Pan left extra
-  ] -- Pan right
-  } -- Pan right extra

Adjusting Rows/Columns
~~~~~~~~~~~~~~~~~~~~~~

-  = -- Move to next page of lines
-  - -- Move to previous page of lines
-  C -- Add a column
-  c -- Remove a column
-  K -- Add a row
-  k -- Remove a row

--------------

Modifying Absorption lines
--------------------------

Limits and Blends
~~~~~~~~~~~~~~~~~
