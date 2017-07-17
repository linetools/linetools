.. _AbsLine:

*************
AbsLine Class
*************

.. index:: AbsLine

Notebooks
=========

.. toctree::
   :maxdepth: 1

   Simple Examples <AbsLine_examples>


Overview
========

This Class is a child of the abstract
:ref:`SpectralLine` Class.  See that
documentation for the base methods.

AbsLine is designed to organize and analyze an absorption line.
In addition to the attributes defaulted to SpectralLine,
this class has:

================ ================= ========= ========================================
Property         Variable          Type      Description
================ ================= ========= ========================================
Doppler param.   attrib['b']       Quantity  Doppler parameter
b sigma          attrib['sig_b']   Quantity  1 sigma uncertainty in b
Column density   attrib['N']       Quantity  Column density
N sigma          attrib['sig_N']   Quantity  1 sigma uncertainty in N
N flag           attrib['flag_N']  int       Column density flag
================ ================= ========= ========================================


Instantiation
=============

from_dict
---------

Instantiate from a dict.  The keys *ltype* ('Abs')
and *trans* are required.

fill_data
---------

Attributes
==========

See the Table above.  logN and sig_logN are commonly used
as well.

Plots
=====

Methods
=======

get_Wr_from_N_b
---------------
It returns the rest-frame equivalent width for a given `N` and `b`. It uses the approximation given
by Draine 2011 book (eq. 9.27), which comes from atomic physics considerations.::

    abslin1 = AbsLine('HI 1215')
    N = [10**12.0, 10**12.1, 10**12.2, 10**14.0] / (u.cm*u.cm)
    b = [20, 20, 20, 20] * u.km / u.s
    Wr = abslin1.get_Wr_from_N_b(N, b)
    print(Wr)
    <Quantity [ 0.0053758 , 0.00674464, 0.00845469, 0.21149773] Angstrom>

get_Wr_from_N
---------------
It returns the approximated rest-frame equivalent width for
a given N. It uses the approximation given by Draine 2011 book
(eq. 9.15), which is valid for tau0<<1 where Wr is independent
of Doppler parameter or gamma.::

    abslin1 = AbsLine('HI 1215')
    N = [10**12.0, 10**12.1, 10**12.2, 10**14.0] / (u.cm*u.cm)
    Wr = abslin1.get_Wr_from_N(N)
    print(Wr)
    <Quantity [ 0.00544783, 0.00685842, 0.00863423, 0.5447833 ] Angstrom>

We can see how the first 4 `Wr` estimations are good approximation to the more exact solution given by
get_Wr_from_N_b(), however, the fourth is off by a factor `>2` because the approximation tau0<<1 is not satisfied.

get_N_from_Wr
-------------
It returns the approximated column density N, for a given rest-frame equivalent width
Wr. This is an approximation only valid for tau0 << 1, where
Wr is independent on Doppler parameter and gamma (see eqs. 9.14 and 9.15 of
Draine 2011). This may be useful to put upper limits on non-detections.::

    abslin1 = AbsLine('HI 1215')
    Wr = [ 0.00544783, 0.00685842, 0.00863423, 0.5447833 ] * u.AA
    N = abslin1.get_N_from_Wr(Wr)
    print(np.log10(N.value))
    array([ 11.99999976,  12.10000029,  12.19999983,  14.        ])


generate_voigt
--------------

measure_aodm
------------


Output
======
