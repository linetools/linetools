# Module to run tests on spectra.utils

from __future__ import print_function, absolute_import, \
     division, unicode_literals

import pytest
from astropy.time import Time

from linetools.spectra import utils as lsu


def test_get_COS_LP_from_date():
    dates_good = ["2009-12-11", "2012-12-12", "2014-06-06",
                  "2017-12-12", Time("2017-12-12")]
    dates_bad =["2008-10-10", "2012-07-23","2014-02-09", "2017-10-02", 1, None]

    for date in dates_good:
        lsu.get_COS_LP_from_date(date)

    for date in dates_bad:
        with pytest.raises(ValueError):
            lsu.get_COS_LP_from_date(date)

