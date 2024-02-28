# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
import os
import pytest
from datetime import datetime

from melodies_monet.util import analysis_util


def test_fill_date_template():

    date = datetime.now()
    date_str = date.strftime('%Y-%m-%b-%d-%j')
    print(date_str)

    template_str = 'Year YYYY, Month MM, Month Name M_ABBR, Day DD'
    filled_str = analysis_util.fill_date_template(template_str, date_str)
    print(filled_str)
    assert(filled_str == date.strftime('Year %Y, Month %m, Month Name %b, Day %d'))

    template_str = 'Year YYYY, Julian Day DDD'
    filled_str = analysis_util.fill_date_template(template_str, date_str)
    print(filled_str)
    assert(filled_str == date.strftime('Year %Y, Julian Day %j'))


def test_find_file(tmpdir):

    test_file = os.path.join(tmpdir, 'test.txt')
    f = open(test_file, 'w')
    f.close()

    filename = analysis_util.find_file(tmpdir, 'test*')
    print(filename)
    assert(filename == test_file)
