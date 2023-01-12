# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
import os
import pytest

from melodies_monet.util import analysis_util


def test_find_file(tmpdir):

    print(tmpdir)

    test_file = os.path.join(tmpdir, 'test.txt')
    f = open(test_file, 'w')
    f.close()

    filename = analysis_util.find_file(tmpdir, 'test*')
    print(filename)
