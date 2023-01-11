# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
import pytest

from pathlib import Path

from melodies_monet import util


def test_find_file(tmpdir):

    save_dir = Path(tmpdir)
