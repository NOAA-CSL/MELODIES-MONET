# SPDX-License-Identifier: Apache-2.0
#
import warnings

import numpy as np
import pytest
import xarray as xr

from melodies_monet.util.write_util import compress_variable


def _da(values):
    return xr.DataArray(np.asarray(values, dtype="float64"), dims="x", name="v")


def test_compress_variable_warns_inexact_fill():
    # Min below -1, so the filled -1 doesn't land on a packed integer
    da = _da(np.linspace(-10, 100, 11))
    da[3] = np.nan

    with pytest.warns(UserWarning, match=r"Variable 'v': NaN values filled with -1 will decode as"):
        compress_variable(da)


def test_compress_variable_no_warning_exact_fill():
    # Min is the fill value itself, which packs to the lowest integer and decodes exactly
    da = _da(np.linspace(0, 100, 11))
    da[3] = np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        compress_variable(da)


def test_compress_variable_no_warning_without_nan():
    da = _da(np.linspace(-10, 100, 11))

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        compress_variable(da)
