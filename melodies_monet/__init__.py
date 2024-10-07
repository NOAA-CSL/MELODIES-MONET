# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
"""
MELODIES MONET
"""
import sys

__version__ = "0.1"

_submodules = [
    "driver",
    "plots",
    "stats",
    "util",
    "tutorial",
]

__all__ = [__version__] + _submodules


if sys.version_info < (3, 7):
    from . import driver, plots, stats, util, tutorial

else:
    # Lazy imports
    # https://peps.python.org/pep-0562/
    # https://github.com/scipy/scipy/blob/fdc31f3d2bfa90a2c214a398668f0e153632d2bb/scipy/__init__.py#L181-L195
    def __dir__():
        return __all__

    import importlib as _importlib

    def __getattr__(name):
        if name in _submodules:
            return _importlib.import_module(f"melodies_monet.{name}")
        else:
            try:
                return globals()[name]
            except KeyError:
                raise AttributeError(f"Module 'melodies_monet' has no attribute '{name}'")
