# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pytest

from melodies_monet.plots import savefig

mpl.use("Agg")

def test_savefig(tmpdir):
    save_dir = Path(tmpdir)

    fig = plt.figure()

    fp = save_dir / "asdf.png"
    
    # Currently must be str, not Path
    with pytest.raises(AttributeError, match="has no attribute 'split'"):
        savefig(fp)

    savefig(fp.as_posix())

    plt.close(fig)
