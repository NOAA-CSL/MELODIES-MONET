import pytest

from melodies_monet import tutorial

try:
    import pooch
except ImportError:
    has_pooch = False
else:
    has_pooch = True


@pytest.mark.xfail(not has_pooch, reason="no pooch")
def test_fetch_example():
    cache_dir = pooch.os_cache("pooch")

    tutorial.fetch_example("wrfchem:racm_esrl")
    assert any("wrfout_d01_tutorial" in p.name for p in cache_dir.glob("*"))

    tutorial.fetch_example("airnow:2019-09")
    assert any("AIRNOW_20190901_20190930.nc" in p.name for p in cache_dir.glob("*"))

    with pytest.raises(ValueError, match="invalid example choice"):
        tutorial.fetch_example("asdf")
