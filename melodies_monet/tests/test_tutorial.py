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

    def _check_in_cache(s):
        # TODO: option to check exact number?
        assert any(s in p.name for p in cache_dir.glob("*"))

    # TODO: make these a parametrize test
    tutorial.fetch_example("wrfchem:racm_esrl")
    _check_in_cache("wrfout_d01_tutorial")

    tutorial.fetch_example("airnow:2019-09")
    _check_in_cache("AIRNOW_20190901_20190930.nc")

    tutorial.fetch_example("aeronet:2019-09")
    _check_in_cache("AERONET_L15_20190901_20190930.nc")

    tutorial.fetch_example("csn:2019_daily")
    _check_in_cache("CSN_DAILY_2019.nc")

    tutorial.fetch_example("improve:2019_daily")
    _check_in_cache("IMPROVE_DAILY_2019.nc")

    with pytest.raises(ValueError, match="invalid example choice"):
        tutorial.fetch_example("asdf")
