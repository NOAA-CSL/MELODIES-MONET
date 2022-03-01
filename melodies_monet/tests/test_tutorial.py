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
    ...
