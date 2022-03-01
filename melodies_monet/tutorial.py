"""
Retrieve the docs/examples tutorial datasets, caching with pooch
"""
# Inspired by xarray.tutorial


_base_url = "https://csl.noaa.gov/groups/csl4/modeldata/melodies-monet/data/"

_examples = {
    "wrfchem": {
        "racm_esrl": (
            "example_model_data/wrfchem_example/racm_esrl/wrfout_d01_tutorial",
            "sha256:1ac5e4a7a4311345a52f78a55ee5630d0c1b30c0df18a5ff1f745425baa31915",
        ),
        "racm_esrl_vcp": (
            "example_model_data/wrfchem_example/racm_esrl_vcp/wrfout_d01_tutorial",
            "sha256:67e0f126031dba9775a1baaeec377d04144da66d040fd27909e418c3be31a0f9",
        )
    },
    "airnow": {
        "2019-09": (
            "example_observation_data/surface/AIRNOW_20190901_20190930.nc",
            "sha256:618e7ddf5609a3f552caf0c830756bef530f75d5ca380aae932cea234f6e3753",
        ),
        "2019-08": (
            "example_observation_data/surface/AIRNOW_20190801_20190831.nc",
            "sha256:da0a864d554218e1242cb7be4b5b34e13251db6558815de990846ca0a83f618e",
        ),
        "2019-07": (
            "example_observation_data/surface/AIRNOW_20190701_20190731.nc",
            "sha256:67da806bcdca90289147254a3fa63a746b8a63ef973a497d0b47c789d48291c2",
        )
    }
}
"""Files to fetch for a certain example, paths relative to the FTP site."""

_examples_flat = {
    f"{a}:{b}": tup
    for a, dct in _examples.items()
    for b, tup in dct.items()
}

example_ids = list(_examples_flat)


def fetch_example(example: str) -> str:
    """Fetch file pertaining to a certain example model or obs dataset.

    Returns
    -------
    fp : str
        Path to the downloaded/cached example file.
    """
    try:
        import pooch
    except ImportError as e:
        raise RuntimeError(
            "downloading and caching the docs/examples datasets requires `pooch`"
        ) from e

    if example not in example_ids:
        raise ValueError(f"invalid example choice. Valid options are {example_ids}")
    else:
        ftp_path, hash = _examples_flat[example]

    url = f"{_base_url}{ftp_path}"
    fp = pooch.retrieve(url, known_hash=hash)

    return fp
