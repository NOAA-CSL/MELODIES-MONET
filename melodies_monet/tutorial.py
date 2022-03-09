"""
Retrieve the docs/examples tutorial datasets, caching with pooch

The data are retrieved from
https://csl.noaa.gov/groups/csl4/modeldata/melodies-monet/
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
        ),
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
        ),
    },
    "camchem": {
        "fv": (
            "example_model_data/cesmfv_example/CAM_chem_merra2_FCSD_1deg_QFED_world_201909-01-09_small_sfc.nc",
            "sha256:b1efe4475796c5a507466c84a101efeb18ff142ddc40a21d0c4691527e7378c4",
        ),
        "se": (
            "example_model_data/musica_example/Sample_MUSICAv0_CONUS_2019-09-05.nc",
            "sha256:0bf2437bdc8c9a6f686bb8f8f1c38f5f7f877d996d31f7c15bbd51f6c4c100db",
        ),
        "se_scrip": (
            "example_model_data/musica_example/ne0CONUS_ne30x8_np4_SCRIP.nc",
            "sha256:890e1e98f52a5687c57cb15e52f481aa17c6b9eea2bf8e860cae2301697cc027",
        ),
    },
    "aeronet": {
        "2019-09": (
            "example_observation_data/surface/AERONET_L15_20190901_20190930.nc",
            "4b32e679ce3d5698ff6aecae4488082dc4bd6f96dd6e02be972f4e0f8a101523",
        ),
        "2019-08": (
            "example_observation_data/surface/AERONET_L15_20190801_20190831.nc",
            "0b29761de0e38b542d6da0a0b52c6131254fa8a5ecb0981a4393eb7ded2eef36",
        ),
    },
    "csn": {
        "2019_daily": (
            "example_observation_data/surface/CSN_DAILY_2019.nc",
            "b430a56da37e69bd2adbb5f17918f9cb14583c652433b43c42c7627bd06edf55",
        ),
    },
    "improve": {
        "2019_daily" : (
            "example_observation_data/surface/IMPROVE_DAILY_2019.nc",
            "599afc11238c30345f9f756e379d527b2e67b5674cdcf5074df1a93677fe7f9d",
        )
    },
    "ncore": {
        "2019_daily" : (
            "example_observation_data/surface/NCORE_DAILY_2019.nc",
            "888fb70f7f6cd9af8b49398b56240fcb70ebe886152c143e6a6016074d4f0bfe",
        )
    },
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
