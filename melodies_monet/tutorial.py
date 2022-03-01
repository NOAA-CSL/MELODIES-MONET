"""
Retrieve the docs/examples tutorial datasets, caching with pooch
"""
# Inspired by xarray.tutorial


_base_url = "https://csl.noaa.gov/groups/csl4/modeldata/melodies-monet/data/"

examples = {
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
}
"""Files to fetch for a certain example, paths relative to the FTP site."""


def fetch_example(example: str):
    """Fetch files pertaining to a certain example."""
    try:
        import pooch
    except ImportError as e:
        raise RuntimeError(
            "downloading and caching the docs/examples datasets requires `pooch`"
        ) from e

    ftp_paths = examples.get(example)
    if ftp_paths is None:
        raise ValueError(f"invalid example choice. Valid options are {list(examples)}")

    for id_, (ftp_path, hash) in ftp_paths.items():
        url = f"{_base_url}{ftp_path}"
        fp = pooch.retrieve(url, known_hash=hash)
