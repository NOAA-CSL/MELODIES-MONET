"""
Retrieve the docs/examples tutorial datasets, caching with pooch
"""
# Inspired by xarray.tutorial


_base_url = "https://csl.noaa.gov/groups/csl4/modeldata/melodies-monet/data/"

examples = {
    "wrfchem": (
        "example_model_data/wrfchem_example/racm_esrl/wrfout_d01_tutorial",
        "example_model_data/wrfchem_example/racm_esrl_vcp/wrfout_d01_tutorial",
    )
}
"""Files to fetch for a certain example, paths relative to the FTP site."""
# TODO: MD5 hashes


def fetch_example(example: str):
    try:
        import pooch
    except ImportError as e:
        raise RuntimeError(
            "downloading and caching the docs/examples datasets requires `pooch`"
        ) from e

    ftp_paths = examples.get(example)
    if ftp_paths is None:
        raise ValueError(f"invalid example choice. Valid options are {list(examples)}")

    for ftp_path in ftp_paths:
        url = f"{_base_url}{ftp_path}"
        fp = pooch.retrieve(url, known_hash=None)

