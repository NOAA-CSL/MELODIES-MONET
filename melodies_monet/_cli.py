# Copyright (C) 2022 National Center for Atmospheric Research and National Oceanic and Atmospheric Administration
# SPDX-License-Identifier: Apache-2.0
#
"""
melodies-monet -- MELODIES MONET CLI
"""
import time
from contextlib import contextmanager
from pathlib import Path

try:
    import typer
except ImportError as e:
    print(
        "The MELODIES MONET CLI requires the module 'typer'. "
        "You can install it with `conda install -c conda-forge typer` or "
        "`pip install typer`. "
        f"The error message was: {e}"
    )
    raise SystemExit(1)

DEBUG = False
INFO_COLOR = typer.colors.CYAN
ERROR_COLOR = typer.colors.BRIGHT_RED
SUCCESS_COLOR = typer.colors.GREEN

HEADER = """
------------------
| MELODIES MONET |
------------------    
""".strip()


@contextmanager
def _timer(desc=""):
    start = time.perf_counter()

    tpl = f"{desc} {{status}} in {{elapsed:.3g}} seconds"

    typer.secho(f"{desc} ...", fg=INFO_COLOR)
    try:
        yield
    except Exception as e:
        typer.secho(
            tpl.format(status="failed", elapsed=time.perf_counter() - start),
            fg=ERROR_COLOR
        )
        typer.secho(f"Error message: {e}", fg=ERROR_COLOR)
        if DEBUG:
            raise
        else:
            typer.echo("(Use the '--debug' flag to see more info.)")
            raise typer.Exit(1)
    else:
        typer.secho(
            tpl.format(status="succeeded", elapsed=time.perf_counter() - start),
            fg=SUCCESS_COLOR
        )


def _version_callback(value: bool):
    from . import __version__

    if value:
        typer.echo(f"melodies-monet {__version__}")
        # TODO: monet/monetio versions?
        raise typer.Exit()


app = typer.Typer()


@app.callback()
def main(
    version: bool = typer.Option(
        False, "--version/", help="Print version.", callback=_version_callback, is_eager=True
    ),
):
    """MELODIES MONET"""


@app.command()
def run(
    control: str = typer.Argument(
        ...,
        help="Path to the control file to use.", 
    ),
    debug: bool = typer.Option(
        False, "--debug/", help="Print more messages (including full tracebacks)."
    ),
):
    """Run MELODIES MONET as described in the control file CONTROL."""

    global DEBUG

    DEBUG = debug

    p = Path(control)
    if not p.is_file():
        typer.echo(f"Error: control file {control!r} does not exist")
        raise typer.Exit(2)

    typer.echo(HEADER)
    typer.secho(f"Using control file: {control!r}", fg=INFO_COLOR)
    typer.secho(f"with full path: {p.absolute().as_posix()}", fg=INFO_COLOR)

    with _timer("Importing the driver"):
        from .driver import analysis
    
    with _timer("Reading control file and initializing"):
        an = analysis()
        an.control = control
        an.read_control()

    with _timer("Opening model(s)"):
        an.open_models()

    # Note: currently MM expects having at least model and at least one obs
    # but in the future, model-to-model only might be an option
    with _timer("Opening observations(s)"):
        an.open_obs()

    with _timer("Pairing"):
        an.pair_data()

    if an.control_dict.get("plots") is not None:
        with _timer("Plotting and saving the figures"):
            an.plotting()

    if an.control_dict.get("stats") is not None:
        with _timer("Computing and saving statistics"):
            an.stats()


@app.command()
def get_aeronet(
    start_date: str = typer.Option(..., "-s", "--start-date", help="Start date."),
    end_date: str = typer.Option(..., "-e", "--end-date", help="End date."),
    daily: bool = typer.Option(False, help="Whether to retrieve the daily averaged data product."),
    freq: str = typer.Option("H", "-f", "--freq", help="Frequency to resample to."),
    interp_to: str = typer.Option(None, "--interp-to", help=(
            "Wavelength(s) to interpolate the AOD values to (micron). "
            "Separate with commas to specify multiple. "
            "Note that this functionality requires pytspack."
        )
    ),
    out_name: str = typer.Option(None, "-o",
        help=(
            "Output file name (or full/relative path). "
            "By default the name is generated like 'AERONET_<product>_<start-date>_<end-date>.nc'"
        )
    ),
    dst: Path = typer.Option(".", "-d", "--dst", help=(
            "Destination directory (to control output location "
            "if using default `out_name`)."
        )
    ),
    compress: bool = typer.Option(True, help=(
            "If true, pack float to int and apply compression using zlib with complevel 7. "
            "This can take time if the dataset is large, but can lead to "
            "significant space savings."
        )
    ),
    num_workers: int = typer.Option(1, "-n", "--num-workers", help="Number of download workers."),
    verbose: bool = typer.Option(False),
    debug: bool = typer.Option(
        False, "--debug/", help="Print more messages (including full tracebacks)."
    ),
):
    """Download AERONET data using monetio and reformat for MM usage."""
    import monetio as mio
    import numpy as np
    import pandas as pd

    from .util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    typer.echo(HEADER)

    start_date = pd.Timestamp(start_date)
    end_date = pd.Timestamp(end_date)
    dates = pd.date_range(start_date, end_date, freq="D")

    # Set destination and file name
    fmt = r"%Y%m%d"
    if out_name is None:
        out_name = f"AERONET_L15_{start_date:{fmt}}_{end_date:{fmt}}.nc"
    else:
        p = Path(out_name)
        if p.name == out_name:
            # `out_name` is just the file name
            out_name = p.name
        else:
            # `out_name` has path
            if dst != Path("."):
                typer.echo(f"warning: overriding `dst` setting {dst.as_posix()!r} with `out_name` {p.as_posix()!r}")
            dst = p.parent
            out_name = p.name

    if interp_to is not None:
        interp_to = np.array([float(x.strip()) for x in interp_to.strip().split(",")])
        interp_to *= 1000  # um -> nm

    with _timer("Fetching data with monetio"):
        try:
            df = mio.aeronet.add_data(
                dates,
                interp_to_aod_values=interp_to,
                daily=daily,
                freq=freq,
                n_procs=num_workers,
                verbose=1 if verbose else 0,
            )
        except ValueError:
            if daily and interp_to is not None:
                typer.echo("Note that using interp with the daily product requires monetio >0.2.2")
            raise
  
    site_vns = [
        "siteid",
        "latitude",
        "longitude",
        "aeronet_instrument_number",
        "elevation",
    ]

    with _timer("Forming xarray Dataset"):
        df = df.dropna(subset=["latitude", "longitude"])

        # Site-specific variables should only vary in x.
        # Here we take the first non-NaN value (should all be same).
        ds_site = (
            df[site_vns]
            .groupby("siteid")
            .first()  # TODO: would be nice to confirm unique-ness
            .to_xarray()
            .rename_dims(siteid="x")
        )

        ds = (
            df
            .set_index(["time", "siteid"])
            .to_xarray()
            .rename_dims(siteid="x")
            .drop_vars(site_vns)
            .merge(ds_site)
            .set_coords(site_vns)
            .assign(x=range(ds_site.dims["x"]))
            .expand_dims("y")
            .transpose("time", "y", "x")
        )

    with _timer("Writing netCDF file"):
        if compress:
            write_ncf(ds, dst / out_name, verbose=verbose)
        else:
            ds.to_netcdf(dst / out_name)


@app.command()
def get_airnow(
    start_date: str = typer.Option(..., "-s", "--start-date", help="Start date."),
    end_date: str = typer.Option(..., "-e", "--end-date", help="End date."),
    daily: bool = typer.Option(False, help=(
            "Whether to retrieve the daily averaged data product. "
            "By default, the hourly data is fetched."
        )
    ),
    out_name: str = typer.Option(None, "-o",
        help=(
            "Output file name (or full/relative path). "
            "By default the name is generated like 'AirNow_<start-date>_<end-date>.nc'"
        )
    ),
    dst: Path = typer.Option(".", "-d", "--dst", help=(
            "Destination directory (to control output location "
            "if using default `out_name`)."
        )
    ),
    compress: bool = typer.Option(True, help=(
            "If true, pack float to int and apply compression using zlib with complevel 7. "
            "This can take time if the dataset is large, but can lead to "
            "significant space savings."
        )
    ),
    num_workers: int = typer.Option(1, "-n", "--num-workers", help="Number of download workers."),
    verbose: bool = typer.Option(False),
    debug: bool = typer.Option(
        False, "--debug/", help="Print more messages (including full tracebacks)."
    ),
):
    """Download AirNow data using monetio and reformat for MM usage.
    
    The date range used is closed on both sides.
    """
    import warnings

    import monetio as mio
    import pandas as pd

    from .util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    typer.echo(HEADER)

    start_date = pd.Timestamp(start_date)
    end_date = pd.Timestamp(end_date)
    dates = pd.date_range(start_date, end_date, freq="H" if not daily else "D")
    if verbose:
        print("Dates:")
        print(dates)

    # Set destination and file name
    fmt = r"%Y%m%d"
    if out_name is None:
        out_name = f"AirNow_{start_date:{fmt}}_{end_date:{fmt}}.nc"
    else:
        p = Path(out_name)
        if p.name == out_name:
            # `out_name` is just the file name
            out_name = p.name
        else:
            # `out_name` has path
            if dst != Path("."):
                typer.echo(f"warning: overriding `dst` setting {dst.as_posix()!r} with `out_name` {p.as_posix()!r}")
            dst = p.parent
            out_name = p.name

    with _timer("Fetching data with monetio"):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="The (error|warn)_bad_lines argument has been deprecated"
            )
            df = mio.airnow.add_data(
                dates,
                download=False,
                wide_fmt=True,  # column for each variable
                n_procs=num_workers,
                daily=daily,
            )

    with _timer("Forming xarray Dataset"):
        df = df.dropna(subset=["latitude", "longitude"])

        site_vns = [
            "site",
            "siteid",
            "utcoffset",
            "latitude",
            "longitude",
            "cmsa_name",
            "msa_code",
            "msa_name",
            "state_name",
            "epa_region",
        ]
        # NOTE: time_local not included since it varies in time as well
        if daily:
            site_vns.remove("utcoffset")  # not present in the daily data product

        # site_vn_str = [
        #     "site",  # site name
        #     "siteid",  # site code (9 or 12 digits/chars)
        #     #
        #     "cmsa_name",
        #     "msa_code",
        #     "msa_name",
        #     "state_name",
        #     "epa_region",
        # ]

        # df[site_vn_str] = df[site_vn_str].astype("string")

        ds_site = (
            df[site_vns]
            # .replace(["", " ", None], pd.NA)  # TODO: monetio should do?
            .groupby("siteid")
            .first()
            .to_xarray()
            .rename_dims(siteid="x")
        )

        # Extract units info so we can add as attrs
        unit_suff = "_unit"
        unit_cols = [n for n in df.columns if n.endswith(unit_suff)]
        assert (df[unit_cols].nunique() == 1).all()
        units = df[unit_cols][~df[unit_cols].isnull()].iloc[0].to_dict()

        cols = [n for n in df.columns if not n.endswith(unit_suff)]
        ds = (
            df[cols]
            .set_index(["time", "siteid"])
            .to_xarray()
            .rename_dims(siteid="x")
            .drop_vars(site_vns)
            .merge(ds_site)
            .set_coords(["latitude", "longitude"])
            .assign(x=range(ds_site.dims["x"]))
        )

        # Add units
        for k, u in units.items():
            vn = k[:-len(unit_suff)]
            ds[vn].attrs.update(units=u)

        # Fill in local time array
        # (in the df, not all sites have rows for all times, so we have NaTs at this point)
        if not daily:
            ds["time_local"] = ds.time + ds.utcoffset.astype("timedelta64[h]")

        # Expand
        ds = (
            ds
            .expand_dims("y")
            .transpose("time", "y", "x")
        )

    with _timer("Writing netCDF file"):
        if compress:
            write_ncf(ds, dst / out_name, verbose=verbose)
        else:
            ds.to_netcdf(dst / out_name)


cli = app

_typer_click_object = typer.main.get_command(app)  # for sphinx-click in docs


if __name__ == "__main__":
    cli()
