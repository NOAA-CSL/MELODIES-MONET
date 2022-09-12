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


@app.command(name="run")
def main(
    control: str = typer.Argument(
        ...,
        help="Path to the control file to use.", 
    ),
    debug: bool = typer.Option(
        False, "--debug/", help="Print more messages (including full tracebacks)."
    ),
    version: bool = typer.Option(
        False, "--version/", help="Print version.", callback=_version_callback, is_eager=True
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
    freq: str = typer.Option(None, "-f", "--freq", help="Frequency to resample to."),
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
    num_workers: int = typer.Option(1, "-n", "--num-workers", help="Number of download workers."),
    verbose: bool = typer.Option(False),
):
    """Download AERONET data using monetio and reformat for MM usage."""
    import monetio as mio
    import numpy as np
    import pandas as pd
    import xarray as xr

    from .util.write_util import write_ncf

    def expand_dims(ds, index=0, site_variable=None):
        from numpy import unique

        # First set a new index for the siteid
        ds["x"] = index
        ds = ds.expand_dims(["x"])
        ds = ds.set_coords(["x"])

        # Now reduce the site variables to single element variables
        for sv in site_variable:
            tmp = [unique(ds[sv])[0]]
            ds[sv] = (("x",), tmp)

        return ds

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
                print(f"warning: overriding `dst` setting {dst.as_posix()!r} with `out_name` {p.as_posix()!r}")
            dst = p.parent
            out_name = p.name

    #standard_wavelengths = np.array([0.34, 0.44, 0.55, 0.66, 0.86, 1.63, 11.1]) * 1000.0
    # ^ some of these overlap with existing wls in the dataset
    #   so the later `dfp.to_xarray()` fails since we get duplicate column names
    standard_wavelengths = np.array([0.55]) * 1000.0
    # ^ only really need 550 nm for the comparison to UFS-Aerosol
    if daily:
        standard_wavelengths = None
        # TODO: currently doesn't work with the daily data due to differences
        # in column names
    df = mio.aeronet.add_data(
        dates,
        interp_to_aod_values=standard_wavelengths,
        daily=daily,
        freq=freq,
        n_procs=num_workers,
        verbose=1 if verbose else 0,
    )

    dfp = df.rename({"siteid": "x"}, axis=1).set_index(["time", "x"])
    columns = dfp.columns.to_list()
    good_columns = []
    remove_columns = []
    for column in columns:
        good_columns.append(column)
        try:
            dfp[good_columns].to_xarray()
            if verbose:
                print("COLUMN SUCCESS:", column)
        except:
            if verbose:
                print("COLUMN FAILURE:", column)
                remove_columns.append(column)
                good_columns.pop()

    dfp = dfp.drop(columns=remove_columns).dropna(subset=["latitude", "longitude"])
    dft = df.drop(columns=remove_columns)
    dsets = [
        dft.loc[df.siteid == s].set_index(["time"]).to_xarray()
        for s in df.siteid.unique()
    ]

    site_variable = [
        "siteid",
        "latitude",
        "longitude",
        "aeronet_instrument_number",
        "elevation",
    ]

    # Now `site_variable` are the single element attributes of the individual sites
    # so lets simplify this
    for index, d in enumerate(dsets):
        dsets[index] = expand_dims(d, index=index, site_variable=site_variable)

    # Now combine all the datasets for each site into a single dataset
    ds = xr.concat(dsets, dim="x").set_coords(site_variable)
    
    # Write the file
    t = ds.expand_dims("y").transpose("time", "y", "x")
    write_ncf(t, dst / out_name, verbose=verbose)


cli = app

_typer_click_object = typer.main.get_command(app)  # for sphinx-click in docs


if __name__ == "__main__":
    cli()
