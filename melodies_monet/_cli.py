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
    freq: str = typer.Option(None, "-f", "--freq", help="Frequency to resample to."),
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
    num_workers: int = typer.Option(1, "-n", "--num-workers", help="Number of download workers."),
    verbose: bool = typer.Option(False),
):
    """Download AERONET data using monetio and reformat for MM usage."""
    import monetio as mio
    import numpy as np
    import pandas as pd

    from .util.write_util import write_ncf

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

    if interp_to is not None:
        interp_to = np.array([float(x.strip()) for x in interp_to.strip().split(",")])
        interp_to *= 1000  # um -> nm

    try:
        df = mio.aeronet.add_data(
            dates,
            interp_to_aod_values=interp_to,
            daily=daily,
            freq=freq,
            n_procs=num_workers,
            verbose=1 if verbose else 0,
        )
    except ValueError as e:
        print(f"Error loading AERONET: {e}")
        if daily and interp_to is not None:
            print("Note that using interp with the daily product requires monetio >0.2.2")
        raise typer.Exit(1)
    except Exception as e:
        print(f"Error loading AERONET: {e}")
        raise typer.Exit(1)
  
    site_vns = [
        "siteid",
        "latitude",
        "longitude",
        "aeronet_instrument_number",
        "elevation",
    ]

    # Site-specific variables should only vary in x
    ds_site = (
        df[site_vns]
        .groupby("siteid")
        .first()  # TODO: would be nice to confirm unique-ness
        .to_xarray()
        .rename_dims(siteid="x")
    )

    ds = (
        df
        .dropna(subset=["latitude", "longitude"])
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

    write_ncf(ds, dst / out_name, verbose=verbose)


cli = app

_typer_click_object = typer.main.get_command(app)  # for sphinx-click in docs


if __name__ == "__main__":
    cli()
