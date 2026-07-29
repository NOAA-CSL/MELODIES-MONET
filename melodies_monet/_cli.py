# SPDX-License-Identifier: Apache-2.0
#
"""
melodies-monet -- MELODIES MONET CLI
"""
import os
import time
from contextlib import contextmanager
from pathlib import Path
from typing import List, Tuple

_LOGGING_LEVEL = os.environ.get("MM_LOGGING_LEVEL", None)
if _LOGGING_LEVEL is not None:
    import logging

    logging.basicConfig(level=_LOGGING_LEVEL.upper())

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


def _get_full_name(obj):
    """Get the full name of a function or type,
    including the module name if not builtin."""
    import builtins
    import inspect

    mod = inspect.getmodule(obj)
    name = obj.__qualname__
    if mod is None or mod is builtins:
        return name
    else:
        return f"{mod.__name__}.{name}"


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
        typer.secho(f"Error message (type: {_get_full_name(type(e))}): {e}", fg=ERROR_COLOR)
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


@contextmanager
def _ignore_pandas_numeric_only_futurewarning():
    """Disable pandas `numeric_only` FutureWarning"""
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            category=FutureWarning,
            message=(
                "The default value of numeric_only in DataFrameGroupBy.mean is deprecated. "
                "In a future version, numeric_only will default to False. "
                "Either specify numeric_only or select only columns "
                "which should be valid for the function."
            ),
        )
        yield


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
        from melodies_monet.driver import analysis
    
    with _timer("Reading control file and initializing"):
        an = analysis()
        an.control = control
        an.read_control()
        if debug and not an.debug:
            typer.secho(
                f"Setting `analysis.debug` (was {an.debug}) to True since --debug used.",
                fg=INFO_COLOR,
            )
            an.debug = True

    with _timer("Opening model(s)"):
        an.open_models()

    # Note: currently MM expects having at least model and at least one obs
    # but in the future, model-to-model only might be an option
    with _timer("Opening observations(s)"):
        an.open_obs()

    with _timer("Pairing"):
        if an.read is not None:
            an.read_analysis()
        else:
            an.pair_data()

    if an.save is not None:
        with _timer("Saving paired datasets"):
            an.save_analysis()

    if an.control_dict.get("plots") is not None:
        with _timer("Plotting and saving the figures"), _ignore_pandas_numeric_only_futurewarning():
            an.plotting()

    if an.control_dict.get("stats") is not None:
        with _timer("Computing and saving statistics"), _ignore_pandas_numeric_only_futurewarning():
            an.stats()


_DATE_FMT_NOTE = (
    "Date can be in any format accepted by `pandas.date_range()`, "
    "e.g., 'YYYY-MM-DD', or 'M/D/YYYY'. "
    "Time other than 0 UTC can be specified by adding trailing ' HH[:MM[:SS]]', "
    "but this might not have an effect on the output."
)
_DATE_END_NOTE = (
    "As not specifying time implies 0 UTC, "
    "to get the full last day for hourly data, you should specify hour, e.g., append ' 23' "
    "or increase end date by one day. "
    "For daily data, this is not necessary."
)


@app.command()
def get_aeronet(
    start_date: str = typer.Option(..., "-s", "--start-date", help=f"Start date. {_DATE_FMT_NOTE}"),
    end_date: str = typer.Option(..., "-e", "--end-date", help=f"End date. {_DATE_FMT_NOTE} {_DATE_END_NOTE}"),
    daily: bool = typer.Option(False, help="Whether to retrieve the daily averaged data product."),
    freq: str = typer.Option("h", "-f", "--freq", help=(
            "Frequency to resample to. "
            "Mean is used to reduce the time groups (as opposed to nearest, e.g.)."
        )
    ),
    interp_to: str = typer.Option(None, "--interp-to", help=(
            "Wavelength(s) to interpolate the AOD values to (unit: micron). "
            "Separate with commas to specify multiple. "
            "Examples: '0.55' (550 nm), '0.55,0.7,1.0'. "
            "Note that this functionality requires pytspack "
            "(https://github.com/noaa-oar-arl/pytspack)."
        )
    ),
    out_name: str = typer.Option(None, "-o",
        help=(
            "Output file name (or full/relative path). "
            "By default the name is generated like 'AERONET_<product>_<start-date>_<end-date>.nc'."
        )
    ),
    dst: Path = typer.Option(".", "-d", "--dst", help=(
            "Destination directory (to control output location "
            "if using default output file name)."
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

    from melodies_monet.util.write_util import write_ncf

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
            .swap_dims(siteid="x")
        )

        ds = (
            df
            .set_index(["time", "siteid"])
            .to_xarray()
            .swap_dims(siteid="x")
            .drop_vars(site_vns)
            .merge(ds_site)
            .set_coords(site_vns)
            .assign(x=range(ds_site.sizes["x"]))
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
    start_date: str = typer.Option(..., "-s", "--start-date", help=f"Start date. {_DATE_FMT_NOTE}"),
    end_date: str = typer.Option(..., "-e", "--end-date", help=f"End date. {_DATE_FMT_NOTE} {_DATE_END_NOTE}"),
    daily: bool = typer.Option(False, help=(
            "Whether to retrieve the daily averaged data product. "
            "By default, the hourly data is fetched."
        )
    ),
    out_name: str = typer.Option(None, "-o",
        help=(
            "Output file name (or full/relative path). "
            "By default the name is generated like 'AirNow_<start-date>_<end-date>.nc'."
        )
    ),
    dst: Path = typer.Option(".", "-d", "--dst", help=(
            "Destination directory (to control output location "
            "if using default output file name)."
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
    """Download AirNow data using monetio and reformat for MM usage."""
    import warnings

    import monetio as mio
    import pandas as pd

    from melodies_monet.util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    if verbose:
        from dask.diagnostics import ProgressBar

        ProgressBar().register()

    typer.echo(HEADER)

    start_date = pd.Timestamp(start_date)
    end_date = pd.Timestamp(end_date)
    dates = pd.date_range(start_date, end_date, freq="h" if not daily else "D")
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
        # NOTE: time_local not included since it varies in time as well as by site
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
            .swap_dims(siteid="x")
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
            .swap_dims(siteid="x")
            .drop_vars(site_vns)
            .merge(ds_site)
            .set_coords(["latitude", "longitude"])
            .assign(x=range(ds_site.sizes["x"]))
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


@app.command()
def get_ish_lite(
    start_date: str = typer.Option(..., "-s", "--start-date", help=f"Start date. {_DATE_FMT_NOTE}"),
    end_date: str = typer.Option(..., "-e", "--end-date", help=f"End date. {_DATE_FMT_NOTE} {_DATE_END_NOTE}"),
    country: str = typer.Option(None, "--country",
        help=(
            "Two-letter country code (e.g., in order of site count, "
            "US, RS, CA, AS, BR, IN, CH, NO, JA, UK, FR, ...)."
        )
    ),
    state: str = typer.Option(None, "--state", help="Two-letter state code (e.g., MD, ...)."),
    box: Tuple[float, float, float, float] = typer.Option((None, None, None, None), "--box",
        help=(
            "Bounding box for site selection. "
            "(latmin, lonmin, latmax, lonmax) in [-180, 180) format. "
            "Can't be used if specifying country or state."
        )
    ),
    out_name: str = typer.Option(None, "-o",
        help=(
            "Output file name (or full/relative path). "
            "By default the name is generated like 'ISH-Lite_<start-date>_<end-date>.nc'."
        )
    ),
    dst: Path = typer.Option(".", "-d", "--dst", help=(
            "Destination directory (to control output location "
            "if using default output file name)."
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
    """Download ISH-Lite data using monetio and reformat for MM usage.
    
    Note that the data are stored in yearly files by site, so the runtime
    mostly depends on the number of unique years that your date range includes,
    as well as any site selection narrowing.
    You can use --country or --state or --box to select groups of sites.
    ISH-Lite is an hourly product.
    """
    import warnings

    import monetio as mio
    import pandas as pd

    from melodies_monet.util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    if verbose:
        from dask.diagnostics import ProgressBar

        ProgressBar().register()

    typer.echo(HEADER)

    start_date = pd.Timestamp(start_date)
    end_date = pd.Timestamp(end_date)
    dates = pd.date_range(start_date, end_date, freq="h")
    if verbose:
        print("Dates:")
        print(dates)

    if box == (None, None, None, None):
        box = None

    # Set destination and file name
    fmt = r"%Y%m%d"
    if out_name is None:
        out_name = f"ISH-Lite_{start_date:{fmt}}_{end_date:{fmt}}.nc"
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
            df = mio.ish_lite.add_data(
                dates,
                box=box,
                state=state,
                country=country,
                resample=False,
                n_procs=num_workers,
                verbose=verbose,
            )

    with _timer("Computing UTC offset for selected ISH-Lite sites"):
        import datetime

        from timezonefinder import TimezoneFinder
        from pytz import timezone, utc

        tf = TimezoneFinder(in_memory=True)
        ref_date = datetime.datetime(2022, 1, 1, 0, 0)

        def get_utc_offset(*, lat, lon):
            s = tf.timezone_at(lng=lon, lat=lat)
            assert s is not None

            tz_target = timezone(s)
            ref_date_tz_target = tz_target.localize(ref_date)
            ref_date_utc = utc.localize(ref_date)
            uo_h = (ref_date_utc - ref_date_tz_target).total_seconds() / 3600

            return uo_h


        locs = df[["siteid", "latitude", "longitude"]].groupby("siteid").first().reset_index()
        locs["utcoffset"] = locs.apply(lambda r: get_utc_offset(lat=r.latitude, lon=r.longitude), axis="columns")

        df = df.merge(locs[["siteid", "utcoffset"]], on="siteid", how="left")


    with _timer("Forming xarray Dataset"):
        df = df.dropna(subset=["latitude", "longitude"])

        df = df.rename(
            columns={
                "station name": "station_name",
                "elev(m)": "elevation",
            },
            errors="ignore",
        )

        site_vns = [
            "siteid",
            "latitude",
            "longitude",
            "country",
            "state",
            "station_name",
            "usaf",
            "wban",
            "icao",
            "elevation",
            "utcoffset",
            "begin",
            "end",
        ]
        # NOTE: time_local not included since it varies in time as well as by site

        ds_site = (
            df[site_vns]
            .groupby("siteid")
            .first()
            .to_xarray()
            .swap_dims(siteid="x")
        )

        # TODO: units?
        units = {}

        cols = list(df.columns)
        ds = (
            df[cols]
            .set_index(["time", "siteid"])
            .to_xarray()
            .swap_dims(siteid="x")
            .drop_vars(site_vns)
            .merge(ds_site)
            .set_coords(["latitude", "longitude"])
            .assign(x=range(ds_site.sizes["x"]))
        )

        # Add units
        for k, u in units.items():
            vn = k
            ds[vn].attrs.update(units=u)

        # Fill in local time array
        # (in the df, not all sites have rows for all times, so we have NaTs at this point)
        ds["time_local"] = ds.time + (ds.utcoffset * 60).astype("timedelta64[m]")

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


@app.command()
def get_ish(
    start_date: str = typer.Option(..., "-s", "--start-date", help=f"Start date. {_DATE_FMT_NOTE}"),
    end_date: str = typer.Option(..., "-e", "--end-date", help=f"End date. {_DATE_FMT_NOTE} {_DATE_END_NOTE}"),
    freq: str = typer.Option("h", "-f", "--freq", help=(
            "Frequency to resample to. "
            "Mean is used to reduce the time groups (as opposed to nearest, e.g.)."
        )
    ),
    country: str = typer.Option(None, "--country",
        help=(
            "Two-letter country code (e.g., in order of site count, "
            "US, RS, CA, AS, BR, IN, CH, NO, JA, UK, FR, ...)."
        )
    ),
    state: str = typer.Option(None, "--state", help="Two-letter state code (e.g., MD, ...)."),
    box: Tuple[float, float, float, float] = typer.Option((None, None, None, None), "--box",
        help=(
            "Bounding box for site selection. "
            "(latmin, lonmin, latmax, lonmax) in [-180, 180) format. "
            "Can't be used if specifying country or state."
        )
    ),
    out_name: str = typer.Option(None, "-o",
        help=(
            "Output file name (or full/relative path). "
            "By default the name is generated like 'ISH_<start-date>_<end-date>.nc'."
        )
    ),
    dst: Path = typer.Option(".", "-d", "--dst", help=(
            "Destination directory (to control output location "
            "if using default output file name)."
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
    """Download ISH data using monetio and reformat for MM usage.
    
    Note that the data are stored in yearly files by site, so the runtime
    mostly depends on the number of unique years that your date range includes,
    as well as any site selection narrowing.
    You can use --country or --state or --box to select groups of sites.
    Time resolution may be sub-hourly, depending on site,
    thus we resample to hourly by default.
    """
    import warnings

    import monetio as mio
    import pandas as pd

    from melodies_monet.util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    if verbose:
        from dask.diagnostics import ProgressBar

        ProgressBar().register()

    typer.echo(HEADER)

    start_date = pd.Timestamp(start_date)
    end_date = pd.Timestamp(end_date)
    dates = pd.date_range(start_date, end_date, freq="h")
    if verbose:
        print("Dates:")
        print(dates)

    if box == (None, None, None, None):
        box = None

    # Set destination and file name
    fmt = r"%Y%m%d"
    if out_name is None:
        out_name = f"ISH_{start_date:{fmt}}_{end_date:{fmt}}.nc"
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

    with _timer("Fetching data with monetio"), _ignore_pandas_numeric_only_futurewarning():
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="The (error|warn)_bad_lines argument has been deprecated"
            )
            df = mio.ish.add_data(
                dates,
                box=box,
                state=state,
                country=country,
                resample=True,
                window=freq,
                n_procs=num_workers,
                verbose=verbose,
            )

    with _timer("Computing UTC offset for selected ISH sites"):
        import datetime

        from timezonefinder import TimezoneFinder
        from pytz import timezone, utc

        tf = TimezoneFinder(in_memory=True)
        ref_date = datetime.datetime(2022, 1, 1, 0, 0)

        def get_utc_offset(*, lat, lon):
            s = tf.timezone_at(lng=lon, lat=lat)
            assert s is not None

            tz_target = timezone(s)
            ref_date_tz_target = tz_target.localize(ref_date)
            ref_date_utc = utc.localize(ref_date)
            uo_h = (ref_date_utc - ref_date_tz_target).total_seconds() / 3600

            return uo_h


        locs = df[["siteid", "latitude", "longitude"]].groupby("siteid").first().reset_index()
        locs["utcoffset"] = locs.apply(lambda r: get_utc_offset(lat=r.latitude, lon=r.longitude), axis="columns")

        df = df.merge(locs[["siteid", "utcoffset"]], on="siteid", how="left")


    with _timer("Forming xarray Dataset"):
        df = (
            df.dropna(subset=["latitude", "longitude"])
            .rename(
                columns={
                    "station name": "station_name",
                    "elev(m)": "elevation",
                },
            errors="ignore",
            )
            .drop(columns=["elev"], errors="ignore")  # keep just elevation from the site meta file
        )

        site_vns = [
            "siteid",
            "latitude",
            "longitude",
            "country",
            "state",
            "station_name",
            "usaf",
            "wban",
            "icao",
            "elevation",
            "utcoffset",
            "begin",
            "end",
        ]
        # NOTE: time_local not included since it varies in time as well as by site

        ds_site = (
            df[site_vns]
            .groupby("siteid")
            .first()
            .to_xarray()
            .swap_dims(siteid="x")
        )

        # TODO: units?
        units = {}

        cols = list(df.columns)
        ds = (
            df[cols]
            .set_index(["time", "siteid"])
            .to_xarray()
            .swap_dims(siteid="x")
            .drop_vars(site_vns)
            .merge(ds_site)
            .set_coords(["latitude", "longitude"])
            .assign(x=range(ds_site.sizes["x"]))
        )

        # Add units
        for k, u in units.items():
            vn = k
            ds[vn].attrs.update(units=u)

        # Fill in local time array
        # (in the df, not all sites have rows for all times, so we have NaTs at this point)
        ds["time_local"] = ds.time + (ds.utcoffset * 60).astype("timedelta64[m]")

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


@app.command()
def get_aqs(
    start_date: str = typer.Option(..., "-s", "--start-date", help=f"Start date. {_DATE_FMT_NOTE}"),
    end_date: str = typer.Option(..., "-e", "--end-date", help=f"End date. {_DATE_FMT_NOTE} {_DATE_END_NOTE}"),
    daily: bool = typer.Option(False, help=(
            "Whether to retrieve the daily averaged data product. "
            "By default, the hourly data is fetched."
        )
    ),
    param: List[str] = typer.Option(["O3", "PM2.5", "PM10"], "-p", "--params", help=(
            "Parameter groups. "
            "Use '-p' more than once to get multiple groups. "
            "Other examples: 'SPEC' (speciated PM2.5), 'PM10SPEC' (speciated PM10), "
            "'VOC', 'NONOxNOy', 'SO2', 'NO2', 'CO', 'PM2.5_FRM'."
        )
    ),
    # TODO: add network selection option once working in monetio
    out_name: str = typer.Option(None, "-o",
        help=(
            "Output file name (or full/relative path). "
            "By default the name is generated like 'AQS_<start-date>_<end-date>.nc'."
        )
    ),
    dst: Path = typer.Option(".", "-d", "--dst", help=(
            "Destination directory (to control output location "
            "if using default output file name)."
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
    """Download EPA AQS data using monetio and reformat for MM usage.

    These are archived data, stored in per-year per-parameter-group files, described at
    https://aqs.epa.gov/aqsweb/airdata/download_files.html

    Recent-past data are generally not available from this source.
    """
    import warnings

    import monetio as mio
    import pandas as pd

    from melodies_monet.util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    if verbose:
        from dask.diagnostics import ProgressBar

        ProgressBar().register()

    typer.echo(HEADER)

    start_date = pd.Timestamp(start_date)
    end_date = pd.Timestamp(end_date)
    dates = pd.date_range(start_date, end_date, freq="h" if not daily else "D")
    if verbose:
        print("Dates:")
        print(dates)
        print("Params:")
        print(param)

    # Set destination and file name
    fmt = r"%Y%m%d"
    if out_name is None:
        out_name = f"AQS_{start_date:{fmt}}_{end_date:{fmt}}.nc"
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
            try:
                df = mio.aqs.add_data(
                    dates,
                    param=param,
                    daily=daily,
                    network=None,
                    download=False,
                    local=False,
                    wide_fmt=True,  # column for each variable
                    n_procs=num_workers,
                    meta=False,  # TODO: enable or add option once monetio fixes released
                )
            except KeyError as e:
                if daily and str(e) == "'time'":
                    typer.echo("Note that the daily option currently requires monetio >0.2.5")
                raise

    with _timer("Fetching site metadata"):
        # Need UTC offset in order to compute local time
        # But currently the `meta=True` option doesn't work
        meta0 = pd.read_csv(
            "https://aqs.epa.gov/aqsweb/airdata/aqs_sites.zip",
            encoding="ISO-8859-1",
            usecols=[0, 1, 2, 17, 24, 25],
            dtype=str,
        )
        meta = (
            meta0.copy()
            .assign(
                siteid=meta0["State Code"] + meta0["County Code"] + meta0["Site Number"],
                utcoffset=meta0["GMT Offset"].astype(int),
            )
            .drop(
                columns=["Site Number", "GMT Offset"],
            )
            .rename(
                columns={
                    "State Code": "state_code",
                    "County Code": "county_code",
                    "City Name": "city_name",
                    "CBSA Name": "cbsa_name",
                }
            )
        )
        meta.loc[meta["city_name"] == "Not in a City", "city_name"] = "Not in a city"  # normalize

        counties0 = pd.read_csv(
            "https://aqs.epa.gov/aqsweb/documents/codetables/states_and_counties.csv",
            encoding="ISO-8859-1",
            dtype=str,
        )
        counties = (
            counties0.copy()
            .rename(
                columns={
                    "State Code": "state_code",
                    "State Name": "state_name",
                    "State Abbreviation": "state_abbr",
                    "County Code": "county_code",
                    "County Name": "county_name",
                    "EPA Region": "epa_region",  # note without R prefix
                }
            )
        )
        counties["epa_region"] = "R" + counties["epa_region"].str.lstrip("0")

        meta = meta.merge(counties, on=["state_code", "county_code"], how="left")

        if daily:
            meta = meta.drop(columns=["utcoffset"])

    with _timer("Forming xarray Dataset"):
        # Select requested time period (older monetio doesn't do this)
        df = df[df.time.between(dates[0], dates[-1], inclusive="both")]

        df = df.dropna(subset=["latitude", "longitude"])

        # Variables associated with a measurement,
        # currently not properly useful in the wide format.
        if daily:
            v_vns = [
                "parameter_code",
                "poc",
                "parameter_name",
                "sample_duration",
                "pollutant_standard",
                "event_type",
                "observation_count",
                "observation_percent",
                "1st_max_value",
                "1st_max_hour",
                "aqi",
                "method_code",
                "method_name",
            ]
        else:
            v_vns = [
                "parameter_code",
                "poc",  # parameter occurrence code
                "parameter_name",
                "mdl",  # method detection limit
                "uncertainty",
                "method_type",
                "method_code",
                "method_name",
            ]
        df = df.drop(columns=v_vns).drop_duplicates()
        # TODO: may be better to get long fmt and drop these first and then pivot
        # TODO: option to average duplicate measurements at same site instead of keeping first?
        if "datum" in df:
            df = df.drop(columns=["datum"])

        site_vns = [
            "siteid",
            "state_code",
            "state_name",
            "state_abbr",
            "county_code",
            "county_name",
            "city_name",
            "cbsa_name",
            "site_num",
            "epa_region",
            "latitude",
            "longitude",
        ]
        if daily:
            site_vns.extend(["local_site_name", "address", "msa_name"])
        else:
            site_vns.append("utcoffset")
        # NOTE: time_local not included since it varies in time as well

        df = df.merge(meta, on="siteid", how="left", suffixes=(None, "_meta"))

        ds_site = (
            df[site_vns]
            .groupby("siteid")
            .first()
            .to_xarray()
            .swap_dims(siteid="x")
        )

        # Extract units info so we can add as attrs
        unit_suff = "_unit"
        unit_cols = [n for n in df.columns if n.endswith(unit_suff)]
        assert (df[unit_cols].nunique() == 1).all()
        units = df[unit_cols][~df[unit_cols].isnull()].iloc[0].to_dict()

        cols = [n for n in df.columns if not n.endswith(unit_suff)]
        ds = (
            df[cols]
            .drop(columns=[vn for vn in site_vns if vn != "siteid"])
            .drop(columns=[col for col in df.columns if col.endswith("_meta")])
            .drop_duplicates(["time", "siteid"], keep="first")
            .set_index(["time", "siteid"])
            .to_xarray()
            .swap_dims(siteid="x")
            .merge(ds_site)
            .set_coords(["latitude", "longitude"])
            .assign(x=range(ds_site.sizes["x"]))
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

        # Can't have `/` in variable name for netCDF
        to_rename = [vn for vn in ds.data_vars if "/" in vn]
        ds = ds.rename_vars({vn: vn.replace("/", "_") for vn in to_rename})

    with _timer("Writing netCDF file"):
        if compress:
            write_ncf(ds, dst / out_name, verbose=verbose)
        else:
            ds.to_netcdf(dst / out_name)


@app.command()
def get_openaq(
    start_date: str = typer.Option(..., "-s", "--start-date", help=f"Start date. {_DATE_FMT_NOTE}"),
    end_date: str = typer.Option(..., "-e", "--end-date", help=f"End date. {_DATE_FMT_NOTE} {_DATE_END_NOTE}"),
    out_name: str = typer.Option(None, "-o",
        help=(
            "Output file name (or full/relative path). "
            "By default the name is generated like 'OpenAQ_<start-date>_<end-date>.nc'."
        )
    ),
    dst: Path = typer.Option(".", "-d", "--dst", help=(
            "Destination directory (to control output location "
            "if using default output file name)."
        )
    ),
    param: List[str] = typer.Option(["o3", "pm25", "pm10"], "-p", "--param", help=(
            "Parameters. "
            "Use '-p' more than once to get multiple parameters. "
            "Other examples: 'no', 'no2', 'nox', 'so2', 'co', 'bc'. "
            "Only applicable to the web API methods ('api-v*')."
        )
    ),
    reference_grade: bool = typer.Option(True, help="Include reference-grade sensors."),
    low_cost: bool = typer.Option(False, help="Include low-cost sensors."),
    country: List[str] = typer.Option(None, "-c", "--country",
        help=(
            "Two-letter country code(s). (US, CA, MX, ...). "
            "Use more than once to specify multiple countries."
        )
    ),
    method: str = typer.Option("api-v3", "-m", "--method", help=(
            "Method (reader) to use for fetching data. "
            "Options: 'api-v3', 'api-v2', 'openaq-fetches'."
        )
    ),
    sensor_limit: int = typer.Option(None,
        help=(
            "Limit the number of sensors to fetch data for. "
            "This is useful for testing or debugging. "
            "Only applicable to the 'api-v3' method."
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
    """Download hourly OpenAQ data using monetio and reformat for MM usage."""
    import warnings

    import monetio as mio
    import pandas as pd

    from melodies_monet.util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    typer.echo(HEADER)

    if method not in {"api-v3", "api-v2", "openaq-fetches"}:
        typer.secho(f"Error: method {method!r} not recognized", fg=ERROR_COLOR)
        raise typer.Exit(2)

    start_date = pd.Timestamp(start_date)
    end_date = pd.Timestamp(end_date)

    if method in {"openaq-fetches"}:
        dates = pd.date_range(start_date, end_date, freq="D")
    elif method.startswith("api-v"):
        dates = pd.date_range(start_date, end_date, freq="h")
    else:
        raise AssertionError
    if verbose:
        print("Dates:")
        print(dates)

    if not country:
        country = None

    if verbose and method.startswith("api-v"):
        print("Params:", param)
    if verbose and method == "api-v3" and country is not None:
        print("Country:", country)
    if verbose and method == "api-v3":
        print("Sensor limit:", sensor_limit)

    # Set destination and file name
    fmt = r"%Y%m%d"
    if out_name is None:
        out_name = f"OpenAQ_{start_date:{fmt}}_{end_date:{fmt}}.nc"
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

    sensor_types = []
    if reference_grade:
        sensor_types.append("reference grade")
    if low_cost:
        sensor_types.append("low-cost sensor")
    if not sensor_types and method.startswith("api-v"):
        typer.secho(
            "Error: no sensor types selected. Use --reference-grade and/or --low-cost",
            fg=ERROR_COLOR,
        )
        raise typer.Exit(2)

    if verbose and method in {"openaq-fetches"}:
        from dask.diagnostics import ProgressBar

        ProgressBar().register()

    with _timer("Fetching data with monetio"):
        if method == "openaq-fetches":
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="The (error|warn)_bad_lines argument has been deprecated"
                )
                df = mio.openaq.add_data(
                    dates,
                    n_procs=num_workers,
                    # wide_fmt=True,
                )

            # Address time-wise non-unique site IDs
            # Some (most?) are just slightly different lat/lon
            # But seems like a few are actual time-wise lat/lon duplicates
            df = df.drop_duplicates(["time", "siteid"])

        elif method.startswith("api-v"):
            kws = dict(
                parameters=param,
                sensor_type=sensor_types,
                wide_fmt=True,
                timeout=60,
                retry=15,
                threads=num_workers if num_workers > 1 else None,
            )
            if method == "api-v3":
                kws.update(
                    hourly=True,
                    country=country,
                    sensor_limit=sensor_limit,
                )
                func = mio.obs.openaq_v3.add_data
            elif method == "api-v2":
                func = mio.obs.openaq_v2.add_data
            else:
                raise AssertionError
            df = func(
                dates,
                **kws,
            )

            dupes = df[df.duplicated(["time", "siteid"], keep=False)]
            if not dupes.empty:
                typer.echo(
                    f"warning: {len(dupes)} unexpected time-siteid duplicated rows:"
                )
                if verbose:
                    typer.echo(dupes)
                df = df.drop_duplicates(["time", "siteid"])
        else:
            raise AssertionError

        if df.empty:
            raise RuntimeError("No data found")

        if method == "api-v2":
            # Drop times not on the hour
            good = df.time == df.time.dt.floor("H")
            typer.echo(f"Dropping {(~good).sum()}/{len(good)} rows that aren't on the hour.")
            df = df[good]

    with _timer("Forming xarray Dataset"):
        df = df.drop(columns=["index"], errors="ignore")
        df = df.dropna(subset=["latitude", "longitude"])

        if method == "openaq-fetches":
            site_vns = [
                "siteid",  # based on country and lat/lon
                "latitude",
                "longitude",
                "utcoffset",
                #
                "city",
                "country",  # 2-char codes
                #
                "sourceName",
                "sourceType",  # "government"
            ]
            # NOTE: time_local not included since it varies in time as well
        elif method == "api-v2":
            site_vns = [
                "siteid",  # real OpenAQ location ID
                "latitude",
                "longitude",
                "utcoffset",
                #
                "location",
                "city",
                "country",
                #
                "entity",
                "sensor_type",
                "is_mobile",
                "is_analysis",
            ]
            for vn in ["city", "is_analysis"]:  # may have been dropped for being all null
                if vn not in df.columns:
                    site_vns.remove(vn)
        elif method == "api-v3":
            df = df.drop(
                columns=[
                    "sensor_id",
                    "period_label",
                ],
                errors="ignore",
            )
            site_vns = [
                "siteid",  # real OpenAQ location ID
                "latitude",
                "longitude",
                "utcoffset",
                #
                "country",
                #
                "is_mobile",
                "is_monitor",
            ]
        else:
            raise AssertionError

        ds_site = (
            df[site_vns]
            .groupby("siteid")
            .first()
            .to_xarray()
            .swap_dims(siteid="x")
        )

        ds = (
            df.drop(columns=[vn for vn in site_vns if vn not in ["siteid"]])
            .set_index(["time", "siteid"])
            .to_xarray()
            .swap_dims(siteid="x")
            .merge(ds_site)
            .set_coords(["latitude", "longitude"])
            .assign(x=range(ds_site.dims["x"]))
        )

        # Rename species vars and add units as attr
        nice_us = {"ppm": "ppmv", "ugm3": "ug m-3", "ppb": "pbbv"}
        for vn0 in [n for n in df.columns if n.endswith(("_ppm", "ppb", "_ugm3", "_umg3"))]:
            i_last_underscore = vn0.rfind("_")
            vn, u = vn0[:i_last_underscore], vn0[i_last_underscore + 1:]
            if u == "umg3":
                u = "ugm3"
            nice_u = nice_us[u]
            ds[vn0].attrs.update(units=nice_u)
            ds = ds.rename_vars({vn0: vn})

        # Fill in local time array
        # (in the df, not all sites have rows for all times, so we have NaTs at this point)
        ds["time_local"] = ds.time + ds.utcoffset

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
