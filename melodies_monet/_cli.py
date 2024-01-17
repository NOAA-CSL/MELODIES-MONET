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

from typing import Tuple

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
        from .driver import analysis
    
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
        an.pair_data()

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
    freq: str = typer.Option("H", "-f", "--freq", help=(
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

    from .util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    if verbose:
        from dask.diagnostics import ProgressBar

        ProgressBar().register()

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
            "[latmin, lonmin, latmax, lonmax] in [-180, 180) format. "
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

    from .util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    if verbose:
        from dask.diagnostics import ProgressBar

        ProgressBar().register()

    typer.echo(HEADER)

    start_date = pd.Timestamp(start_date)
    end_date = pd.Timestamp(end_date)
    dates = pd.date_range(start_date, end_date, freq="H")
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
            .assign(x=range(ds_site.dims["x"]))
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
    freq: str = typer.Option("H", "-f", "--freq", help=(
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
    You can use --country or --state to select groups of sites.
    Time resolution may be sub-hourly, depending on site,
    thus we resample to hourly by default.
    """
    import warnings

    import monetio as mio
    import pandas as pd

    from .util.write_util import write_ncf

    global DEBUG

    DEBUG = debug

    if verbose:
        from dask.diagnostics import ProgressBar

        ProgressBar().register()

    typer.echo(HEADER)

    start_date = pd.Timestamp(start_date)
    end_date = pd.Timestamp(end_date)
    dates = pd.date_range(start_date, end_date, freq="H")
    if verbose:
        print("Dates:")
        print(dates)

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
            .assign(x=range(ds_site.dims["x"]))
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



cli = app

_typer_click_object = typer.main.get_command(app)  # for sphinx-click in docs


if __name__ == "__main__":
    cli()
