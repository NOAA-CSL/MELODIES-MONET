"""
melodies-monet -- MELODIES MONET CLI
"""
import time
from contextlib import contextmanager
from pathlib import Path

# import io
# import threading
# from itertools import cycle

import typer

DEBUG = False
INFO_COLOR = typer.colors.CYAN
ERROR_COLOR = typer.colors.BRIGHT_RED
SUCCESS_COLOR = typer.colors.GREEN


# class Spinner(object):
#     # Based on https://gist.github.com/cevaris/79700649f0543584009e
#     spinner_cycle = cycle(['-', '/', '|', '\\'])

#     def __init__(self, desc="", *, speed=0.25):
#         self.stop_running = threading.Event()
#         self.spin_thread = threading.Thread(target=self.init_spin)
#         self._desc = desc
#         self._speed = speed

#     def start(self):
#         self.spin_thread.start()

#     def stop(self):
#         self.stop_running.set()
#         self.spin_thread.join()

#     def init_spin(self):
#         import sys

#         print(f"{self._desc} ", end="")
#         while not self.stop_running.is_set():
#             print(next(self.spinner_cycle), flush=True, end="")
#             time.sleep(self._speed)
#             print("\b", end="")


# @contextmanager
# def _spinner(desc=""):
#     start = time.perf_counter()

#     print(f"{desc} ...")
#     try:
#         spinner = Spinner(desc, speed=0.2)
#         spinner.start()
#         yield
#     finally:
#         spinner.stop()
#         elapsed = time.perf_counter() - start


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


def main(
    control: str = typer.Argument(
        ...,
        help="Path to the control file to use.", 
    ),
    debug: bool = typer.Option(
        False,
        "--debug/",
        help="Print more messages (including full tracebacks).",
    )
):
    """Run MELODIES MONET as described in the control file CONTROL."""

    global DEBUG

    DEBUG = debug

    p = Path(control)
    if not p.is_file():
        typer.echo(f"Error: control file {control!r} does not exist")
        raise typer.Exit(2)

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


def cli():
    typer.run(main)


if __name__ == "__main__":
    cli()
