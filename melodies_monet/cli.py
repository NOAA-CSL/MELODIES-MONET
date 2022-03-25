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

    typer.echo(f"{desc} ...")
    try:
        yield
    finally:
        elapsed = time.perf_counter() - start
        typer.echo(f"{desc} took {elapsed:.3g} seconds")



def main(
    control: str = typer.Argument(
        ...,
        help="Path to the control file to use.", 
    ),
):
    """Run MELODIES MONET as described in the control file CONTROL."""

    p = Path(control)
    if not p.is_file():
        typer.echo(f"Error: control file {control!r} does not exist")
        raise typer.Exit(2)

    typer.echo(f"Using control file: {control!r}")
    typer.echo(f"with full path: {p.absolute().as_posix()}")

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
