"""
melodies-monet -- MELODIES MONET CLI
"""
from pathlib import Path

import typer


def main(
    control: str = typer.Argument(
        ...,
        help="Path to the control file to use.", 
    ),
):
    """Run MELODIES MONET as described in the control file CONTROL."""

    if not Path(control).is_file():
        typer.echo(f"Error: control file {control!r} does not exist")
        raise typer.Exit(2)

    from .driver import analysis
    
    an = analysis()
    an.control = control
    an.read_control()

    an.open_models()

    an.open_obs()

    an.pair_data()

    an.plotting()

    an.stats()


def cli():
    typer.run(main)


if __name__ == "__main__":
    cli()
