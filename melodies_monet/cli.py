"""
melodies-monet -- MELODIES MONET CLI
"""
from pathlib import Path
from typing import Optional

import typer


def main(
    control: Optional[str] = typer.Argument(
        None,
        help="Path to control file to run.", 
        show_default="./control.yaml"
    ),
):
    """Run MELODIES MONET as described in the control file CONTROL.
    If not provided, 'control.yaml' is looked for in the current directory.    
    """
    if control is None:
        control = "./control.yaml"

    if not Path(control).is_file():
        typer.echo(f"error: {control!r} does not exist")
        return 1

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
