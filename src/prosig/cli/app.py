from importlib.metadata import PackageNotFoundError, version as package_version

import typer

APP_NAME = "prosig"
DEVELOPER_NAME = "Junjun Mao"
DEVELOPER_EMAIL = "junjun.mao@gmail.com"

app = typer.Typer(
    help="Protein signature discovery and function inference.",
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
)


def get_version() -> str:
    try:
        return package_version(APP_NAME)
    except PackageNotFoundError:
        return "0+unknown"


@app.callback()
def cli() -> None:
    """Protein signature discovery and function inference."""


@app.command()
def version() -> None:
    """Show the installed ProSig version and developer information."""
    typer.echo(f"ProSig version: {get_version()}")
    typer.echo(f"Developer: {DEVELOPER_NAME} <{DEVELOPER_EMAIL}>")


def main() -> None:
    app()
