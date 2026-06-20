from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version

import typer

from prosig.cli.build_library import build_library
from prosig.cli.logging import configure_logging, log_level_names
from prosig.cli.setup_data import setup_data

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
def cli(
    log_level: str = typer.Option(
        "INFO",
        "--log-level",
        case_sensitive=False,
        help="Set the log level.",
    ),
) -> None:
    """Protein signature discovery and function inference."""
    if log_level.upper() not in log_level_names():
        valid_levels = ", ".join(log_level_names())
        raise typer.BadParameter(f"choose one of: {valid_levels}")
    configure_logging(log_level)


@app.command()
def version() -> None:
    """Show the installed ProSig version and developer information."""
    typer.echo(f"ProSig version: {get_version()}")
    typer.echo(f"Developer: {DEVELOPER_NAME} <{DEVELOPER_EMAIL}>")


app.command(name="setup-data")(setup_data)
app.command(name="build-library")(build_library)


def main() -> None:
    app()
