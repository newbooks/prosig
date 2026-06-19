from typer.testing import CliRunner

from prosig.cli.app import app


def test_version_command() -> None:
    result = CliRunner().invoke(app, ["version"])

    assert result.exit_code == 0
    assert "ProSig version:" in result.stdout
    assert "Developer: Junjun Mao <junjun.mao@gmail.com>" in result.stdout


def test_short_help_option() -> None:
    result = CliRunner().invoke(app, ["-h"])

    assert result.exit_code == 0
    assert "Usage:" in result.stdout
    assert "version" in result.stdout


def test_log_level_option_suppresses_info_logs() -> None:
    result = CliRunner().invoke(app, ["--log-level", "WARNING", "fetch", "--dry-run"])

    assert result.exit_code == 0
    assert "[INFO]:" not in result.output
