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
    assert "build-library" in result.stdout
    assert "inspect" in result.stdout


def test_build_library_help_includes_options() -> None:
    result = CliRunner().invoke(app, ["build-library", "-h"])

    assert result.exit_code == 0
    assert "Build the compact GO graph" in result.stdout
    assert "--go-obo" in result.stdout
    assert "--swissprot" in result.stdout
    assert "--go-out" in result.stdout
    assert "--prosite-dat" in result.stdout
    assert "--motif-out" in result.stdout
    assert "--write-report" in result.stdout
    assert "--role-map" in result.stdout
    assert "--cluster-out" in result.stdout
    assert "--cluster-config" in result.stdout
    assert "--cluster-neighbors" not in result.stdout
    assert "--cluster-resolution" not in result.stdout
    assert "--cluster-stats-out" not in result.stdout
    assert "--cluster-progress-interval" not in result.stdout
    assert "--cluster-term-cache" not in result.stdout
    assert "--cluster-profile-cache" not in result.stdout
    assert "--cluster-min-informative-ic" not in result.stdout
    assert "--cluster-max-posting-fraction" not in result.stdout
    assert "--cluster-max-posting-size" not in result.stdout
    assert "--force" in result.stdout
    assert "-f" in result.stdout
    assert "--namespace" not in result.stdout
    assert "--include-part-of" not in result.stdout
    assert "--ic-log-base" not in result.stdout
    assert "--min-count" not in result.stdout


def test_log_level_option_suppresses_info_logs() -> None:
    result = CliRunner().invoke(
        app, ["--log-level", "WARNING", "setup-data", "--dry-run"]
    )

    assert result.exit_code == 0
    assert "[INFO]:" not in result.output
