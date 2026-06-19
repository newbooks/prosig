import logging

LOG_FORMAT = "[%(levelname)s]: %(message)s"
LOGGER_NAME = "prosig"


def configure_logging(level: str) -> None:
    numeric_level = getattr(logging, level.upper(), None)
    if not isinstance(numeric_level, int):
        valid_levels = ", ".join(log_level_names())
        raise ValueError(f"invalid log level {level!r}; choose one of: {valid_levels}")

    logging.basicConfig(
        format=LOG_FORMAT,
        level=numeric_level,
        force=True,
    )


def get_logger() -> logging.Logger:
    return logging.getLogger(LOGGER_NAME)


def log_level_names() -> list[str]:
    return ["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
