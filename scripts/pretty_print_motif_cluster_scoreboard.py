#!/usr/bin/env python3
"""Pretty-print a motif-cluster scoreboard pickle artifact."""

from __future__ import annotations

import argparse
import pickle
import pprint
from pathlib import Path
from typing import Any

def main() -> None:
    args = _parse_args()
    scoreboard = _load_pickle(args.scoreboard)
    pprint.pp(scoreboard, sort_dicts=True, width=args.width)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pretty-print a motif-cluster scoreboard pickle artifact."
    )
    parser.add_argument(
        "scoreboard",
        type=Path,
        help="scoreboard pickle to print",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=100,
        help="maximum output width in characters (default: 100)",
    )
    args = parser.parse_args()
    if args.width < 1:
        parser.error("--width must be at least 1")
    return args


def _load_pickle(path: Path) -> Any:
    if not path.is_file():
        raise SystemExit(f"Scoreboard pickle does not exist: {path}")
    with path.open("rb") as handle:
        return pickle.load(handle)


if __name__ == "__main__":
    main()
