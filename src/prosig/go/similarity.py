"""Molecular Function GO semantic similarity using Lin score."""

from __future__ import annotations

import pickle
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from prosig.go.build import MF_NAMESPACE


@dataclass(frozen=True)
class GoTermRecord:
    """Display and scoring metadata for one GO term."""

    go_id: str
    name: str
    freq: float | None
    ic: float | None
    depth: int | None
    parents: tuple[str, ...]
    children: tuple[str, ...]


@dataclass(frozen=True)
class GoLinSimilarityResult:
    """Detailed Lin similarity result for two GO terms."""

    go1: str
    go2: str
    similarity: float | None
    mica: str | None
    ic_go1: float | None
    ic_go2: float | None
    ic_mica: float | None
    status: str
    reason: str
    common_ancestors: tuple[str, ...] = ()


class GoSimilarity:
    """Reusable MF-only GO similarity index over a ProSig GO graph artifact."""

    def __init__(self, artifact: dict[str, Any]) -> None:
        meta = artifact.get("meta", {})
        namespace = meta.get("namespace")
        if namespace is not None and namespace != MF_NAMESPACE:
            raise ValueError(f"GO artifact namespace must be {MF_NAMESPACE!r}")
        terms = artifact.get("terms")
        if not isinstance(terms, dict):
            raise ValueError("GO artifact is missing a 'terms' mapping")
        self.meta = meta
        self.terms = terms

    @classmethod
    def from_pickle(cls, path: str | Path) -> GoSimilarity:
        """Load a ProSig GO graph pickle and return a similarity index."""
        with Path(path).open("rb") as handle:
            artifact = pickle.load(handle)
        if not isinstance(artifact, dict):
            raise ValueError("GO graph pickle must contain a dictionary artifact")
        return cls(artifact)

    def term(self, go_id: str) -> GoTermRecord | None:
        """Return metadata for one GO term, or None if absent."""
        term = self.terms.get(go_id)
        if term is None:
            return None
        return GoTermRecord(
            go_id=go_id,
            name=str(term.get("name", "")),
            freq=_optional_float(term.get("freq")),
            ic=_optional_float(term.get("ic")),
            depth=_optional_int(term.get("depth")),
            parents=tuple(term.get("parents", ())),
            children=tuple(term.get("children", ())),
        )

    def ancestors_including_self(self, go_id: str) -> frozenset[str]:
        """Return precomputed ancestors for a term, including the term itself."""
        term = self.terms.get(go_id)
        if term is None:
            return frozenset()
        return frozenset((go_id, *term.get("ancestors", ())))

    def find_mica(self, go1: str, go2: str) -> str | None:
        """Return the most informative common ancestor for two GO terms."""
        common = self.ancestors_including_self(go1) & self.ancestors_including_self(go2)
        candidates = [
            (go_id, ic)
            for go_id in common
            if (ic := self._ic(go_id)) is not None
        ]
        if not candidates:
            return None
        return max(candidates, key=lambda item: (item[1], item[0]))[0]

    def lin(self, go1: str, go2: str) -> float | None:
        """Return the Lin similarity score without constructing diagnostics."""
        term1 = self.terms.get(go1)
        term2 = self.terms.get(go2)
        if term1 is None or term2 is None:
            return None

        ic1 = _optional_float(term1.get("ic"))
        ic2 = _optional_float(term2.get("ic"))
        if ic1 is None or ic2 is None:
            return None

        denominator = ic1 + ic2
        if denominator == 0.0:
            return None

        mica_ic = self._mica_ic_from_raw_terms(go1, term1, go2, term2)
        if mica_ic is None:
            return None
        return 2 * mica_ic / denominator

    def lin_with_details(self, go1: str, go2: str) -> GoLinSimilarityResult:
        """Return Lin score plus diagnostics for two MF GO terms."""
        term1 = self.term(go1)
        term2 = self.term(go2)
        if term1 is None and term2 is None:
            return _unavailable(go1, go2, "missing_go1_go2")
        if term1 is None:
            return _unavailable(go1, go2, "missing_go1", term2=term2)
        if term2 is None:
            return _unavailable(go1, go2, "missing_go2", term1=term1)
        if term1.ic is None:
            return _unavailable(go1, go2, "missing_ic_go1", term1=term1, term2=term2)
        if term2.ic is None:
            return _unavailable(go1, go2, "missing_ic_go2", term1=term1, term2=term2)

        common = self.ancestors_including_self(go1) & self.ancestors_including_self(go2)
        common_with_ic = tuple(
            sorted(
                (go_id for go_id in common if self._ic(go_id) is not None),
                key=lambda go_id: (self._ic(go_id) or 0.0, go_id),
            )
        )
        if not common_with_ic:
            return _unavailable(
                go1,
                go2,
                "no_common_ancestor",
                term1=term1,
                term2=term2,
            )

        mica = max(common_with_ic, key=lambda go_id: (self._ic(go_id) or 0.0, go_id))
        ic_mica = self._ic(mica)
        denominator = term1.ic + term2.ic
        if denominator == 0.0:
            return _unavailable(
                go1,
                go2,
                "zero_ic_denominator",
                term1=term1,
                term2=term2,
                mica=mica,
                ic_mica=ic_mica,
                common_ancestors=common_with_ic,
            )

        return GoLinSimilarityResult(
            go1=go1,
            go2=go2,
            similarity=2 * (ic_mica or 0.0) / denominator,
            mica=mica,
            ic_go1=term1.ic,
            ic_go2=term2.ic,
            ic_mica=ic_mica,
            status="ok",
            reason="",
            common_ancestors=common_with_ic,
        )

    def _ic(self, go_id: str) -> float | None:
        term = self.terms.get(go_id)
        if term is None:
            return None
        return _optional_float(term.get("ic"))

    def _mica_ic_from_raw_terms(
        self,
        go1: str,
        term1: dict[str, Any],
        go2: str,
        term2: dict[str, Any],
    ) -> float | None:
        ancestors1 = term1.get("ancestors", ())
        ancestors2 = term2.get("ancestors", ())
        if len(ancestors1) <= len(ancestors2):
            return self._max_common_ancestor_ic(
                go1,
                ancestors1,
                go2,
                ancestors2,
            )
        return self._max_common_ancestor_ic(
            go2,
            ancestors2,
            go1,
            ancestors1,
        )

    def _max_common_ancestor_ic(
        self,
        iter_go: str,
        iter_ancestors: Any,
        other_go: str,
        other_ancestors: Any,
    ) -> float | None:
        mica: tuple[float, str] | None = None
        for ancestor in _iter_ancestors_including_self(iter_go, iter_ancestors):
            if ancestor != other_go and ancestor not in other_ancestors:
                continue
            ic = self._ic(ancestor)
            if ic is None:
                continue
            candidate = (ic, ancestor)
            if mica is None or candidate > mica:
                mica = candidate
        if mica is None:
            return None
        return mica[0]


def _optional_float(value: object) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _optional_int(value: object) -> int | None:
    if isinstance(value, int):
        return value
    return None


def _iter_ancestors_including_self(
    go_id: str,
    ancestors: Any,
) -> Iterator[str]:
    yield go_id
    yield from ancestors


def _unavailable(
    go1: str,
    go2: str,
    reason: str,
    *,
    term1: GoTermRecord | None = None,
    term2: GoTermRecord | None = None,
    mica: str | None = None,
    ic_mica: float | None = None,
    common_ancestors: tuple[str, ...] = (),
) -> GoLinSimilarityResult:
    return GoLinSimilarityResult(
        go1=go1,
        go2=go2,
        similarity=None,
        mica=mica,
        ic_go1=term1.ic if term1 is not None else None,
        ic_go2=term2.ic if term2 is not None else None,
        ic_mica=ic_mica,
        status="unavailable",
        reason=reason,
        common_ancestors=common_ancestors,
    )
