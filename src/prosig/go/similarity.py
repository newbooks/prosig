"""Molecular Function GO semantic similarity using Lin score."""

from __future__ import annotations

import logging
import pickle
import re
from collections import OrderedDict
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, TypeAlias

import numpy as np

from prosig.go.build import MF_NAMESPACE

GO_ID_PATTERN = re.compile(r"GO:\d{7}")
GO_ID_TOKEN_PATTERN = re.compile(r"^GO:\d{7}$")
TERM_PAIR_CACHE_ENTRY_BYTES = 256
PROFILE_PAIR_CACHE_ENTRY_BYTES = 512


@dataclass(frozen=True)
class CacheStats:
    """Bounded cache accounting for long-running clustering diagnostics."""

    budget_mb: int
    max_entries: int
    entries: int
    hits: int
    misses: int
    evictions: int


class BoundedSimilarityCache:
    """Approximate memory-bounded LRU cache for scalar similarity scores."""

    def __init__(self, max_bytes: int, *, entry_bytes: int) -> None:
        self.max_entries = max(0, max_bytes // entry_bytes)
        self._scores: OrderedDict[Any, float | None] = OrderedDict()
        self.hits = 0
        self.misses = 0
        self.evictions = 0

    def __len__(self) -> int:
        return len(self._scores)

    def get(self, key: Any) -> tuple[bool, float | None]:
        """Return (found, score), marking found entries recently used."""
        try:
            score = self._scores.pop(key)
        except KeyError:
            self.misses += 1
            return False, None
        self._scores[key] = score
        self.hits += 1
        return True, score

    def put(self, key: Any, score: float | None) -> None:
        """Store a score, evicting least-recently-used entries past the budget."""
        if self.max_entries == 0:
            return
        if key in self._scores:
            self._scores.pop(key)
        self._scores[key] = score
        while len(self._scores) > self.max_entries:
            self._scores.popitem(last=False)
            self.evictions += 1

    def stats(self, *, budget_mb: int) -> CacheStats:
        """Return current cache accounting."""
        return CacheStats(
            budget_mb=budget_mb,
            max_entries=self.max_entries,
            entries=len(self._scores),
            hits=self.hits,
            misses=self.misses,
            evictions=self.evictions,
        )


class BoundedTermPairCache(BoundedSimilarityCache):
    """Bounded cache for scalar GO term-pair Lin scores."""

    def __init__(self, max_bytes: int) -> None:
        super().__init__(max_bytes, entry_bytes=TERM_PAIR_CACHE_ENTRY_BYTES)


class BoundedProfilePairCache(BoundedSimilarityCache):
    """Bounded cache for scalar GO profile-pair AMB scores."""

    def __init__(self, max_bytes: int) -> None:
        super().__init__(max_bytes, entry_bytes=PROFILE_PAIR_CACHE_ENTRY_BYTES)


TermPairCache: TypeAlias = dict[int, float | None] | BoundedTermPairCache
ProfilePairCache: TypeAlias = (
    dict[tuple[tuple[int, ...], tuple[int, ...]], float | None]
    | BoundedProfilePairCache
)


@dataclass(frozen=True)
class FastGoSimilarityIndex:
    """Runtime index for repeated scalar GO similarity calculations."""

    ic_by_term: dict[str, float]
    ancestors_by_term: dict[str, frozenset[str]]
    term_id_by_term: dict[str, int]
    term_by_id: tuple[str, ...]
    ic_by_id: tuple[float, ...]
    ancestor_mask_by_id: tuple[int, ...]
    ancestors_by_desc_ic_by_id: tuple[tuple[int, ...], ...]


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


@dataclass(frozen=True)
class GoBestMatch:
    """One directional GO set best-match row."""

    source: str
    target: str
    score: float


@dataclass(frozen=True)
class GoSetSimilarityResult:
    """Detailed AMB similarity result for two GO term sets."""

    query1: str
    query2: str
    terms1: tuple[str, ...]
    terms2: tuple[str, ...]
    valid_terms1: tuple[str, ...]
    valid_terms2: tuple[str, ...]
    similarity: float | None
    status: str
    reason: str
    best_matches_1_to_2: tuple[GoBestMatch, ...] = ()
    best_matches_2_to_1: tuple[GoBestMatch, ...] = ()
    missing_terms1: tuple[str, ...] = ()
    missing_terms2: tuple[str, ...] = ()


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

    def set_lin_amb(
        self,
        terms1: tuple[str, ...] | list[str] | set[str],
        terms2: tuple[str, ...] | list[str] | set[str],
    ) -> float | None:
        """Return scalar AMB similarity over pairwise Lin scores."""
        valid_terms1 = tuple(
            term for term in _deduplicate(terms1) if term in self.terms
        )
        valid_terms2 = tuple(
            term for term in _deduplicate(terms2) if term in self.terms
        )
        return self._set_lin_amb_for_valid_terms(valid_terms1, valid_terms2)

    def set_lin_amb_with_details(
        self,
        terms1: tuple[str, ...] | list[str] | set[str],
        terms2: tuple[str, ...] | list[str] | set[str],
        *,
        query1: str = "",
        query2: str = "",
    ) -> GoSetSimilarityResult:
        """Return AMB score plus directional best-match diagnostics."""
        all_terms1 = _deduplicate(terms1)
        all_terms2 = _deduplicate(terms2)
        valid_terms1 = tuple(term for term in all_terms1 if term in self.terms)
        valid_terms2 = tuple(term for term in all_terms2 if term in self.terms)
        missing_terms1 = tuple(term for term in all_terms1 if term not in self.terms)
        missing_terms2 = tuple(term for term in all_terms2 if term not in self.terms)

        if not valid_terms1 or not valid_terms2:
            return GoSetSimilarityResult(
                query1=query1,
                query2=query2,
                terms1=all_terms1,
                terms2=all_terms2,
                valid_terms1=valid_terms1,
                valid_terms2=valid_terms2,
                similarity=None,
                status="unavailable",
                reason="empty_cleaned_set",
                missing_terms1=missing_terms1,
                missing_terms2=missing_terms2,
            )

        best_1_to_2, best_2_to_1 = self._directional_best_matches(
            valid_terms1,
            valid_terms2,
        )
        if not best_1_to_2 or not best_2_to_1:
            return GoSetSimilarityResult(
                query1=query1,
                query2=query2,
                terms1=all_terms1,
                terms2=all_terms2,
                valid_terms1=valid_terms1,
                valid_terms2=valid_terms2,
                similarity=None,
                status="unavailable",
                reason="no_valid_pairwise_similarity",
                best_matches_1_to_2=best_1_to_2,
                best_matches_2_to_1=best_2_to_1,
                missing_terms1=missing_terms1,
                missing_terms2=missing_terms2,
            )

        return GoSetSimilarityResult(
            query1=query1,
            query2=query2,
            terms1=all_terms1,
            terms2=all_terms2,
            valid_terms1=valid_terms1,
            valid_terms2=valid_terms2,
            similarity=_amb_from_best_matches(best_1_to_2, best_2_to_1),
            status="ok",
            reason="",
            best_matches_1_to_2=best_1_to_2,
            best_matches_2_to_1=best_2_to_1,
            missing_terms1=missing_terms1,
            missing_terms2=missing_terms2,
        )

    def _ic(self, go_id: str) -> float | None:
        term = self.terms.get(go_id)
        if term is None:
            return None
        return _optional_float(term.get("ic"))

    def _set_lin_amb_for_valid_terms(
        self,
        valid_terms1: tuple[str, ...],
        valid_terms2: tuple[str, ...],
    ) -> float | None:
        if not valid_terms1 or not valid_terms2:
            return None
        best_1_to_2, best_2_to_1 = self._directional_best_matches(
            valid_terms1,
            valid_terms2,
        )
        if not best_1_to_2 or not best_2_to_1:
            return None
        return _amb_from_best_matches(best_1_to_2, best_2_to_1)

    def _directional_best_matches(
        self,
        valid_terms1: tuple[str, ...],
        valid_terms2: tuple[str, ...],
    ) -> tuple[tuple[GoBestMatch, ...], tuple[GoBestMatch, ...]]:
        best_1_to_2_by_term: dict[str, GoBestMatch] = {}
        best_2_to_1_by_term: dict[str, GoBestMatch] = {}
        for term1 in valid_terms1:
            for term2 in valid_terms2:
                score = self.lin(term1, term2)
                if score is None:
                    continue
                current_1_to_2 = best_1_to_2_by_term.get(term1)
                if current_1_to_2 is None or score > current_1_to_2.score:
                    best_1_to_2_by_term[term1] = GoBestMatch(term1, term2, score)
                current_2_to_1 = best_2_to_1_by_term.get(term2)
                if current_2_to_1 is None or score > current_2_to_1.score:
                    best_2_to_1_by_term[term2] = GoBestMatch(term2, term1, score)

        return (
            tuple(
                best_1_to_2_by_term[term]
                for term in valid_terms1
                if term in best_1_to_2_by_term
            ),
            tuple(
                best_2_to_1_by_term[term]
                for term in valid_terms2
                if term in best_2_to_1_by_term
            ),
        )

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


def build_fast_go_similarity_index(artifact: dict[str, Any]) -> FastGoSimilarityIndex:
    """Build a compact runtime index for repeated scalar GO similarity calls."""
    terms = artifact.get("terms")
    if not isinstance(terms, dict):
        raise ValueError("GO artifact is missing a 'terms' mapping")

    ic_by_term = {
        go_id: float(ic)
        for go_id, term in terms.items()
        if isinstance(term, dict) and isinstance((ic := term.get("ic")), (int, float))
    }
    term_by_id = tuple(sorted(ic_by_term))
    term_id_by_term = {
        go_id: term_id for term_id, go_id in enumerate(term_by_id)
    }
    ic_by_id = tuple(ic_by_term[go_id] for go_id in term_by_id)
    ancestors_by_term = {
        go_id: frozenset(
            ancestor
            for ancestor in (go_id, *terms.get(go_id, {}).get("ancestors", ()))
            if ancestor in ic_by_term
        )
        for go_id in term_by_id
    }

    ancestor_mask_by_id: list[int] = []
    ancestors_by_desc_ic_by_id: list[tuple[int, ...]] = []
    for go_id in term_by_id:
        ancestor_ids = [
            term_id_by_term[ancestor]
            for ancestor in ancestors_by_term[go_id]
            if ancestor in term_id_by_term
        ]
        mask = 0
        for ancestor_id in ancestor_ids:
            mask |= 1 << ancestor_id
        ancestor_mask_by_id.append(mask)
        ancestors_by_desc_ic_by_id.append(
            tuple(
                sorted(
                    ancestor_ids,
                    key=lambda ancestor_id: (ic_by_id[ancestor_id], ancestor_id),
                    reverse=True,
                )
            )
        )

    return FastGoSimilarityIndex(
        ic_by_term=ic_by_term,
        ancestors_by_term=ancestors_by_term,
        term_id_by_term=term_id_by_term,
        term_by_id=term_by_id,
        ic_by_id=ic_by_id,
        ancestor_mask_by_id=tuple(ancestor_mask_by_id),
        ancestors_by_desc_ic_by_id=tuple(ancestors_by_desc_ic_by_id),
    )


def valid_go_profile(
    index: FastGoSimilarityIndex,
    terms: tuple[str, ...] | list[str] | set[str],
) -> tuple[str, ...]:
    """Return sorted deduplicated terms present in a fast similarity index."""
    return tuple(
        sorted(term for term in _deduplicate(terms) if term in index.ic_by_term)
    )


def lin_fast(
    index: FastGoSimilarityIndex,
    go1: str,
    go2: str,
    *,
    term_pair_cache: TermPairCache | None = None,
) -> float | None:
    """Compute scalar Lin similarity using the fast index."""
    term_id1 = index.term_id_by_term.get(go1)
    term_id2 = index.term_id_by_term.get(go2)
    if term_id1 is None or term_id2 is None:
        return None

    cache_key = _canonical_term_pair_id(term_id1, term_id2, len(index.ic_by_id))
    if isinstance(term_pair_cache, BoundedTermPairCache):
        found, score = term_pair_cache.get(cache_key)
        if found:
            return score
    elif isinstance(term_pair_cache, dict) and cache_key in term_pair_cache:
        return term_pair_cache[cache_key]

    ic1 = index.ic_by_id[term_id1]
    ic2 = index.ic_by_id[term_id2]
    denominator = ic1 + ic2
    if denominator == 0.0:
        similarity = None
    else:
        ancestor_mask1 = index.ancestor_mask_by_id[term_id1]
        ancestor_mask2 = index.ancestor_mask_by_id[term_id2]
        if not ancestor_mask1 & ancestor_mask2:
            similarity = None
        else:
            ancestors = index.ancestors_by_desc_ic_by_id[term_id1]
            other_ancestor_mask = ancestor_mask2
            if len(index.ancestors_by_desc_ic_by_id[term_id2]) < len(ancestors):
                ancestors = index.ancestors_by_desc_ic_by_id[term_id2]
                other_ancestor_mask = ancestor_mask1
            mica_ic = 0.0
            for ancestor_id in ancestors:
                if other_ancestor_mask & (1 << ancestor_id):
                    mica_ic = index.ic_by_id[ancestor_id]
                    break
            similarity = 2 * mica_ic / denominator

    if isinstance(term_pair_cache, BoundedTermPairCache):
        term_pair_cache.put(cache_key, similarity)
    elif isinstance(term_pair_cache, dict):
        term_pair_cache[cache_key] = similarity
    return similarity


def build_lin_similarity_matrix(
    index: FastGoSimilarityIndex,
    *,
    logger: logging.Logger | None = None,
    progress_interval_seconds: float = 60.0,
) -> np.ndarray:
    """Build a dense float32 pairwise Lin matrix from a fast GO index."""
    import time

    term_count = len(index.term_by_id)
    matrix = np.full((term_count, term_count), np.nan, dtype=np.float32)
    total_pairs = term_count * (term_count + 1) // 2
    computed_pairs = 0
    last_log_time = time.monotonic()
    if logger is not None:
        logger.info(
            "Building GO term Lin similarity matrix for %s terms "
            "(%s upper-triangle pairs)",
            f"{term_count:,}",
            f"{total_pairs:,}",
        )

    for term_id1 in range(term_count):
        for term_id2 in range(term_id1, term_count):
            similarity = lin_fast_by_id(index, term_id1, term_id2)
            if similarity is not None:
                matrix[term_id1, term_id2] = similarity
                matrix[term_id2, term_id1] = similarity
            computed_pairs += 1
        if (
            logger is not None
            and time.monotonic() - last_log_time >= progress_interval_seconds
        ):
            last_log_time = time.monotonic()
            logger.info(
                "Built GO term Lin matrix rows %s/%s; computed %s/%s pairs",
                f"{term_id1 + 1:,}",
                f"{term_count:,}",
                f"{computed_pairs:,}",
                f"{total_pairs:,}",
            )

    if logger is not None:
        logger.info(
            "Built GO term Lin similarity matrix for %s terms (%s pairs)",
            f"{term_count:,}",
            f"{total_pairs:,}",
        )
    return matrix


def lin_fast_by_id(
    index: FastGoSimilarityIndex,
    term_id1: int,
    term_id2: int,
) -> float | None:
    """Compute scalar Lin similarity using integer term IDs."""
    ic1 = index.ic_by_id[term_id1]
    ic2 = index.ic_by_id[term_id2]
    denominator = ic1 + ic2
    if denominator == 0.0:
        return None

    ancestor_mask1 = index.ancestor_mask_by_id[term_id1]
    ancestor_mask2 = index.ancestor_mask_by_id[term_id2]
    if not ancestor_mask1 & ancestor_mask2:
        return None

    ancestors = index.ancestors_by_desc_ic_by_id[term_id1]
    other_ancestor_mask = ancestor_mask2
    if len(index.ancestors_by_desc_ic_by_id[term_id2]) < len(ancestors):
        ancestors = index.ancestors_by_desc_ic_by_id[term_id2]
        other_ancestor_mask = ancestor_mask1
    mica_ic = 0.0
    for ancestor_id in ancestors:
        if other_ancestor_mask & (1 << ancestor_id):
            mica_ic = index.ic_by_id[ancestor_id]
            break
    return 2 * mica_ic / denominator


def set_lin_amb_fast_for_valid_profiles(
    index: FastGoSimilarityIndex,
    profile1: tuple[str, ...],
    profile2: tuple[str, ...],
    *,
    term_pair_cache: TermPairCache | None = None,
    profile_pair_cache: ProfilePairCache | None = None,
    lin_similarity_matrix: np.ndarray | None = None,
) -> float | None:
    """Compute scalar AMB Lin similarity for prevalidated sorted GO profiles."""
    profile_ids1 = _profile_ids(index, profile1)
    profile_ids2 = _profile_ids(index, profile2)
    cache_key = _canonical_profile_pair_id(profile_ids1, profile_ids2)
    if isinstance(profile_pair_cache, BoundedProfilePairCache):
        found, score = profile_pair_cache.get(cache_key)
        if found:
            return score
    elif isinstance(profile_pair_cache, dict) and cache_key in profile_pair_cache:
        return profile_pair_cache[cache_key]

    if not profile1 or not profile2:
        similarity = None
    else:
        best_1_to_2: dict[str, float] = {}
        best_2_to_1: dict[str, float] = {}
        for term1 in profile1:
            term_id1 = index.term_id_by_term[term1]
            for term2 in profile2:
                term_id2 = index.term_id_by_term[term2]
                if lin_similarity_matrix is None:
                    score = lin_fast(
                        index,
                        term1,
                        term2,
                        term_pair_cache=term_pair_cache,
                    )
                else:
                    score_value = lin_similarity_matrix[term_id1, term_id2]
                    score = None if np.isnan(score_value) else float(score_value)
                if score is None:
                    continue
                if term1 not in best_1_to_2 or score > best_1_to_2[term1]:
                    best_1_to_2[term1] = score
                if term2 not in best_2_to_1 or score > best_2_to_1[term2]:
                    best_2_to_1[term2] = score

        if not best_1_to_2 or not best_2_to_1:
            similarity = None
        else:
            mean_1_to_2 = sum(best_1_to_2.values()) / len(best_1_to_2)
            mean_2_to_1 = sum(best_2_to_1.values()) / len(best_2_to_1)
            similarity = (mean_1_to_2 + mean_2_to_1) / 2

    if isinstance(profile_pair_cache, BoundedProfilePairCache):
        profile_pair_cache.put(cache_key, similarity)
    elif isinstance(profile_pair_cache, dict):
        profile_pair_cache[cache_key] = similarity
    return similarity


def set_lin_amb_fast(
    index: FastGoSimilarityIndex,
    terms1: tuple[str, ...] | list[str] | set[str],
    terms2: tuple[str, ...] | list[str] | set[str],
    *,
    term_pair_cache: TermPairCache | None = None,
    profile_pair_cache: ProfilePairCache | None = None,
    lin_similarity_matrix: np.ndarray | None = None,
) -> float | None:
    """Compute scalar AMB Lin similarity between two GO-term sets."""
    return set_lin_amb_fast_for_valid_profiles(
        index,
        valid_go_profile(index, terms1),
        valid_go_profile(index, terms2),
        term_pair_cache=term_pair_cache,
        profile_pair_cache=profile_pair_cache,
        lin_similarity_matrix=lin_similarity_matrix,
    )


def _canonical_term_pair_id(term_id_a: int, term_id_b: int, term_count: int) -> int:
    lower, upper = (
        (term_id_a, term_id_b)
        if term_id_a <= term_id_b
        else (term_id_b, term_id_a)
    )
    return lower * term_count + upper


def _profile_ids(
    index: FastGoSimilarityIndex,
    profile: tuple[str, ...],
) -> tuple[int, ...]:
    return tuple(index.term_id_by_term[term] for term in profile)


def _canonical_profile_pair_id(
    profile_ids_a: tuple[int, ...],
    profile_ids_b: tuple[int, ...],
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    return (
        (profile_ids_a, profile_ids_b)
        if profile_ids_a <= profile_ids_b
        else (profile_ids_b, profile_ids_a)
    )


def _optional_float(value: object) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    return None


def _optional_int(value: object) -> int | None:
    if isinstance(value, int):
        return value
    return None


def _deduplicate(terms: tuple[str, ...] | list[str] | set[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(terms))


def _amb_from_best_matches(
    best_1_to_2: tuple[GoBestMatch, ...],
    best_2_to_1: tuple[GoBestMatch, ...],
) -> float:
    mean_1_to_2 = sum(match.score for match in best_1_to_2) / len(best_1_to_2)
    mean_2_to_1 = sum(match.score for match in best_2_to_1) / len(best_2_to_1)
    return (mean_1_to_2 + mean_2_to_1) / 2


def parse_go_term_set(value: str) -> tuple[str, ...]:
    """Parse a GO term set, with optional parentheses around separators."""
    query = value.strip()
    has_open = query.startswith("(")
    has_close = query.endswith(")")
    if has_open != has_close:
        raise ValueError(f"GO term set has unbalanced parentheses: {value}")
    if has_open and has_close:
        query = query[1:-1].strip()
    body = query
    if not body:
        raise ValueError(f"No GO terms found in set input: {value}")
    terms = [term.strip() for term in re.split(r"[;,]", body) if term.strip()]
    invalid_terms = [term for term in terms if not GO_ID_TOKEN_PATTERN.fullmatch(term)]
    if invalid_terms:
        invalid = ", ".join(invalid_terms)
        raise ValueError(f"Malformed GO term(s) in set input: {invalid}")
    return _deduplicate(terms)


def is_go_term_set_input(value: str) -> bool:
    query = value.strip()
    return (
        query.startswith("(")
        or query.endswith(")")
        or ";" in query
        or "," in query
        or bool(GO_ID_TOKEN_PATTERN.fullmatch(query))
    )


def load_accession_mf_go_terms(path: str | Path) -> dict[str, tuple[str, ...]]:
    """Load headerless accession-to-MF-GO TSV produced by build-library."""
    accession_terms: dict[str, tuple[str, ...]] = {}
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line.strip():
                continue
            columns = line.split("\t")
            if len(columns) < 2:
                raise ValueError(
                    f"{path} line {line_number} must contain accession and GO terms"
                )
            accession = columns[0].strip()
            if not accession:
                raise ValueError(f"{path} line {line_number} has empty accession")
            accession_terms[accession] = _deduplicate(GO_ID_PATTERN.findall(columns[1]))
    return accession_terms


def resolve_go_set_query(
    query: str,
    accession_terms: dict[str, tuple[str, ...]],
) -> tuple[str, ...]:
    """Resolve one parenthesized GO set or accession query."""
    if is_go_term_set_input(query):
        return parse_go_term_set(query)
    terms = accession_terms.get(query.strip())
    if not terms:
        raise ValueError(f"No GO terms found for accession: {query}")
    return terms


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
