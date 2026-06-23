"""Rule-based Molecular Function GO natural-language descriptions."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from prosig.go.build import MF_ROOT

BINDING_ROLES = {
    "binding",
    "binding_cofactor",
    "binding_nucleic_acid",
    "binding_generic",
}
HEAD_ROLE_ORDER = {
    "catalytic": 0,
    "transporter": 1,
    "receptor": 2,
    "transcription_factor": 3,
    "regulator": 4,
    "structural": 5,
    "motor": 6,
    "binding_cofactor": 7,
    "binding_nucleic_acid": 8,
    "binding_generic": 9,
    "binding": 10,
    "unknown": 11,
}
WEAK_BINDING_NAMES = {
    "binding",
    "protein binding",
    "ion binding",
    "metal ion binding",
    "small molecule binding",
}


@dataclass(frozen=True)
class GoFunctionTerm:
    """One resolved term in a function description."""

    go_id: str
    name: str
    role: str
    used: bool
    dropped: bool
    missing: bool = False


@dataclass(frozen=True)
class GoFunctionDescription:
    """Structured result for GO-set-to-function composition."""

    query: str
    terms: tuple[GoFunctionTerm, ...]
    summary: str
    head: str | None
    modifiers: tuple[str, ...]
    supporting_terms: tuple[str, ...]
    dropped_terms: tuple[str, ...]

    def asdict(self) -> dict[str, Any]:
        return asdict(self)


def describe_go_function(
    query: str,
    go_terms: tuple[str, ...] | list[str] | set[str],
    go_graph_terms: dict[str, dict[str, Any]],
    *,
    max_modifiers: int = 3,
) -> GoFunctionDescription:
    """Compose a conservative function summary from MF GO terms."""
    ordered_terms = tuple(dict.fromkeys(go_terms))
    valid_terms = tuple(go_id for go_id in ordered_terms if go_id in go_graph_terms)
    dropped_terms = _dropped_ancestor_terms(valid_terms, go_graph_terms)
    candidate_terms = tuple(
        go_id for go_id in valid_terms if go_id not in dropped_terms and go_id != MF_ROOT
    )

    term_rows = tuple(
        GoFunctionTerm(
            go_id=go_id,
            name=_term_name(go_id, go_graph_terms),
            role=_term_role(go_id, go_graph_terms),
            used=go_id in candidate_terms,
            dropped=go_id in dropped_terms,
            missing=go_id not in go_graph_terms,
        )
        for go_id in ordered_terms
    )

    head = _select_head(candidate_terms, go_graph_terms)
    if head is None:
        return GoFunctionDescription(
            query=query,
            terms=term_rows,
            summary=f"{query} has no resolvable Molecular Function GO terms.",
            head=None,
            modifiers=(),
            supporting_terms=(),
            dropped_terms=dropped_terms,
        )

    modifiers = _binding_modifiers(
        candidate_terms,
        go_graph_terms,
        skip={head},
        max_modifiers=max_modifiers,
    )
    head_phrase = _head_phrase(_term_name(head, go_graph_terms))
    supporting_terms = _supporting_terms(
        valid_terms,
        go_graph_terms,
        skip={head, *candidate_terms},
    )
    summary = _compose_sentence(query, head_phrase, modifiers, supporting_terms)
    return GoFunctionDescription(
        query=query,
        terms=term_rows,
        summary=summary,
        head=head,
        modifiers=modifiers,
        supporting_terms=supporting_terms,
        dropped_terms=dropped_terms,
    )


def _dropped_ancestor_terms(
    go_terms: tuple[str, ...],
    go_graph_terms: dict[str, dict[str, Any]],
) -> tuple[str, ...]:
    dropped: list[str] = []
    for go_id in go_terms:
        if go_id == MF_ROOT:
            dropped.append(go_id)
            continue
        for other_id in go_terms:
            if other_id == go_id:
                continue
            if go_id in go_graph_terms.get(other_id, {}).get("ancestors", ()):
                dropped.append(go_id)
                break
    return tuple(dict.fromkeys(dropped))


def _select_head(
    go_terms: tuple[str, ...],
    go_graph_terms: dict[str, dict[str, Any]],
) -> str | None:
    if not go_terms:
        return None
    return min(
        go_terms,
        key=lambda go_id: (
            HEAD_ROLE_ORDER.get(_term_role(go_id, go_graph_terms), 99),
            -_term_float(go_id, go_graph_terms, "ic"),
            -_term_int(go_id, go_graph_terms, "depth"),
            -len(_term_name(go_id, go_graph_terms)),
            go_id,
        ),
    )


def _binding_modifiers(
    go_terms: tuple[str, ...],
    go_graph_terms: dict[str, dict[str, Any]],
    *,
    skip: set[str],
    max_modifiers: int,
) -> tuple[str, ...]:
    if max_modifiers <= 0:
        return ()

    modifiers: list[str] = []
    for go_id in sorted(
        (term for term in go_terms if term not in skip),
        key=lambda term: (
            HEAD_ROLE_ORDER.get(_term_role(term, go_graph_terms), 99),
            -_term_float(term, go_graph_terms, "ic"),
            -_term_int(term, go_graph_terms, "depth"),
            term,
        ),
    ):
        role = _term_role(go_id, go_graph_terms)
        name = _term_name(go_id, go_graph_terms)
        if role not in BINDING_ROLES or name in WEAK_BINDING_NAMES:
            continue
        modifier = _binding_modifier(name)
        if modifier and modifier not in modifiers:
            modifiers.append(modifier)
        if len(modifiers) >= max_modifiers:
            break
    return tuple(modifiers)


def _supporting_terms(
    go_terms: tuple[str, ...],
    go_graph_terms: dict[str, dict[str, Any]],
    *,
    skip: set[str],
) -> tuple[str, ...]:
    supporting: list[str] = []
    for go_id in go_terms:
        if go_id in skip or go_id == MF_ROOT:
            continue
        role = _term_role(go_id, go_graph_terms)
        if role in BINDING_ROLES:
            continue
        phrase = _head_phrase(_term_name(go_id, go_graph_terms))
        if phrase and phrase not in supporting:
            supporting.append(phrase)
    return tuple(supporting)


def _compose_sentence(
    query: str,
    head_phrase: str,
    modifiers: tuple[str, ...],
    supporting_terms: tuple[str, ...],
) -> str:
    noun_phrase = " ".join((*_merge_binding_modifiers(modifiers), head_phrase)).strip()
    article = _article_for(noun_phrase)
    sentence = f"{query} is annotated as {article} {noun_phrase}"
    if supporting_terms:
        sentence += f" with {_join_phrases(supporting_terms)} activity"
    return f"{sentence}."


def _binding_modifier(name: str) -> str:
    if not name.endswith(" binding"):
        return ""
    ligand = name[: -len(" binding")].strip()
    if ligand.endswith(" ion"):
        ligand = ligand[: -len(" ion")].strip()
    return f"{ligand}-binding" if ligand else ""


def _head_phrase(name: str) -> str:
    for suffix in (" activity", " molecular function"):
        if name.endswith(suffix):
            return name[: -len(suffix)].strip()
    return name


def _merge_binding_modifiers(modifiers: tuple[str, ...]) -> tuple[str, ...]:
    if len(modifiers) <= 1:
        return modifiers
    if not all(modifier.endswith("-binding") for modifier in modifiers):
        return modifiers

    stems = tuple(modifier[: -len("-binding")] for modifier in modifiers)
    if len(stems) == 2:
        return (f"{stems[0]}- and {stems[1]}-binding",)
    prefix = ", ".join(f"{stem}-" for stem in stems[:-1])
    return (f"{prefix}, and {stems[-1]}-binding",)


def _term_name(go_id: str, go_graph_terms: dict[str, dict[str, Any]]) -> str:
    term = go_graph_terms.get(go_id)
    if term is None:
        return "missing"
    return str(term.get("name", ""))


def _term_role(go_id: str, go_graph_terms: dict[str, dict[str, Any]]) -> str:
    term = go_graph_terms.get(go_id)
    if term is None:
        return "missing"
    semantic_role = term.get("semantic_role")
    if isinstance(semantic_role, dict) and semantic_role.get("role"):
        return str(semantic_role["role"])
    name = _term_name(go_id, go_graph_terms)
    if "binding" in name:
        return "binding_generic"
    if "transcription factor" in name:
        return "transcription_factor"
    if "transporter" in name:
        return "transporter"
    if "receptor" in name:
        return "receptor"
    if name.endswith(" activity"):
        return "catalytic"
    return "unknown"


def _term_float(
    go_id: str,
    go_graph_terms: dict[str, dict[str, Any]],
    key: str,
) -> float:
    value = go_graph_terms.get(go_id, {}).get(key)
    if isinstance(value, (int, float)):
        return float(value)
    return -1.0


def _term_int(
    go_id: str,
    go_graph_terms: dict[str, dict[str, Any]],
    key: str,
) -> int:
    value = go_graph_terms.get(go_id, {}).get(key)
    if isinstance(value, int):
        return value
    return -1


def _article_for(phrase: str) -> str:
    if phrase[:1].lower() in {"a", "e", "i", "o", "u"}:
        return "an"
    return "a"


def _join_phrases(phrases: tuple[str, ...]) -> str:
    if len(phrases) <= 1:
        return "".join(phrases)
    if len(phrases) == 2:
        return " and ".join(phrases)
    return f"{', '.join(phrases[:-1])}, and {phrases[-1]}"
