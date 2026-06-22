from __future__ import annotations

import gzip
import json
import logging
import math
import pickle
import shutil
from collections import Counter, defaultdict, deque
from collections.abc import Iterator
from datetime import UTC, datetime
from pathlib import Path
from statistics import mean, median
from typing import Any

MF_ROOT = "GO:0003674"
MF_NAMESPACE = "molecular_function"
EXCLUDED_EVIDENCE = frozenset({"ND", "NAS"})
# This reviewed-evidence policy is intended for Swiss-Prot primary accessions
# from uniprot_sprot.dat.gz. Do not apply it blindly to unreviewed sources.
LOGGER_NAME = "prosig"


def build_go_pkl(
    *,
    go_obo: Path,
    swissprot: Path,
    go_out: Path,
    report_out: Path | None = None,
    role_map: Path | None = None,
) -> dict[str, Any]:
    logger = logging.getLogger(LOGGER_NAME)

    logger.info("Parsing GO OBO file for Molecular Function terms: %s", go_obo)
    terms, obsolete_go_ids = _parse_go_obo_mf(go_obo)
    logger.info(
        "Parsed %s connected Molecular Function GO terms",
        format_log_number(len(terms)),
    )

    logger.info("Parsing Swiss-Prot MF GO annotations: %s", swissprot)
    logger.info("Propagating MF GO annotations and calculating IC values")
    annotation_stats = apply_ic_from_swissprot(
        terms,
        swissprot,
        obsolete_go_ids=obsolete_go_ids,
    )
    logger.info(
        "Parsed %s primary accessions; %s have high-quality MF annotations",
        format_log_number(annotation_stats["n_accessions_provided"]),
        format_log_number(annotation_stats["n_accessions_with_hq_mf_go"]),
    )
    logger.info(
        "GO annotation accession summary: total=%s; MF=%s; MF high-quality=%s; "
        "BP=%s; BP high-quality=%s; CC=%s; CC high-quality=%s",
        format_log_number(annotation_stats["n_accessions_provided"]),
        format_log_number(annotation_stats["n_accessions_with_any_mf_go"]),
        format_log_number(annotation_stats["n_accessions_with_hq_mf_go"]),
        format_log_number(annotation_stats["n_accessions_with_any_bp_go"]),
        format_log_number(annotation_stats["n_accessions_with_hq_bp_go"]),
        format_log_number(annotation_stats["n_accessions_with_any_cc_go"]),
        format_log_number(annotation_stats["n_accessions_with_hq_cc_go"]),
    )
    logger.info(
        "Calculated IC values using %s accessions; skipped %s HQ MF GO "
        "assignments not in graph",
        format_log_number(annotation_stats["n_accessions_used_for_ic"]),
        format_log_number(annotation_stats["n_hq_mf_go_assignments_not_in_graph"]),
    )
    logger.info(
        "Skipped %s HQ MF GO assignments because the GO term is obsolete",
        format_log_number(annotation_stats["n_hq_mf_go_assignments_obsolete"]),
    )
    frequency_metadata = format_go_frequency_metadata(annotation_stats)
    logger.info(
        "%s GO terms did not receive valid IC because no accession matched them",
        format_log_number(count_terms_without_valid_ic(terms)),
    )
    log_top_frequency_terms(logger, terms, limit=10)

    role_stats = None
    if role_map is not None:
        unknown_role_out = go_out.parent / "go_terms_unknown_role.txt"
        logger.info("Loading GO semantic role map: %s", role_map)
        logger.info(
            "Assigning GO semantic roles to %s non-root GO terms",
            format_log_number(len(terms) - 1),
        )
        logger.info("Applying Layer 1 GO anchor/ancestor role matching")
        logger.info("Applying Layer 2 keyword role matching to remaining terms")
        role_stats = assign_semantic_roles_from_file(
            terms,
            role_map,
            unknown_role_out=unknown_role_out,
        )
        logger.info(
            "Processed %s GO terms for semantic role assignment",
            format_log_number(role_stats["processed"]),
        )
        logger.info(
            "GO semantic role layer summary:\n%s",
            format_semantic_role_layer_summary(role_stats),
        )
        logger.info(
            "GO semantic role stats:\n%s",
            format_semantic_role_stats(role_stats["role_counts"]),
        )
        logger.info("Wrote GO terms with unknown semantic role: %s", unknown_role_out)

    logger.info("Assembling GO graph metadata")
    meta = {
        "schema_version": "1.0",
        "namespace": MF_NAMESPACE,
        "source_obo": str(go_obo),
        "annotation_source": "Swiss-Prot",
        "ic_formula": "-log(freq)",
        "frequency_denominator": "accessions with at least one valid MF graph term",
        "propagated_counts": True,
        "obsolete_terms_removed": True,
        "n_terms": len(terms),
        "n_accessions_provided": annotation_stats["n_accessions_provided"],
        "n_accessions_with_hq_mf_go": annotation_stats["n_accessions_with_hq_mf_go"],
        "n_accessions_with_any_mf_go": annotation_stats[
            "n_accessions_with_any_mf_go"
        ],
        "n_hq_mf_go_assignments_not_in_graph": annotation_stats[
            "n_hq_mf_go_assignments_not_in_graph"
        ],
        "n_hq_mf_go_assignments_obsolete": annotation_stats[
            "n_hq_mf_go_assignments_obsolete"
        ],
        **frequency_metadata,
        "created_at": datetime.now(UTC).date().isoformat(),
    }
    if role_stats is not None:
        meta["semantic_role_assignment"] = {
            "role_map": str(role_map),
            "unknown_role_report": str(go_out.parent / "go_terms_unknown_role.txt"),
            "n_processed": role_stats["processed"],
            "n_anchor": role_stats["anchor"],
            "n_keyword": role_stats["keyword"],
            "n_unknown": role_stats["unknown"],
            "role_counts": dict(role_stats["role_counts"]),
        }
    artifact = {"meta": meta, "terms": terms}
    logger.info("Assembled GO graph metadata")

    logger.info("Writing GO graph and IC artifact: %s", go_out)
    go_out.parent.mkdir(parents=True, exist_ok=True)
    with go_out.open("wb") as handle:
        pickle.dump(artifact, handle, protocol=pickle.HIGHEST_PROTOCOL)
    logger.info("Wrote GO graph and IC artifact: %s", go_out)

    if report_out is not None:
        logger.info("Writing GO build validation report: %s", report_out)
        report_out.parent.mkdir(parents=True, exist_ok=True)
        report_out.write_text(format_go_report(artifact), encoding="utf-8")
        logger.info("Wrote GO build validation report: %s", report_out)

    go_json_out = go_out.parent / "go_graph.json"
    logger.info("Writing diagnostic GO graph JSON: %s", go_json_out)
    write_go_graph_json(go_json_out, artifact)
    logger.info("Wrote diagnostic GO graph JSON: %s", go_json_out)

    excluded_mf_annotations_out = go_out.parent / "excluded_mf_annotations.tsv"
    logger.info(
        "Writing excluded MF annotation diagnostics: %s",
        excluded_mf_annotations_out,
    )
    write_excluded_mf_annotation_diagnostics(
        excluded_mf_annotations_out,
        swissprot,
    )
    logger.info(
        "Wrote excluded MF annotation diagnostics: %s",
        excluded_mf_annotations_out,
    )

    accession_mf_go_out = go_out.parent / "accession_mf_go.tsv"
    logger.info("Writing accession MF GO terms: %s", accession_mf_go_out)
    write_accession_mf_go_tsv(accession_mf_go_out, swissprot)
    logger.info("Wrote accession MF GO terms: %s", accession_mf_go_out)

    return artifact


def parse_go_obo_mf(path: Path) -> dict[str, dict[str, Any]]:
    terms, _obsolete_go_ids = _parse_go_obo_mf(path)
    return terms


def _parse_go_obo_mf(path: Path) -> tuple[dict[str, dict[str, Any]], set[str]]:
    all_terms: dict[str, dict[str, Any]] = {}
    obsolete_go_ids: set[str] = set()
    current: dict[str, Any] | None = None
    in_term = False

    def flush(term: dict[str, Any] | None) -> None:
        if not term or not term.get("id"):
            return
        if term.get("is_obsolete"):
            obsolete_go_ids.add(term["id"])
            return
        all_terms[term["id"]] = {
            "name": term.get("name", ""),
            "namespace": term.get("namespace", ""),
            "parents": set(term.get("parents", set())),
        }

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line == "[Term]":
                flush(current)
                current = {"parents": set()}
                in_term = True
                continue
            if line.startswith("[") and line != "[Term]":
                flush(current)
                current = None
                in_term = False
                continue
            if not in_term or current is None:
                continue
            if line.startswith("id:"):
                current["id"] = line.split(":", 1)[1].strip()
            elif line.startswith("name:"):
                current["name"] = line.split(":", 1)[1].strip()
            elif line.startswith("namespace:"):
                current["namespace"] = line.split(":", 1)[1].strip()
            elif line.startswith("is_a:"):
                current["parents"].add(line.split("is_a:", 1)[1].strip().split()[0])
            elif line == "is_obsolete: true":
                current["is_obsolete"] = True
        flush(current)

    terms = {
        go_id: {
            "name": term["name"],
            "parents": sorted(
                parent
                for parent in term["parents"]
                if parent in all_terms
                and all_terms[parent]["namespace"] == MF_NAMESPACE
            ),
            "children": [],
            "ancestors": set(),
            "depth": 0,
            "count": 0,
            "freq": 0.0,
            "ic": None,
        }
        for go_id, term in all_terms.items()
        if term["namespace"] == MF_NAMESPACE
    }

    children: dict[str, list[str]] = defaultdict(list)
    for go_id, term in terms.items():
        for parent in term["parents"]:
            children[parent].append(go_id)
    for go_id, child_ids in children.items():
        terms[go_id]["children"] = sorted(child_ids)

    if MF_ROOT not in terms:
        raise ValueError(f"Required Molecular Function root is missing: {MF_ROOT}")

    depths = _compute_depths(terms)
    disconnected = set(terms) - set(depths)
    if disconnected:
        for go_id in disconnected:
            del terms[go_id]
        for term in terms.values():
            term["parents"] = [parent for parent in term["parents"] if parent in terms]
            term["children"] = [child for child in term["children"] if child in terms]

    ancestor_cache: dict[str, set[str]] = {}
    for go_id in terms:
        terms[go_id]["ancestors"] = _ancestors(go_id, terms, ancestor_cache)
        terms[go_id]["depth"] = depths[go_id]

    return terms, obsolete_go_ids


def parse_swissprot_mf_go(path: Path) -> tuple[dict[str, set[str]], dict[str, int]]:
    accession_to_terms: dict[str, set[str]] = {}
    stats = {
        "n_accessions_provided": 0,
        "n_accessions_used_for_ic": 0,
        "n_hq_mf_go_assignments_not_in_graph": 0,
        "n_accessions_with_any_mf_go": 0,
        "n_accessions_with_hq_mf_go": 0,
        "n_accessions_with_any_bp_go": 0,
        "n_accessions_with_hq_bp_go": 0,
        "n_accessions_with_any_cc_go": 0,
        "n_accessions_with_hq_cc_go": 0,
    }

    for accession, annotation in iter_swissprot_go_annotations(path):
        mf_terms = annotation.high_quality_mf_terms
        stats["n_accessions_provided"] += 1
        _update_go_namespace_stats(stats, annotation)
        if mf_terms:
            accession_to_terms[accession] = mf_terms

    return accession_to_terms, stats


def iter_swissprot_mf_go(path: Path) -> Iterator[tuple[str, set[str]]]:
    for accession, annotation in iter_swissprot_go_annotations(path):
        yield accession, annotation.high_quality_mf_terms


class SwissProtGoAnnotation:
    def __init__(self) -> None:
        self.mf_go_terms: list[tuple[str, str]] = []
        self.high_quality_mf_terms: set[str] = set()
        self.has_mf_go = False
        self.has_high_quality_mf_go = False
        self.has_bp_go = False
        self.has_high_quality_bp_go = False
        self.has_cc_go = False
        self.has_high_quality_cc_go = False


def iter_swissprot_go_annotations(
    path: Path,
) -> Iterator[tuple[str, SwissProtGoAnnotation]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        primary_accession: str | None = None
        annotation = SwissProtGoAnnotation()

        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line == "//":
                if primary_accession is not None:
                    yield primary_accession, annotation
                primary_accession = None
                annotation = SwissProtGoAnnotation()
                continue

            if line.startswith("AC"):
                if primary_accession is None:
                    accessions = [
                        accession.strip()
                        for accession in line[5:].strip().split(";")
                        if accession.strip()
                    ]
                    if accessions:
                        primary_accession = accessions[0]
                continue

            _update_annotation_from_go_line(annotation, line)

        if primary_accession is not None:
            yield primary_accession, annotation


def parse_swissprot_entry(lines: list[str]) -> tuple[str | None, set[str]]:
    accessions: list[str] = []
    mf_terms: set[str] = set()

    for line in lines:
        if line.startswith("AC"):
            accessions.extend(
                accession.strip()
                for accession in line[5:].strip().split(";")
                if accession.strip()
            )
            continue
        go_id = _parse_high_quality_mf_go_id(line)
        if go_id is not None:
            mf_terms.add(go_id)

    return (accessions[0] if accessions else None), mf_terms


def _parse_high_quality_mf_go_id(line: str) -> str | None:
    parsed = _parse_go_dr_line(line)
    if parsed is None:
        return None
    go_id, namespace, evidence = parsed

    if (
        go_id.startswith("GO:")
        and namespace == "F"
        and is_high_quality_evidence(evidence)
    ):
        return go_id
    return None


def _parse_go_dr_line(line: str) -> tuple[str, str, str] | None:
    if not line.startswith("DR   GO;"):
        return None

    parts = [part.strip() for part in line.split(";")]
    if len(parts) < 4:
        return None

    go_id = parts[1]
    namespace = parts[2].split(":", 1)[0].strip()
    evidence = parts[3].split(":", 1)[0].strip()
    if not go_id.startswith("GO:"):
        return None
    return go_id, namespace, evidence


def _update_annotation_from_go_line(
    annotation: SwissProtGoAnnotation,
    line: str,
) -> None:
    parsed = _parse_go_dr_line(line)
    if parsed is None:
        return

    go_id, namespace, evidence = parsed
    is_high_quality = is_high_quality_evidence(evidence)
    if namespace == "F":
        annotation.has_mf_go = True
        annotation.mf_go_terms.append((go_id, evidence))
        if is_high_quality:
            annotation.has_high_quality_mf_go = True
            annotation.high_quality_mf_terms.add(go_id)
    elif namespace == "P":
        annotation.has_bp_go = True
        if is_high_quality:
            annotation.has_high_quality_bp_go = True
    elif namespace == "C":
        annotation.has_cc_go = True
        if is_high_quality:
            annotation.has_high_quality_cc_go = True


def apply_ic_from_accessions(
    terms: dict[str, dict[str, Any]],
    accession_to_terms: dict[str, set[str]],
) -> dict[str, Any]:
    counts = {go_id: 0 for go_id in terms}
    skipped_go_term_counts: Counter[str] = Counter()
    obsolete_go_term_counts: Counter[str] = Counter()
    used_accessions = 0

    for direct_terms in accession_to_terms.values():
        used_accessions += _apply_accession_terms_to_counts(
            terms,
            counts,
            direct_terms,
            skipped_go_term_counts,
            obsolete_go_term_counts,
            obsolete_go_ids=set(),
        )

    _finalize_ic_terms(terms, counts, used_accessions)

    return {
        "n_accessions_used_for_ic": used_accessions,
        "n_hq_mf_go_assignments_not_in_graph": sum(
            skipped_go_term_counts.values()
        ),
        "n_hq_mf_go_assignments_obsolete": sum(obsolete_go_term_counts.values()),
        "skipped_go_term_counts": skipped_go_term_counts,
        "obsolete_go_term_counts": obsolete_go_term_counts,
    }


def apply_ic_from_swissprot(
    terms: dict[str, dict[str, Any]],
    path: Path,
    *,
    obsolete_go_ids: set[str],
) -> dict[str, Any]:
    counts = {go_id: 0 for go_id in terms}
    direct_go_term_counts: Counter[str] = Counter()
    skipped_go_term_counts: Counter[str] = Counter()
    obsolete_go_term_counts: Counter[str] = Counter()
    stats: dict[str, Any] = {
        "n_accessions_provided": 0,
        "n_accessions_with_any_mf_go": 0,
        "n_accessions_with_hq_mf_go": 0,
        "n_accessions_with_any_bp_go": 0,
        "n_accessions_with_hq_bp_go": 0,
        "n_accessions_with_any_cc_go": 0,
        "n_accessions_with_hq_cc_go": 0,
        "n_accessions_used_for_ic": 0,
        "n_hq_mf_go_assignments_not_in_graph": 0,
        "n_hq_mf_go_assignments_obsolete": 0,
        "direct_go_term_counts": direct_go_term_counts,
        "skipped_go_term_counts": skipped_go_term_counts,
        "obsolete_go_term_counts": obsolete_go_term_counts,
    }

    for _accession, annotation in iter_swissprot_go_annotations(path):
        direct_terms = annotation.high_quality_mf_terms
        stats["n_accessions_provided"] += 1
        _update_go_namespace_stats(stats, annotation)
        if direct_terms:
            direct_go_term_counts.update(direct_terms)
        stats["n_accessions_used_for_ic"] += _apply_accession_terms_to_counts(
            terms,
            counts,
            direct_terms,
            skipped_go_term_counts,
            obsolete_go_term_counts,
            obsolete_go_ids=obsolete_go_ids,
        )

    _finalize_ic_terms(terms, counts, stats["n_accessions_used_for_ic"])
    stats["n_hq_mf_go_assignments_not_in_graph"] = sum(
        skipped_go_term_counts.values()
    )
    stats["n_hq_mf_go_assignments_obsolete"] = sum(
        obsolete_go_term_counts.values()
    )
    return stats


def _update_go_namespace_stats(
    stats: dict[str, Any],
    annotation: SwissProtGoAnnotation,
) -> None:
    if annotation.has_mf_go:
        stats["n_accessions_with_any_mf_go"] += 1
    if annotation.has_high_quality_mf_go:
        stats["n_accessions_with_hq_mf_go"] += 1
    if annotation.has_bp_go:
        stats["n_accessions_with_any_bp_go"] += 1
    if annotation.has_high_quality_bp_go:
        stats["n_accessions_with_hq_bp_go"] += 1
    if annotation.has_cc_go:
        stats["n_accessions_with_any_cc_go"] += 1
    if annotation.has_high_quality_cc_go:
        stats["n_accessions_with_hq_cc_go"] += 1


def format_go_frequency_metadata(annotation_stats: dict[str, Any]) -> dict[str, Any]:
    direct_go_term_counts = annotation_stats["direct_go_term_counts"]
    frequencies = list(direct_go_term_counts.values())
    return _frequency_stats(frequencies)


def count_terms_without_valid_ic(terms: dict[str, dict[str, Any]]) -> int:
    return sum(1 for term in terms.values() if term["ic"] is None)


def log_top_frequency_terms(
    logger: logging.Logger,
    terms: dict[str, dict[str, Any]],
    *,
    limit: int,
) -> None:
    logger.info("Top %d most frequent MF GO terms:", limit)
    for rank, (go_id, term) in enumerate(
        _ranked_terms(terms, key="count", reverse=True, limit=limit),
        start=1,
    ):
        logger.info(
            "%s. %s %s count=%s freq=%s ic=%s",
            format_log_number(rank),
            go_id,
            term["name"],
            format_log_number(term["count"]),
            format_log_number(term["freq"]),
            format_log_number(term["ic"]),
        )


def format_log_number(value: int | float | None) -> str:
    if value is None:
        return "None"
    if isinstance(value, int):
        return f"{value:,}"
    if value == 0:
        return "0"
    return f"{value:,.4f}".rstrip("0").rstrip(".")


def _frequency_stats(frequencies: list[int]) -> dict[str, int | float | str]:
    if not frequencies:
        return {
            "mf_frequency_min": 0,
            "mf_frequency_median": 0,
            "mf_frequency_mean": 0.0,
            "mf_frequency_max": 0,
            "mf_frequency_status": "EMPTY",
        }

    return {
        "mf_frequency_min": min(frequencies),
        "mf_frequency_median": median(frequencies),
        "mf_frequency_mean": mean(frequencies),
        "mf_frequency_max": max(frequencies),
        "mf_frequency_status": "OK",
    }


def write_go_graph_json(path: Path, artifact: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    diagnostic_artifact = {
        "_comment": "Diagnostic only. Use go_graph.pkl as the runtime artifact.",
        **artifact,
    }
    path.write_text(
        json.dumps(_json_ready(diagnostic_artifact), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: _json_ready(item) for key, item in value.items()}
    if isinstance(value, set):
        return sorted(value)
    if isinstance(value, list):
        return [_json_ready(item) for item in value]
    return value


def write_excluded_mf_annotation_diagnostics(
    path: Path,
    swissprot: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("accession\tgo_term\tevidence\n")
        for accession, annotation in iter_swissprot_go_annotations(swissprot):
            for go_id, evidence in annotation.mf_go_terms:
                if evidence in EXCLUDED_EVIDENCE:
                    handle.write(f"{accession}\t{go_id}\t{evidence}\n")


def write_accession_mf_go_tsv(
    path: Path,
    swissprot: Path,
) -> None:
    """Write primary Swiss-Prot accessions and high-quality MF GO terms."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for accession, terms in iter_swissprot_mf_go(swissprot):
            if not terms:
                continue
            handle.write(f"{accession}\t{';'.join(sorted(terms))}\n")


def ensure_role_map_from_template(
    role_map: Path,
    template: Path,
) -> bool:
    """Create role_map from template if missing; return True when created."""
    if role_map.exists():
        return False
    role_map.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(template, role_map)
    return True


def assign_semantic_roles_from_file(
    terms: dict[str, dict[str, Any]],
    role_map: Path,
    *,
    unknown_role_out: Path,
) -> dict[str, Any]:
    role_config = parse_role_map(role_map)
    return assign_semantic_roles(
        terms,
        role_config,
        unknown_role_out=unknown_role_out,
    )


def assign_semantic_roles(
    terms: dict[str, dict[str, Any]],
    role_config: dict[str, dict[str, dict[str, Any]]],
    *,
    unknown_role_out: Path,
) -> dict[str, Any]:
    anchor_index = _build_anchor_index(role_config)
    keyword_rules = _build_keyword_rules(role_config)
    role_counts: Counter[str] = Counter()
    stats: dict[str, Any] = {
        "processed": 0,
        "anchor": 0,
        "keyword": 0,
        "unknown": 0,
        "role_counts": role_counts,
    }
    unknown_terms: list[tuple[str, str]] = []

    for go_id, term in terms.items():
        if go_id == MF_ROOT:
            continue
        stats["processed"] += 1

        role = _semantic_role_from_anchors(go_id, term, terms, anchor_index)
        if role is None:
            role = _semantic_role_from_keywords(term["name"], keyword_rules)

        if role is None:
            stats["unknown"] += 1
            unknown_terms.append((go_id, term["name"]))
            term["semantic_role"] = {
                "role": "unknown",
                "priority": 0,
                "source": "unknown",
                "matched": None,
            }
            role_counts["unknown"] += 1
            continue

        stats[role["source"]] += 1
        role_counts[role["role"]] += 1
        term["semantic_role"] = {
            "role": role["role"],
            "priority": role["priority"],
            "source": role["source"],
            "matched": role["matched"],
        }

    unknown_role_out.parent.mkdir(parents=True, exist_ok=True)
    unknown_role_out.write_text(
        "".join(f"{go_id}: {name}\n" for go_id, name in sorted(unknown_terms)),
        encoding="utf-8",
    )
    return stats


def format_semantic_role_stats(role_counts: Counter[str] | dict[str, int]) -> str:
    parts = [
        f"  {role} = {format_log_number(count)}"
        for role, count in sorted(role_counts.items())
        if role != "unknown"
    ]
    if "unknown" in role_counts:
        parts.append(f"  unknown = {format_log_number(role_counts['unknown'])}")
    return "\n".join(parts)


def format_semantic_role_layer_summary(role_stats: dict[str, Any]) -> str:
    numbers = {
        "processed": format_log_number(role_stats["processed"]),
        "anchor": format_log_number(role_stats["anchor"]),
        "keyword": format_log_number(role_stats["keyword"]),
        "unknown": format_log_number(role_stats["unknown"]),
    }
    width = max(7, *(len(value) for value in numbers.values()))
    return "\n".join(
        [
            f"  total non-root terms = {numbers['processed']:>{width}}",
            f"  anchor assigned      = {numbers['anchor']:>{width}}",
            f"  keyword assigned     = {numbers['keyword']:>{width}}",
            f"  unknown              = {numbers['unknown']:>{width}}",
        ]
    )


def parse_role_map(path: Path) -> dict[str, dict[str, dict[str, Any]]]:
    """Parse the restricted role_map.yaml structure used by ProSig."""
    config: dict[str, dict[str, dict[str, Any]]] = {
        "roles": {},
        "role_rules": {},
    }
    section: str | None = None
    role: str | None = None
    list_key: str | None = None

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        if indent == 0 and stripped.endswith(":"):
            section = stripped[:-1]
            role = None
            list_key = None
            if section not in config:
                raise ValueError(f"Unsupported role map section: {section}")
            continue
        if section is None:
            raise ValueError(f"Role map entry outside a section: {stripped}")
        if indent == 2 and stripped.endswith(":"):
            role = stripped[:-1]
            config[section][role] = {}
            list_key = None
            continue
        if role is None:
            raise ValueError(f"Role map field outside a role: {stripped}")
        if indent == 4 and ":" in stripped:
            key, value = [item.strip() for item in stripped.split(":", 1)]
            if value:
                config[section][role][key] = int(value) if key == "priority" else value
                list_key = None
            else:
                config[section][role][key] = []
                list_key = key
            continue
        if indent == 6 and stripped.startswith("- "):
            if list_key is None:
                raise ValueError(f"Role map list item outside a list: {stripped}")
            config[section][role][list_key].append(_unquote_yaml_value(stripped[2:]))
            continue
        raise ValueError(f"Unsupported role map line: {raw_line}")

    return config


def _build_anchor_index(
    role_config: dict[str, dict[str, dict[str, Any]]],
) -> dict[str, dict[str, Any]]:
    anchor_index: dict[str, dict[str, Any]] = {}
    for role, config in role_config["roles"].items():
        priority = int(config.get("priority", 0))
        for anchor in config.get("anchors", []):
            current = anchor_index.get(anchor)
            if current is None or priority > current["priority"]:
                anchor_index[anchor] = {"role": role, "priority": priority}
    return anchor_index


def _build_keyword_rules(
    role_config: dict[str, dict[str, dict[str, Any]]],
) -> list[dict[str, Any]]:
    rules: list[dict[str, Any]] = []
    for section in ("roles", "role_rules"):
        for role, config in role_config[section].items():
            priority = int(config.get("priority", 0))
            for keyword in config.get("keywords", []):
                rules.append(
                    {
                        "role": role,
                        "priority": priority,
                        "keyword": keyword,
                        "keyword_lower": keyword.lower(),
                    }
                )
    return sorted(rules, key=lambda rule: rule["priority"], reverse=True)


def _semantic_role_from_anchors(
    go_id: str,
    term: dict[str, Any],
    terms: dict[str, dict[str, Any]],
    anchor_index: dict[str, dict[str, Any]],
) -> dict[str, Any] | None:
    for candidate in _rank_self_and_ancestors_by_ic(go_id, term, terms):
        role = anchor_index.get(candidate)
        if role is not None:
            return {
                "role": role["role"],
                "priority": role["priority"],
                "source": "anchor",
                "matched": candidate,
            }
    return None


def _semantic_role_from_keywords(
    name: str,
    keyword_rules: list[dict[str, Any]],
) -> dict[str, Any] | None:
    name_lower = name.lower()
    for rule in keyword_rules:
        if rule["keyword_lower"] in name_lower:
            return {
                "role": rule["role"],
                "priority": rule["priority"],
                "source": "keyword",
                "matched": rule["keyword"],
            }
    return None


def _rank_self_and_ancestors_by_ic(
    go_id: str,
    term: dict[str, Any],
    terms: dict[str, dict[str, Any]],
) -> list[str]:
    candidates = [go_id, *term.get("ancestors", set())]
    return sorted(
        candidates,
        key=lambda candidate: (
            terms.get(candidate, {}).get("ic") is not None,
            terms.get(candidate, {}).get("ic") or float("-inf"),
            candidate,
        ),
        reverse=True,
    )


def _unquote_yaml_value(value: str) -> str:
    value = _strip_yaml_inline_comment(value).strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    return value


def _strip_yaml_inline_comment(value: str) -> str:
    quote: str | None = None
    escaped = False
    for index, char in enumerate(value):
        if escaped:
            escaped = False
            continue
        if char == "\\":
            escaped = True
            continue
        if quote is not None:
            if char == quote:
                quote = None
            continue
        if char in {"'", '"'}:
            quote = char
            continue
        if char == "#" and (index == 0 or value[index - 1].isspace()):
            return value[:index]
    return value


def is_high_quality_evidence(evidence: str) -> bool:
    return evidence not in EXCLUDED_EVIDENCE


def _apply_accession_terms_to_counts(
    terms: dict[str, dict[str, Any]],
    counts: dict[str, int],
    direct_terms: set[str],
    skipped_go_term_counts: Counter[str],
    obsolete_go_term_counts: Counter[str],
    *,
    obsolete_go_ids: set[str],
) -> int:
    propagated_terms: set[str] = set()
    for go_id in direct_terms:
        if go_id in obsolete_go_ids:
            skipped_go_term_counts[go_id] += 1
            obsolete_go_term_counts[go_id] += 1
            continue
        if go_id not in terms:
            skipped_go_term_counts[go_id] += 1
            continue
        propagated_terms.add(go_id)
        propagated_terms.update(terms[go_id]["ancestors"])

    if not propagated_terms:
        return 0
    for go_id in propagated_terms:
        counts[go_id] += 1
    return 1


def _finalize_ic_terms(
    terms: dict[str, dict[str, Any]],
    counts: dict[str, int],
    used_accessions: int,
) -> None:
    if used_accessions == 0:
        raise ValueError("No Swiss-Prot accessions with valid MF GO terms were found.")

    for go_id, term in terms.items():
        count = counts[go_id]
        term["count"] = count
        if count == 0:
            term["freq"] = 0.0
            term["ic"] = None
            continue
        freq = count / used_accessions
        term["freq"] = freq
        term["ic"] = -math.log(freq)

    terms[MF_ROOT]["freq"] = 1.0
    terms[MF_ROOT]["ic"] = 0.0


def build_go_pkl_from_parsed(
    terms: dict[str, dict[str, Any]],
    accession_to_terms: dict[str, set[str]],
) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    return terms, apply_ic_from_accessions(terms, accession_to_terms)


def format_go_report(artifact: dict[str, Any]) -> str:
    terms = artifact["terms"]
    counted_terms = [term for term in terms.values() if term["count"] > 0]
    ic_terms = [term for term in terms.values() if term["ic"] is not None]
    max_depth = max((term["depth"] for term in terms.values()), default=0)
    root_terms = [go_id for go_id, term in terms.items() if not term["parents"]]
    accessions_used_for_ic = terms[MF_ROOT]["count"]
    accessions_skipped_no_valid_mf = (
        artifact["meta"]["n_accessions_provided"] - accessions_used_for_ic
    )

    lines = [
        f"number of MF terms: {len(terms)}",
        f"number of root terms: {len(root_terms)}",
        f"number of provided accessions: {artifact['meta']['n_accessions_provided']}",
        "number of accessions with HQ MF GO annotation: "
        f"{artifact['meta']['n_accessions_with_hq_mf_go']}",
        "number of accessions used for IC: "
        f"{accessions_used_for_ic}",
        "number of accessions skipped because no valid MF term remained after graph "
        f"filtering: {accessions_skipped_no_valid_mf}",
        "number of HQ MF GO assignments skipped because they were not in the MF graph: "
        f"{artifact['meta']['n_hq_mf_go_assignments_not_in_graph']}",
        "number of HQ MF GO assignments skipped because the GO term is obsolete: "
        f"{artifact['meta']['n_hq_mf_go_assignments_obsolete']}",
        f"number of terms with count > 0: {len(counted_terms)}",
        f"number of terms with IC value: {len(ic_terms)}",
        f"maximum depth: {max_depth}",
        "",
        "top 20 most frequent MF terms:",
    ]
    lines.extend(_format_ranked_terms(terms, key="count", reverse=True))
    lines.append("")
    lines.append("top 20 highest-IC terms with nonzero count:")
    lines.extend(
        _format_ranked_terms(
            {
                go_id: term
                for go_id, term in terms.items()
                if term["count"] > 0 and term["ic"] is not None
            },
            key="ic",
            reverse=True,
        )
    )
    return "\n".join(lines) + "\n"


def _compute_depths(terms: dict[str, dict[str, Any]]) -> dict[str, int]:
    depths: dict[str, int] = {}
    queue: deque[tuple[str, int]] = deque([(MF_ROOT, 0)])
    while queue:
        go_id, depth = queue.popleft()
        if go_id in depths and depths[go_id] <= depth:
            continue
        depths[go_id] = depth
        for child in terms[go_id]["children"]:
            queue.append((child, depth + 1))
    return depths


def _ancestors(
    go_id: str,
    terms: dict[str, dict[str, Any]],
    cache: dict[str, set[str]],
) -> set[str]:
    if go_id in cache:
        return cache[go_id]
    ancestors: set[str] = set()
    for parent in terms[go_id]["parents"]:
        ancestors.add(parent)
        ancestors.update(_ancestors(parent, terms, cache))
    cache[go_id] = ancestors
    return ancestors


def _format_ranked_terms(
    terms: dict[str, dict[str, Any]],
    *,
    key: str,
    reverse: bool,
) -> list[str]:
    return [
        f"{go_id}\t{term['name']}\tcount={term['count']}\tfreq={term['freq']:.6g}"
        f"\tic={term['ic'] if term['ic'] is not None else 'None'}"
        for go_id, term in _ranked_terms(terms, key=key, reverse=reverse, limit=20)
    ]


def _ranked_terms(
    terms: dict[str, dict[str, Any]],
    *,
    key: str,
    reverse: bool,
    limit: int,
) -> list[tuple[str, dict[str, Any]]]:
    return sorted(
        terms.items(),
        key=lambda item: (item[1][key] if item[1][key] is not None else -1, item[0]),
        reverse=reverse,
    )[:limit]
