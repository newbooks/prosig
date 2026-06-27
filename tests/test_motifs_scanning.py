import logging
from pathlib import Path

import pytest

import prosig.motifs.scanning as scanning
from prosig.motifs.scanning import (
    Motif,
    compile_prosig_motif,
    extract_motif_features,
    motif_features_complete,
    motif_literal_prefilter,
    motif_present,
    motif_width_bounds,
    read_prosig_motif_library,
    write_motif_features,
)


def test_read_prosig_motif_library_uses_prosig_pattern_column(tmp_path: Path) -> None:
    motif_file = tmp_path / "prosig_motifs.tsv"
    motif_file.write_text(
        "# ProSig motif library\n"
        "name\tprosite_ac\tdescription\tprosite_pattern\tprosig_pattern\tstatus\n"
        "ASN_GLYCOSYLATION\tPS00001\tN-glycosylation site\t"
        "N-{P}-[ST]-{P}\tN!P[ST]!P\tprosite\n",
        encoding="utf-8",
    )

    motifs = read_prosig_motif_library(motif_file)

    assert len(motifs) == 1
    assert motifs[0].name == "ASN_GLYCOSYLATION"
    assert motifs[0].pattern == "N!P[ST]!P"
    assert motifs[0].status == "prosite"
    assert motifs[0].compiled_pattern.match("NATS") is not None
    assert motifs[0].compiled_pattern.match("NPSS") is None


def test_compile_prosig_motif_supports_macros_and_repeats() -> None:
    motif = Motif(
        name="PKC",
        description="",
        pattern="[ST]?{+}",
        status="prosite",
        compiled_pattern=compile_prosig_motif("[ST]?{+}"),
    )

    result = extract_motif_features("AAASAKSSK", motif)

    assert result.count == 2
    assert result.first_position == 4
    assert result.last_position == 3
    assert result.match_fraction == pytest.approx(2 / 9)


def test_read_prosig_motif_library_marks_c_terminal_anchor(
    tmp_path: Path,
) -> None:
    motif_file = tmp_path / "prosig_motifs.tsv"
    motif_file.write_text(
        "name\tdescription\tprosig_pattern\tstatus\n"
        "C_TERM\tC-terminal motif\tA?(0,2)B>\tprosig\n",
        encoding="utf-8",
    )

    motifs = read_prosig_motif_library(motif_file)

    assert motifs[0].c_terminal_anchor
    assert motifs[0].min_width == 2
    assert motifs[0].max_width == 4
    assert motif_present("XXAAB", motifs[0])


def test_motif_present_scans_c_terminal_motif_from_suffix_only() -> None:
    compiled_pattern = _RecordingPattern(match_start=2)
    motif = Motif(
        name="C_TERM",
        description="",
        pattern="ABC>",
        status="prosig",
        compiled_pattern=compiled_pattern,
        c_terminal_anchor=True,
        min_width=3,
        max_width=3,
    )

    assert motif_present("XXABC", motif)
    assert compiled_pattern.starts == [2]


def test_motif_width_bounds_handles_terminal_set_as_zero_or_one_width() -> None:
    assert motif_width_bounds("AB[CD>]") == (2, 3)


def test_motif_literal_prefilter_uses_required_contiguous_literals() -> None:
    assert motif_literal_prefilter("AB[CD]EF?G") == "AB"
    assert motif_literal_prefilter("A[CD]E") is None


def test_motif_present_uses_literal_prefilter_before_regex() -> None:
    compiled_pattern = _RecordingPattern(match_start=0)
    motif = Motif(
        name="PREFILTERED",
        description="",
        pattern="AB[CD]EF",
        status="prosig",
        compiled_pattern=compiled_pattern,
        literal_prefilter="AB",
    )

    assert not motif_present("XXEFXX", motif)
    assert compiled_pattern.starts == []


def test_write_motif_features_scans_cluster_members_against_fasta(
    tmp_path: Path,
) -> None:
    cluster_file = tmp_path / "clusters.tsv"
    motif_file = tmp_path / "prosig_motifs.tsv"
    fasta_file = tmp_path / "accession.fasta"
    index_file = tmp_path / "accession.fasta.idx"
    output_file = tmp_path / "motif_features.tsv"
    cluster_file.write_text(
        "member_id\tcluster_id\nP1\tcluster_0001\nP2\tcluster_0001\n",
        encoding="utf-8",
    )
    motif_file.write_text(
        "name\tdescription\tprosig_pattern\tstatus\n"
        "N_GLY\tN-glycosylation\tN!P[ST]!P\tprosig\n",
        encoding="utf-8",
    )
    _write_indexed_fasta(
        fasta_file,
        index_file,
        {
            "P1": "AAANATSAA",
            "P2": "NNPSS",
        },
    )

    result = write_motif_features(
        cluster_file=cluster_file,
        motif_file=motif_file,
        fasta_file=fasta_file,
        fasta_index_file=index_file,
        output_file=output_file,
    )

    assert result.feature_rows == 1
    assert output_file.read_text(encoding="utf-8") == (
        "accession\tmotif_id\n"
        "P1\tN_GLY\n"
        "# completed\ttrue\n"
    )
    assert motif_features_complete(output_file)


def test_motif_features_complete_rejects_partial_file(tmp_path: Path) -> None:
    output_file = tmp_path / "motif_features.tsv"
    output_file.write_text(
        "accession\tmotif_id\n"
        "P1\tN_GLY\n",
        encoding="utf-8",
    )

    assert not motif_features_complete(output_file)


def test_write_motif_features_logs_scan_progress_on_interval(
    tmp_path: Path,
    monkeypatch,
    caplog,
) -> None:
    cluster_file = tmp_path / "clusters.tsv"
    motif_file = tmp_path / "prosig_motifs.tsv"
    fasta_file = tmp_path / "accession.fasta"
    index_file = tmp_path / "accession.fasta.idx"
    output_file = tmp_path / "motif_features.tsv"
    cluster_file.write_text(
        "member_id\tcluster_id\nP1\tcluster_0001\nP2\tcluster_0001\n",
        encoding="utf-8",
    )
    motif_file.write_text(
        "name\tdescription\tprosig_pattern\tstatus\n"
        "N_GLY\tN-glycosylation\tN!P[ST]!P\tprosig\n",
        encoding="utf-8",
    )
    _write_indexed_fasta(
        fasta_file,
        index_file,
        {
            "P1": "AAANATSAA",
            "P2": "AAANATSAA",
        },
    )
    monotonic_values = iter([0.0, 61.0, 61.0, 122.0, 122.0])
    monkeypatch.setattr(
        scanning.time,
        "monotonic",
        lambda: next(monotonic_values),
    )

    with caplog.at_level(logging.INFO):
        write_motif_features(
            cluster_file=cluster_file,
            motif_file=motif_file,
            fasta_file=fasta_file,
            fasta_index_file=index_file,
            output_file=output_file,
            processes=1,
            progress_interval_seconds=60.0,
        )

    assert "Motif scan progress: 2/2 accessions completed" in caplog.text


def test_write_motif_features_uses_multiple_processes_with_deterministic_output(
    tmp_path: Path,
) -> None:
    cluster_file = tmp_path / "clusters.tsv"
    motif_file = tmp_path / "prosig_motifs.tsv"
    fasta_file = tmp_path / "accession.fasta"
    index_file = tmp_path / "accession.fasta.idx"
    output_file = tmp_path / "motif_features.tsv"
    cluster_file.write_text(
        "member_id\tcluster_id\nP1\tcluster_0001\nP2\tcluster_0001\n",
        encoding="utf-8",
    )
    motif_file.write_text(
        "name\tdescription\tprosig_pattern\tstatus\n"
        "MOTIF_A\tA motif\tAA\tprosig\n"
        "MOTIF_B\tB motif\tBB\tprosig\n",
        encoding="utf-8",
    )
    _write_indexed_fasta(
        fasta_file,
        index_file,
        {
            "P1": "XXBBXXAA",
            "P2": "XXBBXX",
        },
    )

    result = write_motif_features(
        cluster_file=cluster_file,
        motif_file=motif_file,
        fasta_file=fasta_file,
        fasta_index_file=index_file,
        output_file=output_file,
        processes=2,
    )

    assert result.feature_rows == 3
    assert output_file.read_text(encoding="utf-8") == (
        "accession\tmotif_id\n"
        "P1\tMOTIF_A\n"
        "P1\tMOTIF_B\n"
        "P2\tMOTIF_B\n"
        "# completed\ttrue\n"
    )


class _RecordingPattern:
    def __init__(self, *, match_start: int) -> None:
        self.match_start = match_start
        self.starts: list[int] = []

    def match(self, _sequence: str, start: int):
        self.starts.append(start)
        if start == self.match_start:
            return object()
        return None


def _write_indexed_fasta(
    fasta_file: Path,
    index_file: Path,
    sequences: dict[str, str],
) -> None:
    with fasta_file.open("wb") as fasta_handle:
        locations: dict[str, tuple[int, int]] = {}
        for accession, sequence in sequences.items():
            fasta_handle.write(f">{accession}\n".encode("ascii"))
            offset = fasta_handle.tell()
            fasta_handle.write(sequence.encode("ascii"))
            fasta_handle.write(b"\n")
            locations[accession] = (offset, len(sequence))
    fasta_stat = fasta_file.stat()
    with index_file.open("w", encoding="ascii") as index_handle:
        index_handle.write("# Generated by test\n")
        index_handle.write("# version\t1\n")
        index_handle.write(f"# fasta_name\t{fasta_file.name}\n")
        index_handle.write(f"# fasta_size\t{fasta_stat.st_size}\n")
        index_handle.write(f"# fasta_mtime_ns\t{fasta_stat.st_mtime_ns}\n")
        index_handle.write("accession\toffset\tlength\n")
        for accession, (offset, length) in locations.items():
            index_handle.write(f"{accession}\t{offset}\t{length}\n")
