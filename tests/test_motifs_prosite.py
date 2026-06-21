import logging
from pathlib import Path

from prosig.motifs.prosite import (
    read_prosite_patterns,
    translate_prosite_pattern,
    write_prosig_motif_library,
)


def test_translate_prosite_pattern_to_prosig_syntax() -> None:
    assert translate_prosite_pattern("C-x(2)-C-x(10,15)-H-x(2)-H.") == (
        "C?(2)C?(10,15)H?(2)H"
    )
    assert translate_prosite_pattern("N-{P}-[ST]-{P}.") == "N!P[ST]!P"
    assert translate_prosite_pattern("B-[DE]-N-{P}-[ST]-{P}.") == (
        "[DN]{-}N!P[ST]!P"
    )
    assert translate_prosite_pattern("F-[IVFY]-G-[LM]-M-[G>].") == (
        "F[IVFY]G[LM]M[G>]"
    )


def test_read_prosite_patterns_selects_pattern_entries_with_pa(tmp_path: Path) -> None:
    prosite_file = tmp_path / "prosite.dat"
    prosite_file.write_text(
        "\n".join(
            [
                "ID   KEEP_ME; PATTERN.",
                "AC   PS00001;",
                "DE   Valid pattern.",
                "PA   N-{P}-[ST]-{P}.",
                "//",
                "ID   SKIP_NO_PA; PATTERN.",
                "AC   PS00002;",
                "DE   Missing PA.",
                "//",
                "ID   SKIP_PROFILE; MATRIX.",
                "AC   PS00003;",
                "DE   Not a PATTERN entry.",
                "PA   A.",
                "//",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    result = read_prosite_patterns(prosite_file)

    assert result.total_entries == 3
    assert result.pattern_entries == 2
    assert result.skipped_pattern_entries_without_pa == 1
    assert [entry.name for entry in result.entries] == ["KEEP_ME"]
    assert result.entries[0].description == "Valid pattern"
    assert result.entries[0].prosite_pattern == "N-{P}-[ST]-{P}"


def test_write_prosig_motif_library_omits_failed_translations(
    tmp_path: Path,
    caplog,
) -> None:
    prosite_file = tmp_path / "prosite.dat"
    output_file = tmp_path / "prosig_motifs.tsv"
    prosite_file.write_text(
        "\n".join(
            [
                "ID   CAMP_PHOSPHO_SITE; PATTERN.",
                "AC   PS00004;",
                "DE   cAMP- and cGMP-dependent protein kinase phosphorylation site.",
                "PA   B-[DE]-N-{P}-[ST]-{P}.",
                "//",
                "ID   UNSUPPORTED; PATTERN.",
                "AC   PS99999;",
                "DE   Exclusion containing any residue cannot be represented.",
                "PA   A-{X}-G.",
                "//",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    caplog.set_level(logging.WARNING)
    result = write_prosig_motif_library(prosite_file, output_file)

    assert result.stats.total_entries == 2
    assert result.stats.pattern_entries == 2
    assert result.stats.translated_entries == 1
    assert result.stats.failed_entries == 1
    assert result.stats.macro_converted_entries == 1
    assert result.stats.ambiguous_codes_translated == 1
    assert "Skipping unsupported PROSITE motif UNSUPPORTED" in caplog.text
    assert output_file.read_text(encoding="utf-8") == (
        "# ProSig motif library\n"
        "name\tprosite_ac\tdescription\tprosite_pattern\tprosig_pattern\tstatus\n"
        "CAMP_PHOSPHO_SITE\tPS00004\t"
        "cAMP- and cGMP-dependent protein kinase phosphorylation site\t"
        "B-[DE]-N-{P}-[ST]-{P}\t[DN]{-}N!P[ST]!P\tprosite\n"
    )
