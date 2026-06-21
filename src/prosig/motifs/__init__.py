"""Motif library construction and translation helpers."""

from prosig.motifs.prosite import (
    PROSIG_MOTIF_HEADER,
    ProSigMotifRow,
    PrositeMotifTranslationResult,
    PrositeMotifTranslationStats,
    PrositePatternEntry,
    PrositeReadResult,
    translate_prosite_pattern,
    write_prosig_motif_library,
)

__all__ = [
    "PROSIG_MOTIF_HEADER",
    "ProSigMotifRow",
    "PrositeMotifTranslationResult",
    "PrositeMotifTranslationStats",
    "PrositePatternEntry",
    "PrositeReadResult",
    "translate_prosite_pattern",
    "write_prosig_motif_library",
]
