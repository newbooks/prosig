# Motif Discovery Spec

## Goal

Discover protein sequence signatures that distinguish a positive sequence set from a background or negative sequence set.

## Inputs

- Positive protein sequences.
- Background or negative protein sequences.
- Optional labels, metadata, or family annotations.

## Outputs

- Discovered signatures or motifs.
- Scores and enrichment statistics.
- Sequence hits and coordinates.
- Optional generalized motif representation or PWM.

## Initial Algorithm Notes

- Start with k-mer enumeration.
- Score candidate enrichment with an exact or appropriate statistical test.
- Generalize enriched k-mers into motif/signature patterns.
- Optionally refine candidates with PWM-like representations.
- Compare design and behavior with motif handling methods in `../pclass`.

## Open Questions

- What motif representation should be canonical in ProSig?
- How should ambiguous amino acid symbols be handled?
- What background model is required for the first milestone?
- What thresholds should be defaults versus user-configurable parameters?
