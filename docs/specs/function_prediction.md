# Function Prediction Spec

## Goal

Predict protein function from discovered or curated sequence signatures.

## Inputs

- Protein sequences.
- Signature library.
- Signature-to-function annotations.
- Optional taxonomic, family, or domain metadata.

## Outputs

- Predicted functions.
- Supporting signature hits.
- Confidence or ranking scores.
- Explanation metadata suitable for inspection.

## Initial Workflow Notes

- Scan sequences against a signature library.
- Aggregate signature hits into function evidence.
- Rank candidate functions by evidence strength.
- Report predictions with traceable signature support.

## Open Questions

- Should prediction use rule-based scoring first, supervised learning first, or both?
- How should conflicting signatures be resolved?
- Which metrics define success for the first benchmark?
