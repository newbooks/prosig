# Feature–Cluster Enrichment: An Interpretable Evaluation for Feature Quality

## Purpose

This document describes a simple approach for testing whether engineered protein features carry information about functional class. The method asks one direct question: **when a feature is present, is it more common in one functional cluster than in the rest of the reference proteins?** The answer is represented as a feature–cluster enrichment weight. These weights provide an interpretable feature-quality report and can also support a lightweight cluster-level function predictor.

The method is deliberately simple. It does not assume a particular feature technology, model family, or biological mechanism. A feature may describe sequence composition, structure, physicochemical properties, domains, embeddings, conservation, localization, or any other measurable property. This method is most direct for binary features, but continuous and categorical features can be converted into reproducible binary states.

The binary feature is particular useful in detecting the presense a catalytic site or binding site which plays important roles in protein function.

The current evaluation panel contains 10 functional clusters with 10 sequence-diverse proteins per cluster. Functional labels are defined at the cluster level, while sequence diversity reduces the chance that a feature looks useful merely because the panel contains nearly identical proteins.

## Background and intuition

A useful feature for a functional cluster should satisfy two conditions:

1. **Sensitivity within the cluster:** it should occur in several proteins that have the target function.
2. **Specificity outside the cluster:** it should be uncommon in proteins from other functional clusters.

Neither condition is sufficient alone. A feature found in every protein in a cluster but also in every background protein is not discriminative. A feature found in only one target protein and nowhere else may have high apparent specificity, but the evidence is too sparse to trust. The enrichment baseline therefore records both the frequency contrast and the number of supporting examples.

The result is a scoreboard indexed by feature and functional cluster:

```text
Feature F1 -> cluster_01: strong positive evidence
Feature F1 -> cluster_02: weak positive evidence
Feature F1 -> cluster_03: no positive association
```

A feature has no universal weight. Its meaning is always conditional on the cluster being evaluated. In another word, a weight is given to a feature and cluster combination.

## Required inputs

### Cluster membership

Each reference protein must have one known cluster label:

```text
protein_id  cluster_id
P00001      cluster_01
P00002      cluster_01
P00003      cluster_02
```

The cluster assignment is treated as the target label. For a fair comparison, the clusters should be big enough so the frequency calculation won't swing widely because of the small denominator.

### Feature observations

The simplest input is a sparse table of positive feature observations:

```text
protein_id  feature_id
P00001      feature_A
P00001      feature_C
P00003      feature_B
```

Each protein–feature pair is binary: present or absent. Duplicate observations must be collapsed so that repeated detections do not increase support.

Feature values can be adapted as follows:

* **Binary features:** use presence directly.
* **Continuous features:** convert to binary by threshold.

## Mechanism

### 1. Evaluate every feature against every cluster

For a feature `F` and target cluster `C`, treat membership in `C` as the positive class and all other clusters as the background. Count a two-by-two contingency table:

|                 | In cluster C | Outside cluster C |
| --------------- | ------------ | ----------------- |
| Feature present | TP           | FP                |
| Feature absent  | FN           | TN                |

The counts mean:

* `TP`: target-cluster proteins containing the feature;
* `FP`: background proteins containing the feature;
* `FN`: target-cluster proteins without the feature;
* `TN`: background proteins without the feature.

Also record:

```text
support = TP
cluster_size = TP + FN
background_size = FP + TN
```

### 2. Estimate smoothed frequencies

Raw frequencies can produce an infinite enrichment when no background protein contains a feature. Apply a Jeffreys-prior pseudocount of `0.5`:

```text
cluster_frequency = (TP + 0.5) / (TP + FN + 1)
background_frequency = (FP + 0.5) / (FP + TN + 1)
```

This smoothing keeps all weights finite and reduces the apparent strength of associations supported by very few proteins. It does not eliminate small-sample uncertainty, so support must still be reported separately.

### 3. Calculate feature–cluster enrichment

Define the enrichment ratio and its base-2 logarithm:

```text
enrichment = cluster_frequency / background_frequency
weight = log2(enrichment)
```

Interpretation:

| Weight | Interpretation                                 |
| ------ | ---------------------------------------------- |
| `0`    | Equal frequency inside and outside the cluster |
| `1`    | Two-fold enrichment in the cluster             |
| `2`    | Four-fold enrichment                           |
| `3`    | Eight-fold enrichment                          |
| `< 0`  | Feature is depleted in the cluster             |

The weight measures association strength, not confidence. For example, a weight of 4 supported by one protein is much less convincing than a weight of 4 supported by eight proteins.

### 4. Apply minimum-evidence rules

The original production baseline ignores clusters with fewer than 10 members and feature–cluster associations with support below 5. For the balanced 10-protein panel, `support >= 5` means that a retained feature must occur in at least half of the target cluster.

Recommended first-pass rules are:

```text
minimum cluster size = 10
minimum support = 5
store positive weights only for prediction
```

For feature research, retain all weights in the diagnostic report, including zero and negative values. Depletion can reveal redundancy, anticorrelation, or potential negative evidence. A compact prediction artifact may keep only positive associations after the analysis is complete.

## Simplified example

Assume three clusters contain 10 proteins each. Feature `F1` occurs in 8 of the 10 proteins in `cluster_01`, and in 1 of the 20 proteins outside that cluster:

```text
TP = 8
FN = 2
FP = 1
TN = 19
support = 8
```

The smoothed frequencies are:

```text
cluster_frequency    = (8 + 0.5) / (10 + 1) = 0.7727
background_frequency = (1 + 0.5) / (20 + 1) = 0.0714
```

Therefore:

```text
enrichment = 0.7727 / 0.0714 = 10.82
weight = log2(10.82) = 3.44
```

Feature `F1` is about 10.8 times as frequent in `cluster_01` as in the background, giving strong positive evidence for that cluster. Its support of 8 also passes the proposed minimum-support rule.

The same feature must be evaluated separately against `cluster_02` and `cluster_03`. It may have a positive, neutral, or negative weight for each.

### Quality summaries for each engineered feature

At minimum, report:

* the cluster with the highest positive weight;
* weight and support for that cluster;
* target and background frequencies;
* the number of clusters receiving positive weight;
* coverage: fraction of proteins on which the feature is present;
* out-of-fold top-1 and top-3 cluster accuracy when used for prediction.

A high-quality cluster-specific feature has strong enrichment, adequate support, and few competing cluster associations. A broadly present feature may have good coverage but poor specificity. A rare feature may have excellent specificity but insufficient coverage.

## Cluster-level function prediction

Once feature–cluster weights have been learned from reference proteins, a query protein can be scored using its observed features:

1. Extract the query's features using exactly the same preprocessing and thresholds used for the reference panel.
2. Retrieve every cluster association for each present feature.
3. Discard associations below a chosen weight or support threshold.
4. Aggregate evidence by cluster.
5. Rank clusters and report the contributing features with their weights and support.

The simplest aggregation rule is the maximum supporting weight:

```text
query_score(cluster) = max(weight(feature, cluster) for present features)
```

This rule is easy to interpret and prevents many correlated features from artificially inflating a score. A prediction can be explained as “cluster_04 was ranked first because feature_F7 had weight 3.6 with support 7.”

Alternative aggregations—sum, capped sum, top-k mean, or a learned linear combination—may improve accuracy, but correlated features must then be handled carefully. The maximum-weight rule should be retained as the baseline against which more complex models are compared.

Prediction output should include the ranked cluster IDs, cluster-level function descriptions, score, supporting feature IDs, support counts, and the calibration performance associated with the applied threshold. The output should describe cluster membership evidence rather than claim a definitive molecular function.

## Strengths

* **Interpretable:** every score traces to counts and a frequency ratio.
* **Feature agnostic:** any reproducible feature can be represented and tested.
* **Simple to implement:** no optimizer or large training dataset is required.
* **Support aware:** rare observations can be filtered or clearly flagged.
* **Sparse and efficient:** only observed features and useful associations need to be stored for prediction.
* **Auditable:** predictions can name the exact features that supplied evidence.

## Limitations

* **Small-sample instability:** with 10 proteins per cluster, one protein can materially change frequency and weight.
* **Dependence on cluster labels:** noisy or overly broad functional clusters directly affect the learned associations.
* **Binary simplification:** thresholding discards magnitude, location, count, and uncertainty information.
* **Association is not causation:** enrichment does not prove that a feature is mechanistically responsible for the function.

This procedure turns engineered feature evaluation into a common, reproducible experiment: a feature is valuable when it is repeatedly present in a functional cluster, uncommon in the appropriate background, supported by enough proteins, and predictive for proteins that were not used to define or weight it.
