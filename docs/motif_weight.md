# Motif Weight Calculation for Functional Cluster Prediction

## Purpose

This document describes the baseline algorithm for quantifying how informative a sequence motif is for predicting a functional cluster.

The objective is **not** to determine whether a motif is biologically required for a function. Instead, the objective is to measure whether observing a motif significantly increases the probability that a protein belongs to a particular functional cluster.

The resulting **motif-cluster weight** is later used as evidence during function prediction.

---

# Definitions

## Functional Cluster

A functional cluster is a collection of proteins grouped by functional similarity (derived from GO similarity).

Every protein belongs to exactly one functional cluster.

---

## Motif Match

Each motif is scanned against every protein sequence.

For the baseline implementation, motif matching is binary.

* **1** = motif found
* **0** = motif not found

Multiple occurrences, motif position, and motif density are ignored in Version 1.

---

# Motif-Cluster Association

Each motif is evaluated **independently** against **every functional cluster**.

A motif has **no intrinsic predictive weight**.

Instead, every **motif-cluster pair** has its own predictive statistics and weight.

Conceptually,

```
Motif M
    ↓
Cluster A : weight = 4.2
Cluster B : weight = 1.7
Cluster C : no association
```

The same motif may strongly predict one functional cluster while providing little evidence for another.

---

# Contingency Table

For every motif **M** and functional cluster **C**, construct the following contingency table.

|               | Protein in Cluster C | Protein outside Cluster C |
| ------------- | -------------------: | ------------------------: |
| Motif present |                   TP |                        FP |
| Motif absent  |                   FN |                        TN |

where

* **TP** — proteins that belong to Cluster C and contain motif M.
* **FP** — proteins outside Cluster C that contain motif M.
* **FN** — proteins inside Cluster C that do not contain motif M.
* **TN** — proteins outside Cluster C that do not contain motif M.

These values are obtained directly by counting proteins in the dataset.

No prior biological knowledge is assumed.

---

# Derived Statistics

## Smoothed Cluster Frequency

Fraction of proteins in the cluster containing the motif.

```text
cluster_frequency = (TP + 0.5) / (TP + FN + 1)
```

---

## Smoothed Background Frequency

Fraction of proteins outside the cluster containing the motif.

```text
background_frequency = (FP + 0.5) / (FP + TN + 1)
```

---

## Enrichment

Measure how much more frequently the motif occurs inside the cluster than outside.

```text
cluster_frequency = (TP + 0.5) / (TP + FN + 1)
background_frequency = (FP + 0.5) / (FP + TN + 1)
enrichment = cluster_frequency / background_frequency
```

The `0.5` pseudocount is a Jeffreys prior for binary motif presence. This keeps
weights finite when `FP = 0` and makes low-support perfect-background cases less
dominant than high-support cases.

Interpretation

* **1** : motif occurs equally often inside and outside the cluster.
* **>1** : motif is enriched in the cluster.
* **<1** : motif is depleted in the cluster.

---

## Log Enrichment (Motif Weight)

Because enrichment may span several orders of magnitude, the logarithm is used as the motif weight.

```text
weight = log2(enrichment)
```

Interpretation

| Weight | Meaning                                     |
| ------ | ------------------------------------------- |
| 0      | No predictive value                         |
| 1      | Motif is twice as common inside the cluster |
| 2      | Four times as common                        |
| 3      | Eight times as common                       |
| 5      | Thirty-two times as common                  |

Negative values indicate that the motif is less common than expected.

---

# Support

Highly enriched motifs occurring only once or twice are unreliable.

Therefore record

```text
support = TP
```

Support measures the amount of evidence behind the association.

Version 1 discards motif-cluster pairs with insufficient support:

```text
TP < 5
```

This prevents sparse motif hits from becoming high-weight but unreliable
prediction evidence.

---

# Cluster Size Filter

Very small functional clusters do not provide enough positive examples for
stable motif-cluster weights.

Version 1 discards all motif-cluster combinations for clusters with:

```text
cluster_size < 10
```

These combinations are counted in the score board metadata as ignored because
of cluster size.

---

# Stored Statistics

Each **motif-cluster pair** stores the following information.

| Field                | Description                                    |
| -------------------- | ---------------------------------------------- |
| motif_id             | Motif identifier                               |
| cluster_id           | Functional cluster identifier                  |
| TP                   | Proteins in cluster with motif                 |
| FP                   | Proteins outside cluster with motif            |
| FN                   | Proteins in cluster without motif              |
| TN                   | Proteins outside cluster without motif         |
| support              | TP                                             |
| cluster_frequency    | (TP + 0.5) / (TP + FN + 1)                     |
| background_frequency | (FP + 0.5) / (FP + TN + 1)                     |
| weight               | log₂(cluster_frequency / background_frequency) |

This table represents the learned association between motifs and functional
clusters.

For the production score board artifact, only motif-cluster pairs with
positive predictive weight are stored:

```text
weight > 0
```

Zero and negative weights are omitted to keep the machine-readable artifact
compact. Omitted combinations are counted in metadata as non-positive weights.

---

# Artifact Format

The score board is written as a pickle file for machine consumption:

```text
motif_cluster_scoreboard.pkl
```

The pickle contains a dictionary with:

* `schema_version`
* `kind`
* `parameters`
* `weights`

`weights` is nested by motif ID and then cluster ID:

```python
{
    "motif_A": {
        "cluster_0001": {
            "TP": 12,
            "FP": 3,
            "FN": 8,
            "TN": 977,
            "support": 12,
            "cluster_frequency": 0.6,
            "background_frequency": 0.003061,
            "weight": 7.615,
        }
    }
}
```

A companion JSON metadata file is written for diagnostics:

```text
motif_cluster_scoreboard_meta.json
```

Metadata includes counts for:

* total motif-cluster combinations;
* combinations ignored because the cluster has fewer than 10 members;
* combinations ignored because support is below 5;
* combinations ignored because background frequency is unavailable;
* combinations ignored because weight is zero or negative;
* total non-zero positive weights stored;
* internal calibration top-1, top-3, set accuracy, average prediction count, and
  coverage at motif-weight thresholds 2.0 through 8.0 in 0.5 increments.

---

# Function Prediction

Given a query protein,

1. Scan the sequence against the motif library.
2. Identify all matched motifs.
3. For each matched motif:

   * Retrieve all associated motif-cluster records.
   * Retrieve the stored weight.
4. Accumulate weights for every functional cluster.

Example

```
Motif A

Cluster 15 : +3.4
Cluster 22 : +1.1

Motif B

Cluster 15 : +2.2
Cluster 81 : +4.0

Motif C

Cluster 15 : +1.6
```

Accumulated scores

```
Cluster 15 = 7.2
Cluster 81 = 4.0
Cluster 22 = 1.1
```

Clusters are ranked by total score.

The highest-scoring clusters become the predicted functions.

---

# Design Rationale

This baseline intentionally uses a simple and interpretable scoring model.

Advantages

* Completely data-driven.
* No manually assigned motif weights.
* Independent of motif type.
* Independent of clustering algorithm.
* Easily interpretable.
* Computationally inexpensive.
* Naturally extensible to larger motif libraries.

Most importantly, every weight has a direct probabilistic interpretation:

> **How much more frequently is this motif observed inside a functional cluster than outside?**

This makes motif evidence additive and suitable for weighted voting during function prediction.

---

# Future Extensions

This baseline serves as a reference implementation.

Future versions may incorporate

* Multiple motif occurrences.
* First and last motif positions.
* Motif density.
* Motif co-occurrence.
* Statistical significance testing.
* Bayesian weighting.
* Machine learning refinement.

The baseline motif-cluster weight should remain as the reference method against which more sophisticated scoring algorithms are evaluated.
