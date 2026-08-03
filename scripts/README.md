# Scripts

Use this directory for small operational scripts that do not belong in the package API.

## Function label treemap

`plot_function_label_treemap.py` is intended to be run from `work/`. It reads
`clusters.tsv`, `accession_mf_go.tsv`, `clusters_meta.tsv`, and `go_graph.pkl`,
then writes two treemaps plus auditable TSVs:

- `go_label_population_treemap.html`: accessions grouped by their exact
  Molecular Function GO label set.
- `cluster_label_population_treemap.html`: accessions grouped by `cluster_id`,
  with function descriptions composed from each cluster's `composed_go`.

Install the plotting-only dependency outside the package metadata:

```bash
python -m pip install plotly
```

For static PNG/SVG/PDF export, also install:

```bash
python -m pip install kaleido
```

Example from the repository root:

```bash
cd work
python ../scripts/plot_function_label_treemap.py
```

Use `0` to plot all groups instead of only the largest defaults:

```bash
python ../scripts/plot_function_label_treemap.py \
  --top-go-groups 0 \
  --top-cluster-groups 0
```

## Cluster MDS

`plot_cluster_mds.py` reads a tab-separated cluster similarity matrix, converts
functional similarity to distance with `distance = 1 - similarity`, validates
the matrix, embeds clusters into 2D with metric MDS, and writes PNG/PDF plots
plus coordinates.

Install the internal plotting/analysis dependencies outside the package
metadata:

```bash
python -m pip install pandas scikit-learn
```

Install Seaborn as an optional styling dependency for cleaner Matplotlib output:

```bash
python -m pip install seaborn
```

Example:

```bash
python scripts/plot_cluster_mds.py \
  --input work/simulation_panel/cluster_similarity_matrix.tsv \
  --out-prefix work/simulation_panel/cluster_mds \
  --draw-lines
```

The plotted clusters are translucent gradient circles with random visual radii.
The default radius range is `0.05` to `0.25` in MDS coordinate units:

```bash
python scripts/plot_cluster_mds.py \
  --input work/simulation_panel/cluster_similarity_matrix.tsv \
  --out-prefix work/simulation_panel/cluster_mds \
  --min-radius 0.05 \
  --max-radius 0.25
```

## Scan Similarity Comparison

`compare_scan_similarity.py` randomly samples clusters and accessions from
`clusters.tsv`, computes GO set similarity across all selected accession pairs as
the background distribution, computes the same similarity for pairs within each
selected cluster, then scans each accession sequence with the same machinery used
by `prosig scan` and scores the top predicted GO set against the accession's
known Molecular Function GO profile. The random selection only requires valid GO
profiles and FASTA sequences; it does not apply GO-similarity or sequence-
similarity constraints.

Example from the repository root:

```bash
python scripts/compare_scan_similarity.py \
  --clusters work/clusters.tsv \
  --library-dir work \
  --fasta work/accession.fasta \
  --fasta-index work/accession.fasta.idx \
  --cluster-count 10 \
  --accessions-per-cluster 10 \
  --seed 1 \
  --out-dir work/random_cluster_scan_similarity
```

The output directory contains:

- `selected_accessions.tsv`: sampled clusters/accessions and true GO profiles.
- `pairwise_similarity.tsv`: background and in-cluster pairwise GO similarities.
- `prediction_similarity.tsv`: one top scan prediction and similarity per accession.
- `similarity_summary.tsv`: summary statistics for the three distributions.
- `similarity_summary.json`: JSON copy of the same summary.

To reuse an existing panel instead of sampling from `clusters.tsv`, pass
`--selected-accessions path/to/selected_accessions.tsv`.
