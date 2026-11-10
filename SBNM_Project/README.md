# LIONESS: Sample-Specific Gene Co-expression Networks in Pancreatic Cancer

**Course:** CB 5313 - Systems Biology & Network Modeling

## 1. Research question

An ordinary co-expression network is one network for a whole group of samples.
**LIONESS** (Kuijjer et al., 2019) estimates one network *per sample* with
`w_q = N * w_all - (N - 1) * w_(all without q)`, where `w` is any network
method (here: Pearson correlation) and `N` is the number of samples.

We ask:

1. Do pancreatic tumours show different gene-gene co-expression ("wiring") from
   the matched normal pancreas of the same patients, and does LIONESS detect it?
2. How reliable are LIONESS results when the number of patients is small, and
   how do they compare with the classical group-wise differential co-expression
   method?

## 2. Data

* **GEO series GSE15471** (Badea et al., 2008): whole-tissue expression of
  pancreatic ductal adenocarcinoma, Affymetrix Human U133 Plus 2.0 (GPL570).
* 36 patients, each with **one tumour and one matched normal sample** = 72
  samples (6 technical replicate arrays in the original series are removed).
* The series matrix is already RMA-normalised log2 expression, downloaded
  directly from NCBI GEO by `01_download_data.py`. Probes are mapped to gene
  symbols with the GPL570 annotation (ambiguous probes dropped, best probe per
  gene = highest mean expression) -> 20,848 genes.
* For LIONESS we use the **1000 most variable genes** (499,500 gene pairs).

No raw data is stored in the repository; the script downloads ~32 MB from GEO.

## 3. How to run

Tested with Python 3.9. Needs an internet connection for step 2.

```bash
git clone <your-repo-url>
cd <repo-folder>
pip install -r requirements.txt
python 01_download_data.py                # download + preprocess (about 1 minute)
python 02_run_lioness.py                  # LIONESS networks, ~20 s, writes 144 MB to results/
python 03_analyze_results.py              # main analysis + figures
python 04_parameter_sweep_and_comparison.py   # sweeps + method comparison (~40 s)
python -m pytest tests                    # unit tests for the LIONESS implementation
```

All random steps use fixed seeds, and a fresh run in an empty folder gave
exactly the same numbers as those in `results/` (checked with `diff`).

| File | Purpose |
|---|---|
| `lioness_utils.py` | LIONESS, paired test, FDR, shared helpers |
| `01_download_data.py` | download GEO data, probes -> genes |
| `02_run_lioness.py` | sample-specific networks (`--n-genes` option) |
| `03_analyze_results.py` | differential edges, hubs, sample separation |
| `04_parameter_sweep_and_comparison.py` | **modifications beyond the tutorial** |
| `tests/test_lioness.py` | checks the code against the LIONESS definition |
| `results/` | figures (`*.png`), tables (`*.csv`), text summaries (`*.txt`) |

## 4. What was changed beyond the tutorial

1. **Different, real dataset** (the previous version used synthetic random
   numbers, so nothing real could be found).
2. **Parameter sweeps**: number of genes (250-1500), cohort size (9/18/27/36
   patients) and edge-selection cut-off (top 0.1 / 1 / 5 % of edges).
3. **Method comparison**: LIONESS (+ paired test across patients) versus the
   classical differential co-expression `r_tumour - r_normal` (Fisher z), using
   the same genes and the same "top-k edges" and stability criteria.

Evaluation metrics: number of significant edges (Benjamini-Hochberg FDR 5 %),
a label-swap permutation null (empirical FDR), Jaccard stability of the top
edges across random patient subsets, and patient-grouped cross-validated
classification accuracy.

## 5. Results and interpretation

All numbers below come from `results/summary_main_analysis.txt` and
`results/summary_sweep_comparison.txt`.

### 5.1 Tumour vs normal edges (1000 genes)

* 39,074 of 499,500 edges (7.8 %) differ significantly between matched tumour
  and normal samples (paired t-test, FDR 5 %). 15,377 are stronger in tumour and
  23,697 weaker (`fig1`).
* Label-swap null: on average 1,142 edges pass the same cut-off by chance
  (empirical FDR 0.029). But the null is very spread out (median 464, 95th
  percentile 4,428; 95th-percentile FDR = 0.11), because edges that share a gene
  are strongly dependent. **The p-values should be read as a ranking, not as
  exact error rates.**

### 5.2 Is this just differential expression?

* 811 of the 1000 genes are themselves differentially expressed. This is not a
  surprise: we picked genes by variance over all samples, and tumour-vs-normal
  is the biggest source of variance. It means this gene set is biased towards
  disease genes.
* Significant edges are not more likely than other edges to join two DE genes
  (64.2 % vs 65.8 %), and the edge change is only weakly (and negatively)
  related to the product of the DE statistics (Spearman rho = -0.11, `fig2`).
  So there is **no simple "edge change = expression change" relationship**.
  A likely reason is that a LIONESS weight measures how much a sample
  *deviates from the shared trend* of the two genes, not just whether both are
  high. We did not test this explanation further. Because almost every gene is
  DE, this check has low power, so we do not claim that the edges are
  independent of expression.
* Top hub genes of the differential network include *MUC1*, *IGFBP3*, *CP*,
  *MMP12*, *LYVE1* (`results/top_hub_genes.csv`); several are known to be
  changed in pancreatic cancer, but we did **not** run an enrichment analysis,
  so this is only an observation. Three of the top ten hubs (*CLC*, *CCL8*,
  *INS*) are not DE; *INS* (insulin) suggests that the mix of cell types in the
  normal tissue (islets) contributes to the network differences. This is a
  hypothesis, not something we tested.

### 5.3 Do networks separate the samples?

* PCA (`fig3`): tumour networks form a tight cluster, normal networks are much
  more spread out (a few normal samples are extreme). Median mean-|edge weight|
  is 0.83 for normal and 0.61 for tumour samples.
* Patient-grouped 6-fold cross-validation: **expression 0.889, LIONESS edges
  0.764**, shuffled-label control 0.403. The networks do carry the
  tumour/normal signal, but in this dataset they are **less informative than
  the expression values they were built from**. (Classifier settings were not
  tuned; the LIONESS networks are computed from all samples without using
  labels, which is allowed but means it is a transductive, not a fully
  independent, test.)

### 5.4 Parameter sweeps

* **Number of genes** (`sweep_n_genes.csv`, `fig4` left): 6.5 % (250 genes),
  9.1 % (500), 7.8 % (1000), 7.0 % (1500) of edges are significant; empirical
  FDR about 0.03-0.04 each time. Important point: with Pearson LIONESS the
  weight of a gene pair depends only on those two genes, so this sweep changes
  only *which and how many edges are tested*, not the edge values. The result
  is therefore stable but not very informative.
* **Cohort size** (`stability_subsampling.csv`, `fig4` right): with 9 patients
  almost nothing is significant (1 edge on average); 18 patients give 321 ± 730
  edges (very variable); 27 give 3,894 ± 2,503; all 36 give 11,326 (500 genes).
  The overlap of the top 1 % of edges with the full-cohort top 1 % rises from
  0.07 (9 patients) to 0.17 (18) to 0.35 (27). Note that the 27-patient subsets
  share 75 % of their data with the reference, so the real agreement between
  independent cohorts would be lower. **Conclusion: with a few dozen samples,
  individual LIONESS edges are noisy; only the large-scale picture is reliable.**
* **Sparsity cut-off**: the top 5 % of edges are more stable (Jaccard 0.16-0.49)
  than the top 0.1 % (0.02-0.22). Stricter cut-offs give shorter edge lists but
  not more reproducible ones - the opposite of what the first version claimed.

### 5.5 LIONESS vs classical differential co-expression

* On the same 500 genes the two methods agree only moderately: ranks of |statistic|
  correlate (Spearman 0.69) but the top 0.1 % / 1 % / 5 % lists overlap with
  Jaccard 0.004 / 0.10 / 0.24 (`fig5`).
* The classical method was **more stable** under subsampling (top-1 % Jaccard
  0.13 / 0.27 / 0.44 vs 0.07 / 0.17 / 0.35 for 9 / 18 / 27 patients).
* The classical method calls many more edges significant (133,745 of 499,500
  at 1000 genes, empirical FDR 0.066) because it ignores that tumour and normal
  samples come from the same patient; the LIONESS test uses the pairing.

**Take-home message.** LIONESS gives a sample-level network for each patient
(useful for questions about individual heterogeneity), and its group-level
differences are detectable and show tumour/normal structure. But for the
simple question "which edges differ between two groups" it was *not* better
than the classical method in this dataset: it was less stable with a small
cohort. This is a negative result for LIONESS as a differential-network tool,
not a flaw in the code (the implementation is tested against the definition in
`tests/test_lioness.py`).

## 6. Limitations and bias

* **Stability is not accuracy.** There is no gold-standard set of true edges,
  so stability, permutation null and classification are proxies. We cannot say
  which method finds the "true" rewiring.
* Gene selection by variance over all samples favours DE genes (811/1000).
* Only one dataset and one base network (Pearson), so conclusions may not
  generalise. Pearson LIONESS weights are not bounded by [-1, 1] (range about
  -30 to 37 here) and are sensitive to outlying samples.
* Edges share genes and LIONESS networks share data, so the statistical tests
  are only approximately valid (see the wide permutation null).
* Whole-tissue microarray data mix cell types; differences may reflect tissue
  composition (e.g. fewer acinar/islet cells in tumour tissue), not regulation.
  Co-expression is not proof of regulatory interaction.
* Subsampling reference sets overlap with the subsets (see 5.4).

## References

* Kuijjer ML, Tung MG, Yuan G, Quackenbush J, Glass K. Estimating
  sample-specific regulatory networks. *iScience* 2019;14:226-240.
* Badea L, Herlea V, Dima SO, Dumitrascu T, Popescu I. Combined gene expression
  analysis of whole-tissue and microdissected pancreatic ductal adenocarcinoma.
  *Hepatogastroenterology* 2008;55:2016-2027. (GEO: GSE15471)
