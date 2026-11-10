# LIONESS: Sample-Specific Gene Co-expression Networks in Pancreatic Cancer

**Course:** CB 5313 - Systems Biology & Network Modeling

## 1. Research Question

Traditional gene co-expression networks represent relationships between genes across an entire group of samples. However, they do not show how these relationships may vary from one individual sample to another.

To address this, we implemented **LIONESS (Linear Interpolation to Obtain Network Estimates for Single Samples)**, based on the method proposed by Kuijjer et al. (2019). LIONESS estimates an individual gene co-expression network for every sample instead of generating only one network for the complete dataset.

The method calculates a sample-specific network using the following equation:

$$
w_q = N \times w_{\mathrm{all}} - (N-1) \times w_{\mathrm{all\ without\ }q}
$$

Here, $w$ represents the network weight calculated using a network inference method (Pearson correlation in our implementation), and $N$ is the total number of samples.

Through this project, we investigated two main questions:

1. Do pancreatic tumour samples show differences in gene-gene co-expression patterns (network wiring) compared with the matched normal pancreas samples from the same patients? Can LIONESS capture these differences?
2. How reliable are LIONESS networks when the number of patients is relatively small, and how do they compare with the traditional group-wise differential co-expression method?

## 2. Dataset and Preprocessing

For this project, we used the publicly available **GEO dataset GSE15471** from Badea et al. (2008). It contains gene expression data from pancreatic ductal adenocarcinoma (PDAC) tissue, generated using the Affymetrix Human Genome U133 Plus 2.0 microarray platform (GPL570).

The dataset contains:

- 36 patients, with one tumour sample and one matched normal sample for each patient.
- A total of 72 samples after removing 6 technical replicate arrays from the original dataset.
- RMA-normalised log2 gene expression values obtained directly from the GEO series matrix.
- Gene annotation performed using the GPL570 platform information.

During preprocessing, probes were mapped to their corresponding gene symbols. Ambiguous probes were removed, and when multiple probes represented the same gene, we retained the probe with the highest mean expression. This resulted in a final dataset containing **20,848 genes**.

For the LIONESS analysis, we selected the **1,000 most variable genes**, resulting in 499,500 possible gene-gene pairs.

The raw dataset is not included in the repository. Instead, the download and preprocessing script automatically retrieves the required data from NCBI GEO (approximately 32 MB) when the pipeline is executed.

## 3. How to Run the Project

The project was developed and tested using Python 3.9. An active internet connection is required for downloading the dataset during the preprocessing stage.

Follow these steps to run the complete pipeline:

```bash
git clone https://github.com/SaniyaArshadAli/Systems_Biology-Network_Modeling

cd <repo-folder> (E.g. cd C:\Users\HP OMNIBOOK\Downloads\SBNM_Project\SBNM_Project\SBNM_Project)

pip install -r requirements.txt

# Step 1: Download and preprocess the GEO dataset
python 01_download_data.py

# Step 2: Generate sample-specific LIONESS networks
python 02_run_lioness.py

# Step 3: Perform the main analysis and generate figures
python 03_analyze_results.py

# Step 4: Run parameter sweeps and compare methods
python 04_parameter_sweep_and_comparison.py

# Step 5: Run the unit tests
python -m pytest tests
```

The approximate execution times were:

- Data download and preprocessing: 1 minute.
- LIONESS network generation: 20 seconds.
- Parameter sweeps and method comparison: 40 seconds.

The LIONESS network generation step produces approximately 144 MB of output in the `results/` directory.

To make the results reproducible, all steps involving randomness use fixed random seeds. We also tested the pipeline by running it again in a fresh, empty folder and comparing the generated outputs with the existing results using `diff`. The numerical results matched exactly.

### Project File Structure

| File / Directory | Description |
|---|---|
| `lioness_utils.py` | Contains the LIONESS implementation, paired statistical testing, FDR correction, and shared utility functions. |
| `01_download_data.py` | Downloads the GEO dataset and performs probe-to-gene mapping and preprocessing. |
| `02_run_lioness.py` | Generates sample-specific networks. Also supports the `--n-genes` option. |
| `03_analyze_results.py` | Performs differential edge analysis, hub gene identification, and sample separation analysis. |
| `04_parameter_sweep_and_comparison.py` | Contains the additional experiments, parameter sweeps, and comparison with the classical method. |
| `tests/test_lioness.py` | Tests whether the LIONESS implementation follows the original mathematical definition. |
| `results/` | Stores the generated figures (`.png`), result tables (`.csv`), and analysis summaries (`.txt`). |

## 4. Changes and Experiments Beyond the Tutorial

Apart from implementing the basic LIONESS method, we made additional changes to explore how it performs under different conditions.

### 4.1 Using a Real Gene Expression Dataset

The earlier version used synthetic random numbers, which made it difficult to draw meaningful biological conclusions. We replaced these with the real pancreatic cancer dataset GSE15471, allowing us to study actual tumour and normal tissue samples.

### 4.2 Performing Parameter Sweeps

We experimented with different settings to understand how they affect the results. The parameters included:

- **Number of genes:** 250, 500, 1,000, and 1,500.
- **Cohort size:** 9, 18, 27, and 36 patients.
- **Edge-selection thresholds:** Top 0.1%, 1%, and 5% of all edges.

### 4.3 Comparing LIONESS with Classical Differential Co-expression

We also compared LIONESS, using paired statistical testing across patients, with the traditional differential co-expression approach.

The classical method calculates the difference between tumour and normal correlation values:

$$
r_{\mathrm{tumour}} - r_{\mathrm{normal}}
$$

We used Fisher's z-transformation for the correlation comparison.

To keep the comparison consistent, both methods were evaluated using the same set of genes, top-k edge-selection criteria, and stability measurements.

### Evaluation Metrics

We used several metrics to understand the results:

- **Number of significant edges:** Edges passing the Benjamini-Hochberg false discovery rate (FDR) threshold of 5%.
- **Label-swap permutation test:** Used to estimate how many edges could appear significant by chance and calculate empirical FDR.
- **Jaccard stability:** Measures how consistently the top-ranked edges appear across different random subsets of patients.
- **Patient-grouped cross-validation accuracy:** Evaluates how well the network features and gene expression values distinguish tumour samples from normal samples.

## 5. Results and Interpretation

The results discussed in this section are taken from `results/summary_main_analysis.txt` and `results/summary_sweep_comparison.txt`.

### 5.1 Tumour vs Normal Co-expression Differences

We first examined the differences in gene-gene co-expression between tumour and normal samples using the 1,000 selected genes.

Out of the 499,500 possible gene pairs, **39,074 edges (7.8%) showed statistically significant differences** between the matched tumour and normal samples using a paired t-test with an FDR threshold of 5%.

Among these edges:

- 15,377 showed stronger co-expression in tumour samples.
- 23,697 showed weaker co-expression in tumour samples.

These results are presented in `fig1`.

We also performed a label-swap permutation test to check whether the observed differences could simply be due to chance.

On average, 1,142 edges passed the same significance threshold under the permutation null model. This gave an empirical FDR of approximately 0.029.

However, the number of significant edges varied considerably across permutations:

- Median: 464 edges.
- 95th percentile: 4,428 edges.
- 95th-percentile FDR: 0.11.

This variation is important because different edges are not completely independent of each other. Edges that share a common gene can be strongly related.

Therefore, the statistical p-values should mainly be interpreted as a way of ranking potentially important edges rather than as exact estimates of the error rate.

### 5.2 Are the Network Differences Simply Due to Differential Gene Expression?

One important question we investigated was whether the changes in co-expression networks were simply a reflection of changes in individual gene expression.

Out of the 1,000 selected genes, **811 were differentially expressed (DE)** between tumour and normal samples.

This was expected to some extent because we selected genes based on their variance across all samples. Since the tumour-versus-normal difference is the largest source of variation in this dataset, the selected genes are naturally biased towards disease-associated genes.

We then checked whether significant network edges were more likely to connect two differentially expressed genes.

The results showed:

- 64.2% of significant edges connected two DE genes.
- 65.8% of the remaining edges connected two DE genes.

We also observed only a weak negative relationship between the change in edge weights and the product of the corresponding differential expression statistics (Spearman correlation, rho = -0.11). These results are shown in `fig2`.

This suggests that there is no simple relationship in which a change in gene expression automatically leads to a corresponding change in co-expression.

One possible explanation is that a LIONESS edge weight represents how much the relationship between two genes in an individual sample differs from the overall shared trend. It is therefore not simply a measure of whether both genes have high or low expression.

However, we did not investigate this explanation further. Also, since 811 out of 1,000 genes were differentially expressed, this analysis has limited power to establish whether network changes are truly independent of expression changes. We therefore do not claim that the two are independent.

We also examined the hub genes in the differential network. Some of the top hub genes included:

- MUC1
- IGFBP3
- CP
- MMP12
- LYVE1

Several of these genes have previously been reported to show changes in pancreatic cancer. However, we did not perform a formal enrichment analysis, so this is only an observation based on the identified genes.

Interestingly, three of the top ten hub genes (CLC, CCL8, and INS) were not differentially expressed.

The presence of INS (insulin) is particularly interesting because it may indicate that differences in the cellular composition of normal pancreatic tissue, especially the presence of insulin-producing islet cells, contribute to the observed network differences.

This is a possible explanation, not a conclusion established by our analysis.

### 5.3 Can the Networks Distinguish Tumour and Normal Samples?

We next explored whether the sample-specific networks could separate tumour samples from normal samples.

Using PCA (shown in `fig3`), we observed that the tumour networks formed a relatively tight cluster, while the normal networks were more widely spread out. A few normal samples appeared as extreme observations.

We also compared the median mean absolute edge weights:

- Normal samples: 0.83.
- Tumour samples: 0.61.

This suggests that the overall distribution of network weights differs between the two groups.

To evaluate how well the networks could distinguish tumour from normal samples, we performed patient-grouped 6-fold cross-validation.

The classification results were:

| Features used | Classification accuracy |
|---|---:|
| Gene expression values | 0.889 |
| LIONESS network edges | 0.764 |
| Shuffled-label control | 0.403 |

The LIONESS network features performed better than the shuffled-label control, indicating that they contain information related to tumour and normal status.

However, their accuracy was lower than that obtained using the original gene expression values.

This means that, for this dataset, the sample-specific networks were useful for capturing tumour-versus-normal differences, but they were not as informative for classification as the expression values from which they were generated.

There are also some important considerations when interpreting these results. The classifier settings were not tuned, and the LIONESS networks were generated using all samples without using their class labels. Although this does not introduce label leakage through the use of class labels during network construction, it means that the evaluation is transductive rather than a completely independent test.

### 5.4 Parameter Sweep Results

We performed additional experiments to understand how changing the number of genes, the number of patients, and the edge-selection threshold affects the results.

#### A. Effect of the Number of Genes

We tested four different gene counts: 250, 500, 1,000, and 1,500.

The percentage of significant edges was:

| Number of genes | Significant edges (%) |
|---|---:|
| 250 | 6.5% |
| 500 | 9.1% |
| 1,000 | 7.8% |
| 1,500 | 7.0% |

The empirical FDR remained approximately between 0.03 and 0.04 across these settings.

One important observation is that, with Pearson correlation-based LIONESS, the weight of a gene pair depends only on those two genes. Therefore, changing the number of selected genes does not change the calculated weight of an existing edge. It only changes which edges are included in the analysis and how many are tested.

This explains why the results remained relatively stable across the different gene counts. However, it also means that this particular parameter sweep provides limited additional information about the behaviour of individual edge weights.

The results are available in `sweep_n_genes.csv` and `fig4`.

#### B. Effect of Cohort Size

We also investigated how the number of patients affects the reliability of the detected network differences.

The experiment was performed using subsets containing 9, 18, and 27 patients, along with the complete cohort of 36 patients.

The average number of significant edges was:

| Number of patients | Average significant edges |
|---|---:|
| 9 | 1 |
| 18 | 321 ± 730 |
| 27 | 3,894 ± 2,503 |
| 36 | 11,326 |

These results were obtained using 500 genes.

The variation was particularly high for the 18- and 27-patient subsets, showing that the number of detected significant edges can change considerably depending on which patients are included.

We also measured the overlap between the top 1% of edges from each subset and the top 1% from the full cohort.

| Number of patients | Jaccard overlap |
|---|---:|
| 9 | 0.07 |
| 18 | 0.17 |
| 27 | 0.35 |

As the cohort size increased, the overlap with the full-cohort results also increased.

However, there is an important limitation here. The 27-patient subsets share 75% of their samples with the full reference cohort. Therefore, this overlap is not equivalent to the agreement we would expect between two completely independent datasets. The actual agreement across independent cohorts could be lower.

Overall, the results suggest that individual LIONESS edges can be quite noisy when only a small number of patients are available. The broader network-level patterns appear more reliable than individual edge selections.

These results are presented in `stability_subsampling.csv` and `fig4`.

#### C. Effect of Edge-Selection Threshold

Finally, we examined how different edge-selection thresholds affect the stability of the results.

We compared the top 0.1%, 1%, and 5% of edges.

The top 5% of edges showed better stability, with Jaccard scores ranging from 0.16 to 0.49, compared with 0.02 to 0.22 for the top 0.1%.

This was an interesting observation because we initially expected stricter thresholds to produce more reproducible edge lists.

However, the results showed the opposite. Selecting a very small number of top-ranked edges did not necessarily improve reproducibility.

In this dataset, retaining a larger proportion of edges resulted in more stable selections.

### 5.5 LIONESS vs Classical Differential Co-expression

To understand whether LIONESS offers an advantage over the traditional approach, we compared both methods using the same set of 500 genes.

The two methods showed a moderate correlation in their edge rankings, with a Spearman correlation of 0.69 between the absolute statistics.

However, their actual top-ranked edge lists showed relatively limited overlap:

| Edge-selection threshold | Jaccard overlap |
|---|---:|
| Top 0.1% | 0.004 |
| Top 1% | 0.10 |
| Top 5% | 0.24 |

This indicates that although the two methods may rank edges similarly to some extent, they do not necessarily identify the same top edges.

We also compared their stability across different patient subsets.

| Number of patients | LIONESS | Classical method |
|---|---:|---:|
| 9 | 0.07 | 0.13 |
| 18 | 0.17 | 0.27 |
| 27 | 0.35 | 0.44 |

The classical differential co-expression method showed higher stability across all three cohort sizes.

Another difference was the number of statistically significant edges identified by the two methods.

At 1,000 genes, the classical method identified 133,745 significant edges out of 499,500, with an empirical FDR of 0.066.

This is considerably higher than the number detected by the paired LIONESS analysis.

An important reason is that the classical method does not account for the fact that tumour and normal samples come from the same patient. In contrast, the LIONESS analysis uses paired statistical testing, taking the patient-level matching into account.

Therefore, the higher number of significant edges detected by the classical method should not automatically be interpreted as better performance.

Overall, the comparison showed that LIONESS provides sample-specific network information, but it was less stable than the classical method for identifying differential edges in this relatively small cohort.

## 6. Overall Conclusion

Through this project, we implemented LIONESS to generate sample-specific gene co-expression networks and applied it to a real pancreatic cancer gene expression dataset.

The analysis showed that there are detectable differences in gene co-expression patterns between pancreatic tumour and matched normal samples. LIONESS was able to capture these differences and provided an individual network for each sample, which can be useful when studying patient-level variation.

However, our additional experiments also highlighted some limitations. The individual network edges were relatively noisy when the cohort size was small, and their stability improved as more patients were included. We also found that the classical differential co-expression method was more stable under subsampling and identified a different set of top-ranked edges.

For the classification task, the LIONESS network features achieved an accuracy of 0.764, which was lower than the 0.889 achieved using gene expression values directly.

Therefore, for the specific question of identifying edges that differ between two groups, LIONESS did not outperform the classical method in this dataset.

This does not mean that LIONESS is ineffective. Its main advantage is that it provides sample-specific networks, making it useful for exploring individual differences and network heterogeneity. However, when the primary objective is simply to identify differential co-expression between two groups, our results suggest that the traditional method may be a more stable option under the conditions tested.

The implementation was also tested against the original LIONESS mathematical definition using the unit tests in `tests/test_lioness.py`.

An important outcome of this project is that the results highlight both the usefulness and the limitations of sample-specific network analysis rather than assuming that a more advanced method will always perform better.

## 7. Limitations and Potential Sources of Bias

There are several limitations that need to be considered when interpreting the results.

### 7.1 Stability Is Not the Same as Accuracy

We do not have a gold-standard set of true gene-gene interactions or network edges for this dataset. Therefore, stability, permutation testing, and classification accuracy are indirect measures of performance. We cannot determine which method identifies the biologically true network rewiring based on these metrics alone.

### 7.2 Gene Selection Bias

The genes were selected based on their variance across all samples. Since 811 out of the 1,000 selected genes were differentially expressed, the analysis is biased towards genes showing tumour-versus-normal expression differences.

### 7.3 Limited Dataset and Network Method

The experiments were performed using only one dataset and one network inference method, Pearson correlation. As a result, the findings may not generalise to other datasets or network construction approaches.

Also, Pearson-based LIONESS weights are not restricted to the range of -1 to 1. In our analysis, the values ranged approximately from -30 to 37, making them sensitive to outlying samples.

### 7.4 Statistical Dependence Between Edges

Edges that share genes are not independent of one another. Additionally, LIONESS networks are generated using overlapping sample data. These dependencies can affect the validity of the statistical tests and help explain the wide variation observed in the permutation null distribution.

### 7.5 Differences in Tissue Composition

The dataset contains whole-tissue microarray measurements, meaning that the expression values represent a mixture of different cell types.

Consequently, some of the observed network differences may be related to changes in tissue composition rather than changes in gene regulation. For example, tumour tissue may contain fewer acinar or islet cells compared with normal pancreatic tissue.

It is also important to remember that co-expression does not necessarily indicate a direct regulatory interaction between two genes.

### 7.6 Overlap in Subsampling Reference Sets

The patient subsets used in the stability analysis overlap with the full-cohort reference dataset. In particular, the 27-patient subsets share 75% of their samples with the reference set.

This means that the measured stability may be higher than what would be observed when comparing completely independent patient cohorts.

## References

1. Kuijjer ML, Tung MG, Yuan G, Quackenbush J, Glass K. Estimating sample-specific regulatory networks. *iScience*. 2019;14:226-240.

2. Badea L, Herlea V, Dima SO, Dumitrascu T, Popescu I. Combined gene expression analysis of whole-tissue and microdissected pancreatic ductal adenocarcinoma. *Hepatogastroenterology*. 2008;55:2016-2027. GEO accession: GSE15471.
