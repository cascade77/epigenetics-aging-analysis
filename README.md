# Epigenetics and Aging Analysis

This repository contains two fully documented analyses covering different aspects of DNA methylation biology. The first is a whole genome bisulfite sequencing (WGBS) pipeline applied to breast cancer methylomes, following the Galaxy Training Network tutorial and based on the dataset from Lin et al. (2015). The second is an epigenetic aging clock benchmarking analysis using the Biolearn library, evaluating 8 aging clocks across two public EPIC array blood methylation datasets through correlation matrices, age deviation heatmaps, age prediction plots, MAE benchmarking, and predicted age distribution comparisons.

Both analyses are fully documented with biological rationale, tool explanations, and result interpretations. See the individual directory READMEs for detailed workflow descriptions, full code walkthroughs, and result interpretations for each analysis.

---

## Repository Structure

```
epigenetics-aging-analysis/
│
├── wgbs/                                        # WGBS bisulfite sequencing pipeline
│   ├── README.md                                # Detailed workflow and result interpretations
│   └── results/                                 # All output figures from Galaxy
│       ├── falco_qc_subset1.png                 # Per base sequence content, subset_1.fastq
│       ├── falco_qc_subset2.png                 # Per base sequence content, subset_2.fastq
│       ├── methylation_bias_top_strand.png      # MethylDackel mbias original top strand
│       ├── methylation_profile_subset.png       # plotProfile output from subset sample
│       ├── methylation_level_all_samples.png    # plotProfile output across all 6 samples
│       ├── metilene_dmr_mean_difference.png     # Distribution of DMR methylation differences
│       ├── metilene_dmr_length_nt.png           # DMR length distribution in nucleotides
│       ├── metilene_dmr_length_cpg.png          # DMR length distribution in CpG count
│       ├── metilene_dmr_qvalue.png              # Mean methylation difference vs q-value
│       ├── metilene_dmr_group_comparison.png    # Mean methylation group 1 vs group 2
│       └── metilene_dmr_nt_vs_cpg.png          # DMR length in nt vs length in CpGs
│
├── epic_array/                                  # EPIC array aging clock benchmarking
│   ├── biolearn_aging_clocks.ipynb              # Fully documented analysis notebook
│   ├── README.md                                # Detailed workflow and result interpretations
│   └── results/                                 # All output figures from Colab
│       ├── correlation_matrix_GSE120307.png     # Clock correlation matrix, dataset 1
│       ├── correlation_matrix_GSE41169.png      # Clock correlation matrix, dataset 2
│       ├── age_deviation_heatmap_GSE120307.png  # Age deviation heatmap, dataset 1
│       ├── age_deviation_heatmap_GSE41169.png   # Age deviation heatmap, dataset 2
│       ├── age_prediction_GSE120307.png         # Predicted vs chronological age, dataset 1
│       ├── age_prediction_GSE41169.png          # Predicted vs chronological age, dataset 2
│       ├── mae_comparison.png                   # MAE bar chart across both datasets
│       └── predicted_age_distribution.png       # Box plots of predicted age distributions
│
├── LICENSE
└── README.md
```

---

## Analyses

| # | Directory | Platform | Focus |
|---|-----------|----------|-------|
| 1 | [WGBS](wgbs/) | Galaxy Europe (usegalaxy.eu) | Bisulfite QC, alignment, methylation extraction, DMR detection |
| 2 | [EPIC Array](epic_array/) | Google Colab | Aging clock benchmarking, correlation, heatmaps, MAE comparison |

---

## Results

### WGBS: Bisulfite Sequencing Pipeline

Whole genome bisulfite sequencing pipeline applied to breast cancer and normal breast tissue samples from Lin et al. (2015), processed through the Galaxy Training Network methylation-seq tutorial.

![Falco QC subset1](wgbs/results/falco_qc_subset1.png)

Per base sequence content from Falco QC on subset_1.fastq. T is elevated to approximately 50% and C drops to nearly zero, the characteristic bisulfite conversion signature confirming successful treatment.

![Methylation bias top strand](wgbs/results/methylation_bias_top_strand.png)

MethylDackel methylation bias plot for the original top strand across both paired reads. CpG methylation remains consistent at 70 to 75% across all read positions, confirming no significant positional bias and no trimming was required.

![Methylation profile subset](wgbs/results/methylation_profile_subset.png)

plotProfile output showing average methylation levels around CpG islands for the subset sample. The dip at the TSS confirms that CpG island promoters remain hypomethylated in actively transcribed genes.

![Methylation level all samples](wgbs/results/methylation_level_all_samples.png)

plotProfile across all 6 samples (NB1, NB2, BT089, BT126, BT198, MCF7). Cancer samples show higher TSS methylation compared to normal breast tissue, consistent with the aberrant promoter hypermethylation reported in Lin et al.

![Metilene DMR mean difference](wgbs/results/metilene_dmr_mean_difference.png.png)

Distribution of mean methylation differences between normal breast (NB1, NB2) and invasive ductal carcinoma (BT198). The left-skewed distribution reflects predominant hypomethylation in cancer, with a smaller right peak of hypermethylated DMRs.

![Metilene DMR qvalue](wgbs/results/metilene_dmr_qvalue.png)

Mean methylation difference vs q-value for all detected DMRs. The most significantly differentially methylated regions reach q-values below 1e-100.

---

### EPIC Array: Aging Clock Benchmarking

Eight epigenetic aging clocks benchmarked across two blood EPIC array datasets using the Biolearn library.

![Correlation matrix GSE120307](epic_array/results/correlation_matrix_GSE120307.jpeg)

Clock correlation matrix for GSE120307 (34 samples). First generation chronological clocks (Horvathv1, Hannum, Lin) correlate above 0.90. DunedinPACE shows near-zero or negative correlations with all clocks as expected given it measures pace of aging rather than absolute age.

![Correlation matrix GSE41169](epic_array/results/correlation_matrix_GSE41169.jpeg)

Clock correlation matrix for GSE41169 (95 samples). The same clustering structure holds but correlations are generally lower, consistent with greater biological heterogeneity in the larger cohort.

![Age deviation heatmap GSE120307](epic_array/results/age_deviation_heatmap_GSE120307.jpeg)

Age deviation heatmap for GSE120307. DunedinPACE and Zhang_10 show the most extreme deviations. Most other clocks stay within approximately ±10 years of chronological age.

![Age deviation heatmap GSE41169](epic_array/results/age_deviation_heatmap_GSE41169.jpeg)

Age deviation heatmap for GSE41169. With 95 samples, consistent multi-clock age acceleration patterns are visible in individual samples, suggesting genuine biological age differences rather than clock noise.

![MAE comparison](epic_array/results/mae_comparison.jpeg)

Mean Absolute Error per clock across both datasets. Horvathv1 achieves the lowest MAE (~4 years). Zhang_10 has the highest MAE by a wide margin (~38 years on GSE120307), reflecting the cost of reducing a clock to 10 CpGs.

![Predicted age distribution](epic_array/results/predicted_age_distribution.jpeg)

Box plots of predicted age distributions per clock. Most clocks are centered near the mean chronological age. Zhang_10 is a clear outlier with median predictions near zero on GSE120307.

---

## Platform

WGBS analysis was run on Galaxy Europe (usegalaxy.eu). EPIC array analysis was run on Google Colab.

---

## Dependencies
```
WGBS (Galaxy tools)
Falco 1.2.4
bwameth 0.2.7
MethylDackel 0.5.2
Wig/BedGraph-to-bigWig 1.9.1
computeMatrix 3.5.4
plotProfile 3.5.4
Metilene 0.2.6.1
EPIC Array (Python)
pip install biolearn
```

---

## Citations

Lin, I.-H., Chen, D.-T., Chang, Y.-F., Lee, Y.-L., Su, C.-H., et al. (2015). Hierarchical Clustering of Breast Cancer Methylomes Revealed Differentially Methylated and Expressed Breast Cancer Genes. *PLOS ONE* 10(2): e0118453. https://doi.org/10.1371/journal.pone.0118453

Ying, K., Paulson, S., Perez-Guevara, M., Emamifar, M., Martinez, M.C., Kwon, D., Poganik, J.R., Moqri, M., and Gladyshev, V.N. (2023). Biolearn, an open-source library for biomarkers of aging. *bioRxiv*. https://doi.org/10.1101/2023.12.02.569722

Wolff, J., Ryan, D., Moosmann, V. (2017). DNA Methylation data analysis. *Galaxy Training Materials*. https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html

---

## References

[Galaxy Training Network: DNA Methylation data analysis](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/methylation-seq/tutorial.html)

[Galaxy Training Network: Introduction to DNA Methylation](https://training.galaxyproject.org/training-material/topics/epigenetics/tutorials/introduction-dna-methylation/slides-plain.html)

[Biolearn documentation](https://bio-learn.github.io/)

[Biolearn GEO data sources](https://bio-learn.github.io/data.html)

