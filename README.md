# pvalueTex
[![DOI](https://zenodo.org/badge/295142503.svg)](https://zenodo.org/badge/latestdoi/295142503)

Supplementary Material of manuscript, title: "A Global Genome-wide Scan with Optimal Cutoff Mining for Emerging Biomarkers in Head and Neck Squamous Cell Carcinoma"

## A) R script: main and Cutoff Finding engine
```R
pvalueTex_R_lite.tar.gz, including
 Sep  8 2020 15:03 main_marginSFP_HNSCC.R		# main code
 Sep  8 2020 10:45 TCGA_HNSCC_marginSFP.R		# step 1
 May 20 2020 10:55 cutofFinder_func_HNSCC.R		# cutoff finding procedure (FDR)
 Jul  3 2021 10:58 analysis_export.R.			# step 2 (Bonferroni)
 Sep 19 2020 02:59 volcano_plot.R
 Jul  2 2021 22:00 Validation_GSE65858_survival.R  # final validation
 Jun 26 2021 20:45 Bayes_survival.R  # LASSO regression with Cox (RLassoCox)
 ..
2021/07/22 all R script R_script_2020July.zip

```
(first version at manuscript submitted on 2020/09/19)
(final version at manuscript submitted on 2021/07/04) JPM

## B) Rda repository
```R
 Sep  8  2020 15:01 HNSCC.clinical.Fire.Rda		# HNSCC clinical dataset from TCGA
 Aug 14  2019 HNSCC.clinical.RNAseq.Fire.Rda (size: 85031069byte; HNSCC RNA-Seq 20500 genes
 Aug 14  2019 HNSCC.clinical.RNAseq_tobacco.Fire.Rda	# HNSCC clinical dataset (with tobacco exposure)
 Aug 14  2019 whole_genome.Rda				# "gene ID" of 20500 protein coding genes
```

## C) Article and Keynote of the Video Abstract:
```R
Aug 11 2021 jpm-11-00782-v2.pdf # publised article, https://www.mdpi.com/2075-4426/11/8/782
Aug 15 2021 KOMA_Script_MDPI_jpm_1309345_Keynote_VideoAbstract.pdf # keynote for Video abstract
Video is also available at YouTube: https://youtu.be/VK5F43fusDA

```


## D) Supplementary Figures and unpulished raw tables:
```R

Intermediate_data:
 Sep  8  2019 HNSCC.survival.marginS.20500.tar.gz		# resulting survival tables of each gene (.Rda + .xlsx), size 933Mb 
 Sep  8  2019 file20500_list.txt     				# gene ID and filename list of HNSCC.survival.marginS.20500.tar.gz
 Mar  8  2020 HNSCC_OS_marginS_6429.csv   			# survival tables of 6429 genes (uncorrected P-value < 0.05), with FDR correction
 Sep 19 2020 08:18 HNSCC.survival.marginS.candidate20.tar.gz  # survival tables of 20 candidate genes (.xlsx)
 Mar  8  2020 HNSCC_OS_marginS_10bad_candidates.csv		# gene ID list of "bad" candidates
 Mar  8  2020 HNSCC_OS_marginS_10good_candidates.csv		# gene ID list of "good" candidates
 Oct  7  2020 CI_genes.Rda            candidate gene set belongs to the immune system process and immune response
 Jun  6  2021 HNSCC.mRNA.Exp.SLC35E2.Fire.Rda     mRNA data of SLC35E2 (containing SLC35E2A and SLC35E2B)
 
Supplementary Figures:
 Feb  2  2021 cancers-1082127_supplementaryFigure123.pdf
 Aug  4  2021 MDPI_JPM_1309345-10_supplementary.pdf
```

## Next on version 1.0.0
The implementation by Rstudio shiny App, supporting all TCGA diseases, with custom feature.

We use Semantic Versioning to identify each release. (https://semver.org)

## GPL-3.0-or-later
The content of this material is licensed under the GNU General Public License v3.0 or later.
