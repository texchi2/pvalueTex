# pvalueTex
Supplementary Material of manuscript, title: "A Global Genome-wide Scan with Optimal Cutoff Mining for Emerging Biomarkers in Head and Neck Squamous Cell Carcinoma"

## A) R script: main and Cutoff Finding engine
```R
pvalueTex_R_lite.tar.gz, including
 Sep  8 2020 15:03 main_marginSFP_HNSCC.R		# main code
 Sep  8 2020 10:45 TCGA_HNSCC_marginSFP.R		# step 1
 May 20 2020 10:55 cutofFinder_func_HNSCC.R		# cutoff finding procedure
 Sep  8 2020 10:58 analysis_export.R.			# step 2
```
(current version at manuscript submitted on 2020/09/19)
(Sep 19 2020 02:59 volcano_plot.R)

## B) Rda repository
```R
 Sep  8  2020 15:01 HNSCC.clinical.Fire.Rda		# HNSCC clinical dataset from TCGA
 Aug 14  2019 HNSCC.clinical.RNAseq_tobacco.Fire.Rda	# HNSCC RNA-Seq combining clinical dataset (with tobacco exposure)
 Aug 14  2019 whole_genome.Rda				# "gene ID" of 20500 protein coding genes
```

## C) Supplementary raw tables:
```R
intermediate_data:
 Sep  8  2019 HNSCC.survival.marginS.20500.tar.gz		# resulting survival tables of each gene (.Rda + .xlsx), size 933Mb 
 Sep  8  2019 file20500_list.txt     				# gene ID and filename list of HNSCC.survival.marginS.20500.tar.gz
 Mar  8  2020 HNSCC_OS_marginS_6429.csv   			# survival tables of 6429 genes (uncorrected P-value < 0.05), with FDR correction
 Sep 19 2020 08:18 HNSCC.survival.marginS.candidate20.tar.gz  # survival tables of 20 candidate genes (.xlsx)
 Mar  8  2020 HNSCC_OS_marginS_10bad_candidates.csv		# gene ID list of "bad" candidates
 Mar  8  2020 HNSCC_OS_marginS_10good_candidates.csv		# gene ID list of "good" candidates
```

## Next on version 1.0.0
The implementation by Rstudio shiny App, supporting all TCGA diseases, with custom features

We use Semantic Versioning to identify each release. (https://semver.org)

