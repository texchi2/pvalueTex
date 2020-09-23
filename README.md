# pvalueTex
Supplementary Material of manuscript of "A Global Genome-wide Scan with Optimal Cutoff Mining for Emerging Biomarkers in Head and Neck Squamous Cell Carcinoma"

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
 Sep  8  2019 HNSCC.survival.marginS.20500.tar.gz	# resulting survival tables, size 933Mb 
 Sep  8  2019 file20500_list.txt				# gene ID and filename list of above
 HNSCC_OS_marginS_candidates_Bonferroni.xlsx		# resulting survival tables, filtered Bonferroni correction
 Mar  8  2020 candidates_20genes.txt			# gene ID list of candidates
 Sep 19 2020 08:18 HNSCC.survival.marginS.candidate20.tar.gz	# survival tables of candidate genes (.xlsx)
```

## Next version
The implementation by Rstudio shiny App, supporting all TCGA diseases, with custom features

