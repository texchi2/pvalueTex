# pvalueTex
Supplementary Material of manuscript:

## A) R script: main and Cutoff Finding engine (current version at manuscript submitted on 2020/09/19)
```R
pvalueTex_R_lite.tar.gz, including
 Sep  8 2020 15:03 main_marginSFP_HNSCC.R	# main code
 Sep  8 2020 10:45 TCGA_HNSCC_marginSFP.R	# step 1
 May 20 2020 10:55 cutofFinder_func_HNSCC.R	# cutoff finding procedure
 Sep  8 2020 10:58 analysis_export.R.		# step 2
```
(Sep 19 2020 02:59 volcano_plot.R)

## B) Rda repository
```R
-rw-r--r--  1 tex tex   669734 Sep  8  2020 15:01 HNSCC.clinical.Fire.Rda	# HNSCC clinical dataset from TCGA
-rw-r--r--  1 tex tex 85031069 Aug 14  2019 HNSCC.clinical.RNAseq.Fire.Rda	# HNSCC RNA-Seq combining clinical dataset
-rw-r--r--  1 tex tex   281570 Aug 14  2019 whole_genome.Rda				# "gene ID" of 20500 protein coding genes
```

## C) Supplementary raw tables:
```R
intermediate_data:
 Sep  8  2019 HNSCC.survival.marginS.20500.tar.gz 	# resulting survival tables, size 933Mb 
 Sep  8  2019 file20500_list.txt)					# gene ID and filename list of HNSCC.survival.marginS.20500.tar.gz
 HNSCC_OS_marginS_candidates_Bonferroni.xlsx		# resulting survival tables, filtered Bonferroni correction
 
 Mar  8  2020 candidates_20genes.txt				# gene ID list of candidates
 Sep 19 2020 08:18 HNSCC.survival.marginS.candidate20.tar.gz	# survival tables of candidate genes (.xlsx)
```

## Next version: implementation by Rstudio shiny App, supporting all TCGA diseases, with custom features

