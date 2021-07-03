# raw 20500: candidate_sample.Rda -> summary table of HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda
# # "HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS.Rda" saved at "/home/tex/R/HNSCC_Tex_survival/hnscc_github/marginS/"
# n=6624 => n=6429; with Bonferroni and FDR
# # as HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR
# _marginS_ on [2019/07/02]
load(file="marginS/HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR.Rda")

# embase cancer genes removed: HNSCC_OS_marginS_THREE_pvalue005_noCancerGene.Rda
# pubmed_hnscc_genes371.Rda
# GeneSymbol on listed in TCPA_proteins 237 antibodies
#1        CD31
#2        CHK1
#3        DVL3
#4       EEF2K
#5        EGFR
#6        HER2
#7        HER3
#8        IRS1
#9      NOTCH1
#10       PTEN
#11        YAP
#12       BRAF
#13      FOXM1
#14      ERCC1
#15      KEAP1
#16       NRF2
# cancerGenes <- semi_join(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, TCPA_proteins, by = c("gene_id" = "GeneSymbol"))
# n=27, it is on the list of TCPA; could trying TCPA survival by protein array.
# volcano plot: Y axis -log10(p_value), X axis multi_HR
#HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$p_value_adj_bonferroni

#Creates bar graph that shows the high and low sample proportions from the volcano plot
library(ggplot2)
library(ggrepel)
#res3 <- 
#attach(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
hazards <- subset(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, select=c(gene_id, p_value, p_value_adj_bonferroni, multi_HR))
# names -> gene_id; HR -> multi_HR; log10P -> p_value
hazards$log10P <- -log10(as.numeric(hazards$p_value))
cut_pvalue <- 2.871800e-04 # -log10(0.05) = 1.3 or -log10(1.604667e-04) = 3.794615
ggplot(hazards, aes(x = multi_HR, y = log10P)) +     #Creation of ggplot graph with pre-determined aesthetics
  geom_point(shape = 21, alpha = 0.8, aes(size = log10P, fill = multi_HR)) + 
  scale_fill_distiller(palette = "RdYlGn", trans="log", limits = c(min(hazards$multi_HR), max(hazards$multi_HR))) +
#Diverging  BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
  geom_hline(yintercept = -log10(cut_pvalue), lty = 2) + 
  geom_vline(xintercept = 1) +
  guides(size=FALSE, fill=FALSE) + 
  guides(size=FALSE, fill=FALSE) + 
  xlab("Cox Regression Hazard Ratio") + 
  ylab("-log10(P-value)") + 
  ##  geom_label_repel(data=subset(hazards, log10P >= 1.3 & as.numeric(multi_HR) >= 1.5 | log10P >= 1.3 & as.numeric(multi_HR) <= 0.5),
#  geom_label_repel(data=subset(hazards, log10P >= -log10(cut_pvalue)),
  geom_label_repel(data=subset(hazards, gene_id %in% c("DKK1", "CAMK2N1", "STC2", "PGK1", "SURF4", "USP10", "NDFIP1", "FOXA2", "STIP1", "DKC1", 
                                                       "ZNF557", "ZNF266", "IL19", "MYO1H", "FCGBP", "LOC148709", "EVPLL", "PNMA5", "KIAA1683", "NPB")),
                   aes(label=gene_id), size=5, box.padding = unit(0.4, "lines"), segment.alpha=0.8, segment.size = 0.5)  + 
  theme_bw() + 
  theme(text=element_text(family="sans"),
        axis.title=element_text(size=18))
# legands on Our findings suggested 20 candidate biomarkers, 
# DKK1, CAMK2N1, STC2, PGK1, SURF4, USP10, NDFIP1, FOXA2, STIP1, DKC1, 
# as well as ZNF557, ZNF266, IL19, MYO1H, FCGBP, LOC148709, EVPLL, PNMA5, KIAA1683, and NPB, 
#detach(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# save as TCGA_HNSC_Optimal_Overall_allPlot_unKM_P_multiHR.pdf 取代 Figure 3
