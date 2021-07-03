# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
#############################################################0###
# Tex Since 2021/06/01
#   Differential expression analysis with limma
BiocManager::install("GEOquery")
BiocManager::install("BeadArrayUseCases")

library(GEOquery)
library(limma)
library(umap)

library(BeadArrayUseCases)

# http://bioconductor.org/packages/release/data/experiment/vignettes/BeadArrayUseCases/inst/doc/BeadArrayUseCases.pdf

GSE65858 <- getGEO('GSE65858', destdir=".")
# or 
#head(Meta(GSE65858, GSEMatrix=TRUE))
#View(GSE65858) 
# total 287 participants = 275 + 12 (paired 一人採二次)
#275 samples= 3 were taken from lymph node metastases,
# and 2 are cell lines and 270 (included here) are tumor samples. 
# 287- (12+2) = 273 = ?? + 168(alive) + 88(dead) = 256
# Status	Public on Jun 24, 2015
# Title	GW_001 Head and neck squamous cell carcinoma Illumina HT12v4
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1607684
# features of GSM1607684 (participant)
# tumor_type: Primary, 
# os: 1100 days  overall survival
# os_event: FALSE (alive)
# pfs: 1100  days pregression-free survival
# pfs_event: FALSE (cancer-free)

#Total number of rows: 31330 (RNA-Seq probe #)
#
#R script at https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE65858 GSM1607684
#each participant has whole genome RNA-Seq expression value
#ID_REF	VALUE(RNA-Seq)	Detection Pval by log2-transformed and normalized using
#RSN 
#ILMN_1762337	6.4281	0.2051948 
#ILMN_2055271	6.5589	0.02467532


# > GEO data prepare ####
# load series and platform data from GEO
gset <- getGEO(filename=file.choose(), GSEMatrix =TRUE, AnnotGPL=TRUE)
# or
gset <- getGEO("GSE65858", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable ####
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples####
gsms <- paste0("0010110000001010110010000000101111000001101011X110",
               "1001100111000000010001001001001011101000X100010000",
               "00000000110100010101001000100100010110011100010010",
               "00001001010000011100011000011001010001100001000001",
               "01000011100011000110000110000011XXXXXXXX1X11100010",
               "00000000100000X000XX")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")####
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel] # n=256

# log2 transformation, if not yet ####
ex <- exprs(gset) # 31330 (expression) x 256 (participants)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
# R4> qx
# [1]  5.3782  6.4694  6.7642  7.8070 11.9450 14.2618
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix####
gs <- factor(sml)
groups <- make.names(c("censored","OS_event")) # (0 alive, 1 dead)
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset) # dummy: grouping GSMs by 0 1
colnames(design) <- levels(gs) # "censored","OS_event")

## > extract case clinical data from gset ####
##  (GSE is different from GDS)
## https://www.bioconductor.org/packages/release/bioc/vignettes/GEOquery/inst/doc/GEOquery.html
# patient_ID = GSM...

clinical <- phenoData(gset) # *** check list() to find phenoData()

# !tumor_type: Secondary
gset_OS_mRNA <- clinical@data[!(clinical@data$characteristics_ch1.6 == "tumor_type: Secondary"), 
                              c("geo_accession", "characteristics_ch1.1", "characteristics_ch1.2", 
                                  "characteristics_ch1.6", 
                                  "characteristics_ch1.15", "characteristics_ch1.16")]
# x cbind of gset_OS_mRNA and design (withn dummy OS time)
design <- as.data.frame(design)
design$"geo_accession" <- rownames(design)
gset_OS_mRNA <- merge(gset_OS_mRNA, design, by = "geo_accession")



# >> from pvalueTex, select fdr100 top genes  ####
# candidate_sample at analysis_export.R [line 535] 
# # or candidate_fdr100, candidate_fdr2000, or candidate_fdr6429
#x fdr100_FDRPvalue <- read.csv(file="FDR100.csv", header = T)
#colnames(candidate_fdr100)[2] <- "Gene.symbol"
candidate_fdr100 <- candidate_fdr6429 # all candidate_sample here
candidate_fdr100 <- candidate_fdr2000 # with FDR P value from KM cutoff finding
colnames(candidate_fdr100)[2] <- "Gene.symbol"


# ID_REF to gene symbol converting
#exp_gene_symbol <- featureData(gset) # with gene symbol
# (ID) ILMN_1651209  as  (Gene.symbol) SLC35E2 (no SLC35E2A/B)
convert_ILMN_ID_gene_symbol <- 
  featureData(gset)@data[, c("ID", "Gene.symbol")]

# more than 1 probe on each gene
# n=100 genes (2000)(6429) -> n=116 ID list (2292)(7614)
candidate_fdr100 <- merge(candidate_fdr100, subset(convert_ILMN_ID_gene_symbol, select=c("Gene.symbol", "ID")))
# it has ILMN and gene symbol

# attaching expression values (by ILMN ID)
exp31330 <- as.data.frame(exprs(gset)) # 31330 *  256
#exp31330 <- as.data.frame(t(exp31330)) #  256 * 31330(ILMN ID)

# geo_accession (GSM, samples), ILMN ID (gene)
exp31330$"ID" <- rownames(exp31330) # at 257 column; shown "V1"

# now cohort (n=257) has OS event + FDR100 gene expression
candidate_fdr100 <- merge(candidate_fdr100, exp31330, by="ID")
# t transpose
exp257 <- as.data.frame(t(candidate_fdr100[, c(2, 5:260)])) # 2: gene symbol
# colnames(exp257) <- exp257[1, ] => there is duplicated gene symbol (due to 1 to many probes)
#exp257 <- exp257[-1, ]
exp257$"geo_accession" <- rownames(exp257)

# rownames as geo_accession (samples ID)
# to match their nrow:
df <- as.data.frame(matrix(colnames(gset_OS_mRNA)[1:6], ncol=6))
colnames(df) <- colnames(gset_OS_mRNA)[1:6]
gset_OS_mRNA <- rbind(df, gset_OS_mRNA[, 1:6])

#$$$$$$$$
# **** make sure nrow is 257, if 258
gset_OS_mRNA <- gset_OS_mRNA[-1, ]
#$$$$$$$$
rownames(gset_OS_mRNA)[1] <- exp257$geo_accession[1] # it is also a "case"
gset_OS_mRNA$geo_accession <- as.character(gset_OS_mRNA$geo_accession)
gset_OS_mRNA$geo_accession[1] <- exp257$geo_accession[1]
gset_OS_mRNA2 <- merge(gset_OS_mRNA, exp257, by = "geo_accession")

gset_OS_mRNA3 <- gset_OS_mRNA2
gset_OS_mRNA3 <- remove.factors(gset_OS_mRNA3)
colnames(gset_OS_mRNA3)[7:ncol(gset_OS_mRNA3)] <- 
                      as.character(gset_OS_mRNA3[1, 7:ncol(gset_OS_mRNA3)])
# with gene symbol at colname 7:end
gset_OS_mRNA2 <- gset_OS_mRNA3[-1, ]
save(gset_OS_mRNA2, file="gset_OS_mRNA2.Rda")
# going KM survival




#############################################################3###
#  Kaplan-Meier survival analysis####
#    
## generate life table with each mRNA expression ###
## gset_OS_mRNA2: HNSCC cohort with ID(GSM), OS event, OS time, RNA-Seq value
## dim 256 x (4+30330) -> 122

#OS1 from GSE65858
# gset_OS_mRNA2 merged OS time and gene expression (FDR100)
# 
# 

# from pvalueTex
# > fit KM survival curves ####
# cancer type should be defined at TCGA_cohort <- "LUAD" "HNSC"
# OS time (days), OS event
colnames(gset_OS_mRNA2)[c(5,6)] <- c("OS_time", "OS_event")

library(stringr) # extract last word
df <- as.data.frame(matrix(NaN, ncol = 1, nrow = nrow(gset_OS_mRNA2))) # Don't skip 1st row 
df$OS_time <- str_extract(gset_OS_mRNA2$OS_time,"\\w+$")
df$OS_event <- str_extract(gset_OS_mRNA2$OS_event,"\\w+$") 
# if # extract first word: "\\w+" in regular expression

df$OS_time <- as.numeric(df$OS_time) # days
df$OS_event <- ifelse(df$OS_event == "TRUE", 1, 0) # OS event TRUE as 1 (dead)
# 1: 88 vs 0:168
gset_OS_mRNA2[, c("OS_time", "OS_event")] <- df[, c("OS_time", "OS_event")]

coln7 <- 7 # exp of gene ZNF589
coln122 <- ncol(gset_OS_mRNA2) # 2298 now for FDR2000
# probes number: coln122-coln7+1 =116 (100 genes)
library('taRifx')
gset_OS_mRNA3 <- remove.factors(gset_OS_mRNA2)
# unlist or lapply if legnth > 1 to coerce
gset_OS_mRNA3[, c(coln7:coln122)] <- lapply(gset_OS_mRNA3[, c(coln7:coln122)], as.numeric)
#gset_OS_mRNA3 <- as.numeric(gset_OS_mRNA3[, c(coln7:coln122)])
gset_OS_mRNA2 <- gset_OS_mRNA3
save(gset_OS_mRNA2, file="gset_OS_mRNA2.Rda")



# run from here => go ####
# *************
# removal of all generated pdf
library(survival)
# 
mysurv <- Surv(gset_OS_mRNA2$OS_time, gset_OS_mRNA2$OS_event == 1) #1==dead, time by days

# check MASP1, warnings() on more than one probe to a gene
#which("MASP1"== colnames(gset_OS_mRNA2))
# [1] 1503 1580 1928
#
# for loop (100 genes)
fdr100_KMsurvival <- data.frame("Gene Symbol"=NaN, "FDR Pvalue"=NaN, 
                                   "KM Pvalue"=NaN, "cases_High"=NaN, "cases_Low"=NaN)
# a empty 1 x 5 data.frame
#i <- 7
#
for (i in c(coln7:coln122)) {
# expression cutoff at median value -> grouping by dichomolized it
# cutoff <- sapply(gset_OS_mRNA2[, i], median)
cutoff <- median(gset_OS_mRNA2[, i])
# dichotomize
#oscc$ageDx[gset_OS_mRNA2[, i] < cutoff] <- 0 # underexpression
#oscc$ageDx[gset_OS_mRNA2[, i] >= cutoff] <- 1 # overexpression than cutoff
gset_OS_overexpression <- ifelse((gset_OS_mRNA2[, i] >= cutoff), 1, 0)

# Test for difference (log-rank test) between groups (by PMM1_median 0 vs 1)
#tryCatch(surv_OS1 <- survdiff(mysurv ~ as.vector(gset_OS_mRNA2[, gset_OS_mRNA2M_pos], mode="numeric"), data=gset_OS_mRNA2), error = function(e) return(NA)) # PMM1 high or low
# or OS.km ; as.vector(unlist(gset_OS_mRNA2[, osccCleanNAM_pos]), mode="numeric"), data=gset_OS_mRNA2
# survfit(, type= "kaplan-meier", conf.type = "log-log") 
# for fit, plot; 
surv_OS1 <- survfit(mysurv ~ gset_OS_overexpression, type= "kaplan-meier", conf.type = "log-log") # log-rank test

# survdiff for two group with P value
# pchisq gives the distribution function
surv_OSp <- survdiff(mysurv ~ gset_OS_overexpression) # log-rank test
# N Observed Expected (O-E)^2/E (O-E)^2/V
# broom::glance(surv_OS1)$p.value
p_OS1 <- format(pchisq(surv_OSp$chisq, length(surv_OSp$n)-1, lower.tail = FALSE), digits=3)
# (p_OS1 == p_OS0) is TRUE
#cases_OS1 <- surv_OSp$n[1]
# #Value of survdiff => a list with components: help("survdiff")
# n => the number of subjects in each group.
# obs
# the weighted observed number of events in each group. If there are strata, this will be a matrix with one column per stratum.
# exp
# the weighted expected number of events in each group. If there are strata, this will be a matrix with one column per stratum.
# chisq
# the chisquare statistic for a test of equality.
# var
# the variance matrix of the test.
# strata
# optionally, the number of subjects contained in each stratum.

# [survfit] - Kaplan-Meier curve: P-Value ####
# confidence intervals as log hazard or log(-log(survival))
#OS.km <- survfit(mysurv ~ as.vector(unlist(gset_OS_mRNA2[, gset_OS_mRNA2M_pos]), mode="numeric"), data=gset_OS_mRNA2, type= "kaplan-meier", conf.type = "log-log")
# 365.25 days in a year for xscale => 3650 days for 10 years
# maximal followup to 6.8 years => 2483.7 days
# 12 months per year for 5 years => 60 months, 1826.25 days
# xscale=12, xmax=60

pdf_fn <- paste("rplot_GSE65858_", colnames(gset_OS_mRNA2)[i], ".pdf", sep="")
if (file.exists(pdf_fn)) {
  pdf_fn <- 
  paste("rplot_GSE65858_", colnames(gset_OS_mRNA2)[i], "_", i, ".pdf", sep="")
  }
pdf(pdf_fn) # gene symbol
# save PDF, size 746 x 431 pixel(?)
#
plot(surv_OS1, lty=1, xscale=365.25, xmax=1900, col=c("blue","red"), 
     sub=paste("Kaplan-Meier P Value =", p_OS1), 
     main=paste("OS in GSE65858", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", colnames(gset_OS_mRNA2)[i]), ylab="Percent Survival", xlab="Years")
legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:1, col=c("blue","red"))
dev.off()
#plot(OS.km, lty=1, col=c("blue","red"), sub="p= 0.816", main="OS in gset_OS_mRNA2(n=505)/gene level", ylab="Percent Survival", xlab="Days")
#legend("topright", legend=c('Low', 'High'), lty=1:2, col=c("blue","red"))
# summary(OS.km, times = seq(0, 3000, 100))
# 
# save P value
# [1] "Gene.Symbol" "FDR.Pvalue"  "KM.Pvalue"   "cases_High"  "cases_Low" 
fdrpvalue <- candidate_fdr100$p_value[which(candidate_fdr100$Gene.symbol==colnames(gset_OS_mRNA2)[i])]
if (length(fdrpvalue) > 1) {fdrpvalue <- fdrpvalue[1]} # prevent duplicated genes (FDR Pvalue)
fdr100_KMsurvival[i-which(colnames(gset_OS_mRNA2) == "OS_event"), ] <- 
      c(colnames(gset_OS_mRNA2)[i], 
        fdrpvalue, 
        p_OS1, surv_OS1$n[2], surv_OS1$n[1])
print(paste("[", i, "/", coln122, "] Generating median-cut KM plot of", colnames(gset_OS_mRNA2)[i]))
      
            #"(", which(whole_genome==geneName), ")"))

} # end of for loop

View(fdr100_KMsurvival)
# 12/100 or 152/2000 discovery rate by GSE65858
# 530/7614
plot(fdr100_KMsurvival$FDR.Pvalue[fdr100_KMsurvival$KM.Pvalue < 0.05], fdr100_KMsurvival$KM.Pvalue[fdr100_KMsurvival$KM.Pvalue < 0.05], log = "y")
# % scp  tex@35.201.169.0:~/R/HNSCC_Tex_survival/hnscc_github/rplot_GSE65858.tar.gz ./Downloads



# > KM plot of GSE65858 ####
library(ggplot2)
### not 最棒 P-value plot of KM survival analyses 
### (PvalueplotKM_20genes tiff or png 600*500 ** with Bonferroni correction 決定)
### 
### 這張圖最棒 Figure 6??
#要標出 10+10 gene list on this P-value plot
# candidates_bad_guy
# candidates_good_guy
# #number/OS_pvalue: from cutoff finder, the number (frequency) of OS P-values, which KM P-value < 0.05, in this gene;
# higher probability to be "found" significantly on a random selection of "cutoff" value in the tranditional manner.
#detach("package:ggplot2", unload=TRUE)
#x HNSCC_OS_marginS_pvalue_sorted$p_value_adj_bonferroni <- p.adjust(HNSCC_OS_marginS_pvalue_sorted$p_value, method="bonferroni")
# "HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS.Rda" saved at "/home/tex/R/HNSCC_Tex_survival/hnscc_github/marginS/"
# 這才是 Bonferroni?
## [1] "number"                 "gene_id"               
# [3] "p_value"                "p_value_adj_bonferroni"
#HNSCC_OS_marginS_pvalue_sorted_noNA_bonferroni <- HNSCC_OS_marginS_pvalue_sorted[complete.cases(HNSCC_OS_marginS_pvalue_sorted), ] # removal of NAs
#attach(HNSCC_OS_marginS_pvalue_sorted_noNA_bonferroni) # n=20500-13876(NAs)=6624, uncorrected P-value
# Bonferroni_cutoff <- alpha_HNSCC / (LUAD_n * n_percent_Bonferroni) = 0.05/9393 = 5.323113e-06
# plot() axes scale with any transformation, not just logarithmic (for Weibull plots)
# min x as 3.783-6
# 
# Code from [volcano_plot.R] 2020/09/19
#Creates bar graph that shows the high and low sample proportions from the volcano plot
#library(ggplot2)
library(ggrepel)

# raw 20500: candidate_sample.Rda -> summary table of HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda
# # "HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS.Rda" saved at "/home/tex/R/HNSCC_Tex_survival/hnscc_github/marginS/"
# n=6624 => n=6429; with Bonferroni and FDR
# # as HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR
# _marginS_ on [2019/07/02]
#load(file="marginS/HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR.Rda")


fdr100_KMsurvival$FDR.Pvalue <- as.numeric(fdr100_KMsurvival$FDR.Pvalue)
fdr100_KMsurvival$KM.Pvalue <- as.numeric(fdr100_KMsurvival$KM.Pvalue)
hazards$KM.Pvalue[hazards$Gene.Symbol=="IL19"]
#[1] 0.03390 0.00269
#R4> hazards$KM.Pvalue[hazards$Gene.Symbol=="CAMK2N1"]
#[1] 0.00687
# "FCGBP" => [1] 0.0102
#  hazards$Gene.Symbo[hazards$FDR.Pvalue < cut_pvalue]
#[1] "IL19"    "FAM3D"   "CALML5"  "CAMK2N1" "IL19"   
#[6] "MASP1"   "FCGBP"  

# [1] "Gene.Symbol" "FDR.Pvalue" of PvalueTex [3] "KM.Pvalue" of GSE65858
#attach(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
#hazards <- subset(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, select=c(gene_id, p_value, p_value_adj_bonferroni, multi_HR))
hazards <- subset(fdr100_KMsurvival, KM.Pvalue < 0.05,
                  select=c(Gene.Symbol, FDR.Pvalue, KM.Pvalue))

# names -> Gene.Symbol; (x) HR -> multi_HR -> KM.Pvalue; (y) log10P -> FDR.Pvalue
hazards$log10P <- -log10(as.numeric(hazards$FDR.Pvalue)) # -log10 transforming
cut_pvalue <- 2.871800e-04 # -log10(0.05) = 1.3 or -log10(1.604667e-04) = 3.794615
ggplot(hazards, aes(x = KM.Pvalue, y = log10P)) +     #Creation of ggplot graph with pre-determined aesthetics
  geom_point(shape = 21, alpha = 0.8, aes(size = log10P, fill = KM.Pvalue)) + 
  scale_fill_distiller(palette = "RdYlGn", trans="log", limits = c(min(hazards$KM.Pvalue), max(0.05))) +
  #Diverging  BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
  geom_hline(yintercept = -log10(cut_pvalue), lty = 2) + 
  geom_vline(xintercept = 0.05) +
  guides(size=FALSE, fill=FALSE) + 
  guides(size=FALSE, fill=FALSE) + 
  xlab("Kaplan-Meier estimated P value (GSE65858)") + 
  ylab("-log10(FDR-corrected P value)") + 
  ##  geom_label_repel(data=subset(hazards, log10P >= 1.3 & as.numeric(multi_HR) >= 1.5 | log10P >= 1.3 & as.numeric(multi_HR) <= 0.5),
  #  geom_label_repel(data=subset(hazards, log10P >= -log10(cut_pvalue)),
  geom_label_repel(data=subset(hazards, Gene.Symbol %in% c("DKK1", "CAMK2N1", "STC2", "PGK1", "SURF4", "USP10", "NDFIP1", "FOXA2", "STIP1", "DKC1", 
                                                       "ZNF557", "ZNF266", "IL19", "MYO1H", "FCGBP", "LOC148709", "EVPLL", "PNMA5", "KIAA1683", "NPB")),
                   aes(label=Gene.Symbol), size=5, box.padding = unit(0.4, "lines"), segment.alpha=0.8, segment.size = 0.5)  + 
  theme_bw() + 
  theme(text=element_text(family="sans"),
        axis.title=element_text(size=18))
# legands on Our findings suggested 20 candidate biomarkers, 
# DKK1, CAMK2N1, STC2, PGK1, SURF4, USP10, NDFIP1, FOXA2, STIP1, DKC1, 
# as well as ZNF557, ZNF266, IL19, MYO1H, FCGBP, LOC148709, EVPLL, PNMA5, KIAA1683, and NPB, 
#detach(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# save as TCGA_HNSC_Optimal_Overall_allPlot_unKM_P_multiHR.pdf 取代 Figure 3

#R4> hazards$Gene.Symbo[-log10(hazards$FDR.Pvalue) >2.5]
#[1] "ATP13A4"    "PLAU"       "IL19"      
#[4] "ZNF662"     "POMP"       "DOT1L"     
#[7] "EPHX3"      "IL34"       "FAM3D"     
#[10] "LYPD2"      "SERPINE1"   "CALML5"    
#[13] "RASIP1"     "FUT6"       "BCAR3"     
#[16] "ST6GALNAC1" "AQP1"       "CAMK2N1"   
#[19] "AIG1"       "IL19"       "FAM151B"   
#[22] "ABCB1"      "SMPX"       "MASP1"     
#[25] "FCGBP"      "GRIA3"     



########k########k#
#> Cox proportional hazard regression modeling for GSE65858 ####
# 不怕麻煩，勇敢畫圖
library(survival)
# hazard0 <- subset(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, select=c(gene_id, multi_HR)) # p_value, p_value_adj_bonferroni, 
#colnames(hazard0)[1] <- "Gene.Symbol"
#hazards <- merge(fdr100_KMsurvival, hazard0, by = "Gene.Symbol")
# 7614 probes of gene -> Cox modeling
save(fdr100_KMsurvival, file="fdr7614_KMsurvival_GSE65858.Rda")
fdr100_KMsurvival
# gset_OS_mRNA2.Rda; gset_OS_mRNA2 has survival OS event, OS time and expression values (7620-6)
load(file="gset_OS_mRNA2.Rda")
# osccCleanNA => gset_OS_mRNA2

# code from 7)[COXPH modelling]@TCGA_HNSCC_marginSFP.R
#*** [Univariate for OS]  table 3 left panel
# # 9 features in HNSCC, put in coxph one by one
# number of features 13 -> 14 -> 15 (add a feature: margin status, tobacco exposure)
oscox <- 0 # reset
## boundary of features column from first to last # 
#LL <- which(colnames(gset_OS_mRNA2) == "Gender") # "left" feature
#RR <- which(colnames(gset_OS_mRNA2) == "tobacco") # new "right" feature
#pmm1_pos <- which(colnames(gset_OS_mRNA2) == paste("PMM1", "_median", sep=""))
coln7 <- 7 # exp of gene GAPDH; LL
coln122 <- ncol(gset_OS_mRNA2) # 7620 now for FDR6429; RR
# probes number: coln122-coln7+1 =116 (100 genes) => 7614 probes
#gset_OS_mRNA2$OS_time
#gset_OS_mRNA2$OS_event # 1 as dead

# dichotomize gene expression (>= or < median)
for (i in c(coln7:coln122)) {
  cutoff <- median(gset_OS_mRNA2[, i])
  gset_OS_overexpression <- ifelse((gset_OS_mRNA2[, i] >= cutoff), 1, 0)
  gset_OS_mRNA2[, i] <- gset_OS_overexpression
}

features_os <- colnames(gset_OS_mRNA2)[c(coln7:coln122)] # gene names; x pmm1_pos)] # features selection ["Gender" to "tobacco", "PMM1"] 9 out of 15
osHR <- data.frame(matrix(0, nrow=length(features_os), ncol=4))
## *** looping the cox regression model over several features (genes)
## https://stackoverflow.com/questions/13092923/looping-cox-regression-model-over-several-predictor-variables
## as.formula: text to code (class: forumla list)
coxph_func <- function(x) as.formula(paste("Surv(gset_OS_mRNA2$OS_time, gset_OS_mRNA2$OS_event==1)", x, sep="~"))
formlist <- lapply(features_os, coxph_func)
#coxph_func2 <- function(x) as.formula(paste(x, ',', "data=osccCleanNA", sep=""))
#formlist2 <- lapply(formlist, coxph_func2)
#oscox <- coxph(coxph_func(x), data = osccCleanNA, ties="efron")
oscox <- lapply(formlist, coxph, data = gset_OS_mRNA2) # run coxph, "data=" as additional arguments to FUN. 
#-> model oscox (list 8)
# warnings() In FUN(X[[i]], ...) : X matrix deemed to be singular; variable 2
# correction on [2019/08/14]
os_coef_func <- function(i, oscox1) {
  x1 <- data.frame(summary(oscox1[[i]])$conf.int)[-2]
  x2 <- data.frame(summary(oscox1[[i]])$coefficients)[5]
  return(list(x1, x2))
}
#osHR <- rbind(osHR, cbind(unlist(uni_CI)))
uni_CI <- lapply(c(1:length(oscox)), os_coef_func, oscox) # for (i = 1:8)
for (ii in c(1:length(uni_CI))) {
  # ii=7 has 8 values in two rows <- marginNA issue => pick up first row by index [1, ]
  osHR[ii, ] <- t(c(unlist(uni_CI[[ii]][[1]][1,]), unlist(uni_CI[[ii]][[2]][1,])))
}
# > ii<-7; print(unlist(uni_CI[[ii]]))
# exp.coef.1   exp.coef.2   lower..951   lower..952 
# 3.4057062147           NA 1.7000673239           NA 
# upper..951   upper..952    Pr...z..1    Pr...z..2 
# 6.8225738228           NA 0.0005463051           NA
# => marginNA with NA

#skipNA <- rownames(osHR) %in% c(1) #, "marginNA")
#osHR <- osHR[which(!skipNA), ]
#ciUniMti <- c("Features",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value",	"HR",	"CI95%(L)",	"CI95%(H)",	"P-value")
colnames(osHR) <- c("uni_HR",	"CI95%(L)",	"CI95%(H)",	"Cox.P-value")

# correction of rownames
#featuresUni <- c("Gene expression")
# [1] "Gender"                 "Age at diagnosis"      
# [3] "Clinical T Status"      "Clinical N Status"     
# [5] "Clinical M Status"      "Clinical Stage"        
# [7] "Surgical Margin status" "Tobacco Exposure"
featuresUni <- colnames(gset_OS_mRNA2)[c(coln7:coln122)]
# duplicate 'row.names' are not allowed
#x rownames(osHR) <- featuresUni
osHR <- cbind(data.frame(matrix(featuresUni, ncol=1)), osHR)
colnames(osHR)[1] <- "Gene.Symbol"
# round
osHR <- round(osHR, 3)
# P-value notes:
# Significant codes:  0 ???***??? 0.001 ???**??? 0.01 ???*??? 0.05 ???.??? 0.1 ??? ??? 1
# if <0.001 => mark as "***"
osHR$`Cox.P-value`[osHR$`Cox.P-value` < 0.001] <- "***"

View(osHR) #-> (univariate) table 3 left panel

# where is this .Rda??
save(osHR, file="fdr7614_Cox_HR_GSE65858.Rda") # HR; as osHR
save(osHR, file="HNSCC_OS_GSE65858_pvalueCox_HR.Rda") # with Cox's uni_HR; still as osHR

# HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR from pvalueTex
# HNSCC_OS_GSE65858_pvalueCox_HR as osHR in 7614 probes of RNA-Seq






#> ggplot(Cox uni_HR of GSE65858) ####
# GSE65858 Cox, uni_HR
load(file="HNSCC_OS_GSE65858_pvalueCox_HR.Rda") # as osHR$uni_HR, ... etc.
hazard0 <- subset(osHR, select=c(Gene.Symbol, uni_HR)) # Cox.P-value 
load(file="fdr7614_KMsurvival_GSE65858.Rda") # as fdr100_KMsurvival, from from pvalueTex
# go
# FDR Pvalue (TCGA pvalueTex), KM Pvalue/Cox uni_HR (GSE65858)
hazards <- merge(fdr100_KMsurvival, hazard0, by = "Gene.Symbol")

# removal of duplicated in KM P value
#dup_KM <- !duplicated(hazards$KM.Pvalue)
#dup_GS <- !duplicated(hazards$Gene.Symbol)
# nrow = 7610 probes
hazards <- hazards[!duplicated(hazards), ]

# names -> Gene.Symbol; (x) HR -> multi_HR (-> KM.Pvalue); (y) log10P -> FDR.Pvalue -> KM.Pvalue
hazards$log10P <- -log10(as.numeric(hazards$KM.Pvalue)) # positive (for spot size)
hazards$plog10P <- log10(as.numeric(hazards$KM.Pvalue)) # becoming negative (for X-axis)

gene_labeled_golden <- c("CAMK2N1", "CALML5", "IL19", 
                         "FCGBP")

# from pvalueTex without gene_labeled_golden
gene_labeled <- c(
#  "DKK1", 
#"CAMK2N1", 
# "STC2", "PGK1", "SURF4", #"USP10", "NDFIP1",
#"FOXA2", "STIP1", "DKC1", 
#                  "ZNF557", "ZNF266", 
#"IL19", "FCGBP",
# "MYO1H", "LOC148709", "EVPLL", "PNMA5", "KIAA1683", 
#"NPB",

#    "ATP13A4",    
"PLAU",      # it is in Bayes Cox-Lasso
#"ZNF662",     
#"POMP",       
#"DOT1L",     
#"EPHX3" ,     #"IL34",       
"FAM3D",     
"MSMB", "LYPD2", "RBM11", # *** top KM P value in GSE65858 
#"SERPINE1",   #"CALML5",    
# "RASIP1",    # "FUT6",       
#"BCAR3",     
# "ST6GALNAC1", #"AQP1",   
 #"AIG1",        
# "FAM151B"
# "ABCB1",      
#"SMPX"
#"MASP1" ,   
# "GRIA3"
# CoxHR > 1.81:
"ADAMTSL2",
"BMP6" , "DUSP6" ,
"FXYD5" ,   "GPX8",
#"HERPUD2",
  "IFIT2" ,  
#"KANK4",   
"MYO1B" ,
 "NFKBID",   "PDIA5" ,   "PMEPA1",  "RAB43"  , 
# "RET" ,  
"SPP1" , 
"TCP11L1",  "TIMM10"  ,
"TNFAIP6" , # "TNPO1"  ,
#"TPD52L2", 
"TRAPPC10",
 "TTYH3" ,   "XYLT1" 
)

#R4> hazards$Gene.Symbol[hazards$uni_HR>=1.81 & hazards$KM.Pvalue<0.05]
# red spot above HR 1.81
#[1] "ADAMTSL2" "BMP6"     "CAMK2N1"  "DUSP6"  
#[6] "FXYD5"    "FXYD5"    "FXYD5"    "GPX8"     "GPX8"    
#[11] "HERPUD2"  "IFIT2"    "KANK4"    "MYO1B"    "MYO1B"   
#[16] "NFKBID"   "PDIA5"    "PMEPA1"   "PMEPA1"   "RAB43"   
#[21] "RET"      "SPP1"     "SPP1"     "TCP11L1"  "TIMM10"  
#[26] "TNFAIP6"  "TNPO1"    "TNPO1"    "TPD52L2"  "TRAPPC10"
#[31] "TTYH3"    "XYLT1" 

# pickup HR > 1.5 (in y axis)
hazards1 <- hazards # temporary save here
#hazards <- hazards[hazards$uni_HR > 1.5, ] # 266 probes
#hazards$golden <- hazards$Gene.Symbol %in% gene_labeled_golden # TRUE or FALSE, for scale_color_manual
#cols <- c("TRUE" = "red","FALSE" = "grey50") # named vector



# * FDR correction of KM P value of GSE65858 ####
# unadjusted P value < 0.05 (534 probes) -> fdr
#p.adjust(HNSCC_OS_marginS_pvalue_sorted$p_value, method="bonferroni")
hazards$FDR.KMpvalue[hazards$KM.Pvalue < 0.05] <- 
    p.adjust(hazards$KM.Pvalue[hazards$KM.Pvalue < 0.05], method = "fdr")
#R4> hazards[hazards$Gene.Symbol %in% gene_labeled_golden, c(1,3,10,6,2)]

save(hazards, file="HNSCC_GSE65858_FDRKMpvalue_7614probes.Rda")

cut_pvalue <- 0.05
  #2.871800e-04 # -log10(0.05) = 1.3 or -log10(1.604667e-04) = 3.794615
# x-y transposed
# for GSE65858 Cox HR (multi_HR -> uni_HR)
# KM unajdusted P value:
ggplot(hazards, aes(y = uni_HR, x = plog10P)) +   # x = plog10P  #Creation of ggplot graph with pre-determined aesthetics
  geom_point(shape = 21, alpha = 0.8, aes(size = log10P, fill = uni_HR)) + 
  scale_fill_distiller(palette = "RdYlGn", trans="log", limits = c(min(hazards$uni_HR), max(hazards$uni_HR))) + # min(hazards$uni_HR)
  #Diverging  BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
  geom_vline(xintercept = log10(cut_pvalue), color="red", lty=2) + 
  geom_hline(yintercept = 0.6, lty = 2) +
  geom_hline(yintercept = 1.8, lty = 2) +
  guides(size=FALSE, fill=FALSE) + 
  guides(size=FALSE, fill=FALSE) + 
  ylab("Cox Regression Hazard Ratio") + 
  xlab("log10(P value) of KM estimate (GSE65858)") + # -log10(P-value)
  ##  geom_label_repel(data=subset(hazards, log10P >= 1.3 & as.numeric(multi_HR) >= 1.5 | log10P >= 1.3 & as.numeric(multi_HR) <= 0.5),
  #  geom_label_repel(data=subset(hazards, log10P >= -log10(cut_pvalue)),

  scale_x_continuous(
    expand = expansion(mult = 0.3)) +
  
  # for gene_labeled_golden (CAMK2N1, IL19, and FCGBP)
  geom_label_repel(
    data=subset(hazards, Gene.Symbol %in% gene_labeled_golden),
#    data          = subset(dat, wt > 3),
    aes(label=Gene.Symbol),
    color = "red",
    size = 6,
    nudge_x       = -6 - subset(hazards, 
                               Gene.Symbol %in% gene_labeled_golden)$plog10P,
    nudge_y       = 1.3 - subset(hazards, 
                            Gene.Symbol %in% gene_labeled_golden)$uni_HR,
    segment.size  = 0.5,
    segment.color = "red",
    direction     = "y",
    hjust         = 0
  ) +
  
  # gene_labeled (other than CAMK2N1)
  geom_text_repel(data=subset(hazards, Gene.Symbol %in% gene_labeled),
                   aes(label=Gene.Symbol),
                   color = "grey50", #factor(golden)), #  == TRUE
#                   data          = subset(dat, wt < 3),
                   nudge_x       = 1 - subset(hazards, Gene.Symbol %in% gene_labeled)$plog10P,
#                   segment.color = "grey50",
                   direction     = "y",
                   hjust         = 1,
                   
                # x  color = ifelse(hazards$Gene.Symbol %in% gene_labeled_golden, "red", "grey50"),
#                   scale_color_manual(values = cols),
                                     # breaks = c("TRUE", "FALSE"),
                                     # guide = "none"), # better than color=
                   size=3, box.padding = unit(0.4, "lines"), 
                   segment.alpha=0.8, segment.size = 0.1,
                   max.overlaps = Inf)  + 
  theme_bw() + 
  theme(text=element_text(family="sans"),
        axis.title=element_text(size=18))
# thanks https://ggrepel.slowkow.com/articles/examples.html
###

# *** GSE65858 Cox HR (uni_HR) ####
# KM FDR-ajdusted P value: hazards534
load(file="HNSCC_GSE65858_FDRKMpvalue_7614probes.Rda") # hazards
hazards534 <- hazards[hazards$FDR.KMpvalue < 0.05, ]


ggplot(hazards534, aes(y = uni_HR, x = plog10P)) +   # x = plog10P  #Creation of ggplot graph with pre-determined aesthetics
  geom_point(shape = 21, alpha = 0.8, aes(size = log10P, fill = uni_HR)) + 
  scale_fill_distiller(palette = "RdYlGn", trans="log", limits = c(min(hazards534$uni_HR), max(hazards534$uni_HR))) + # min(hazards$uni_HR)
  #Diverging  BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
  geom_vline(xintercept = log10(cut_pvalue), color="red", lty=2) + 
  geom_hline(yintercept = 0.6, lty = 2) +
  geom_hline(yintercept = 1.8, lty = 2) +
  guides(size=FALSE, fill=FALSE) + 
  guides(size=FALSE, fill=FALSE) + 
  ylab("Cox Regression Hazard Ratio") + 
  xlab("log10(FDR-adjusted P value) of KM estimate") + # -log10(P-value)
  annotate(geom="text", x=-6, y=2, 
           label="HNSCC (GSE65858)",
           color="black", size = 5) +
    ##  geom_label_repel(data=subset(hazards, log10P >= 1.3 & as.numeric(multi_HR) >= 1.5 | log10P >= 1.3 & as.numeric(multi_HR) <= 0.5),
  #  geom_label_repel(data=subset(hazards, log10P >= -log10(cut_pvalue)),
  
  scale_x_continuous(
    expand = expansion(mult = 0.3)) +
  
  # for gene_labeled_golden (CAMK2N1, IL19, and FCGBP)
  # hazards534
    geom_label_repel(
    data=subset(hazards534, Gene.Symbol %in% gene_labeled_golden),
    #    data          = subset(dat, wt > 3),
    aes(label=Gene.Symbol),
    color = "red",
    size = 6,
    nudge_x       = -5 - subset(hazards534, 
                                Gene.Symbol %in% gene_labeled_golden)$plog10P,
    nudge_y       = 1.3 - subset(hazards534, 
                                 Gene.Symbol %in% gene_labeled_golden)$uni_HR,
    segment.size  = 0.5,
    segment.color = "red",
    direction     = "y",
    hjust         = 0
  ) +
  
  # gene_labeled (other than CAMK2N1)
  # hazards534
    geom_text_repel(data=subset(hazards534, Gene.Symbol %in% gene_labeled),
                  aes(label=Gene.Symbol),
                  color = "grey50", #factor(golden)), #  == TRUE
                  #                   data          = subset(dat, wt < 3),
                  nudge_x       = 1 - subset(hazards534, Gene.Symbol %in% gene_labeled)$plog10P,
                  #                   segment.color = "grey50",
                  direction     = "y",
                  hjust         = 1,
                  
                  # x  color = ifelse(hazards$Gene.Symbol %in% gene_labeled_golden, "red", "grey50"),
                  #                   scale_color_manual(values = cols),
                  # breaks = c("TRUE", "FALSE"),
                  # guide = "none"), # better than color=
                  size=4, box.padding = unit(0.4, "lines"), 
                  segment.alpha=0.8, segment.size = 0.1,
                  max.overlaps = Inf)  + 
  theme_bw() + 
  theme(text=element_text(family="sans"),
        axis.title=element_text(size=18))

###
# [2021/07/02] updated
save(hazards534, file="GSE65858_HNSCC_CoxHR_FDRKM.Rda")
ggsave("Rplot_GSE65858_CoxHR_CAMK2N1_top3FDRKM.pdf") #, width = 12, height = 8, dpi = 84)
# legands on Our findings suggested 20 candidate biomarkers, 

#hazards$Gene.Symbol[hazards$log10P>3]
#[1] "AGR2"   "CIDEB"  "CRB3"   "CRYM"  
#[5] "DUSP6"  "ERBB2"  "GGT6"   "JMJD8" 
#[9] "LYPD2"  "MSMB"   "MSMB"   "MSMB"  
#[13] "PEPD"   "PTGER3" "RBM11"  "SHMT1" 
#[17] "SNX14" 




# [2021/06/13]
#> Top3 CAMK2N1/IL19/FCGBP in TCGA ### 
# updated [2020/07/02]
# top3 CAMK2N1/CALML5/FCGBP in TCGA ####
#MSMB (22 genes/32 probes) in TCGA HNSCC (pvalueTex) ###
# all of them have TCGA FDR Pvalue just < 0.05 
# (MSMB > 0.03)(LYPD2 P = 0.000669 = log10 as 3.18)
load(file="HNSCC_GSE65858_FDRKMpvalue_7614probes.Rda") # as hazards
#hazards534 <- hazards[hazards$FDR.KMpvalue < 0.05, ]
#GSE65858_top <- hazards534[hazards534$Gene.Symbol %in%  gene_labeled, c(1,2,6,10)] 
# [1] "Gene.Symbol"  *(TCGA)"FDR.Pvalue"   *(TCGA)"multi_HR" 
hazard0 <- subset(hazards, select=c(Gene.Symbol, FDR.Pvalue, multi_HR, FDR.KMpvalue))
# col6: TCGA multi_HR
# col10: GSE65858 FDR-corrected KM P value
    
# (GSE65858)"FDR.KMpvalue"
#GSE65858_top22 <- GSE65858_top[!duplicated(GSE65858_top$Gene.Symbol), ] 
# 
#idx <- candidate_fdr6429$gene_id %in% GSE65858_top22$Gene.Symbol
#View(candidate_fdr6429[idx, ])

#  ggplot(Cox multi_HR of pvalueTex) ##
# the other volcano plot (KM Pvalue of GSE65858 vs Cox HR of PvalueTex)
# merge HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR + fdr100_KMsurvival (7614 probes)
# TCGA HNSCC Cox, multi_HR
#hazard0 <- subset(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, select=c(gene_id, multi_HR)) # p_value, p_value_adj_bonferroni, 
#colnames(hazard0)[1] <- "Gene.Symbol"

# ***colname will be:
# Gene.Symbol, TCGA (Cox multi_HR, FDR pvalue), GSE65858 (Cox uni_HR, FDR KMPvalue)
# adding GSE65858 (Cox uni_HR)
load(file="HNSCC_OS_GSE65858_pvalueCox_HR.Rda") # uni_HR => osHR
hazard_GSE65858 <- subset(osHR, select=c(Gene.Symbol, uni_HR)) # Cox.P-value 
hazards <- merge(hazard0, hazard_GSE65858, by = "Gene.Symbol")
hazards <- hazards[!duplicated(hazards), ]
#hazard0 <- hazards
# adding GSE65858 (FDR KM P value)
#load(file="HNSCC_GSE65858_FDRKMpvalue_7614probes.Rda") # as hazards

# R4> colnames(hazards)
#[1] "Gene.Symbol": TCGA *"FDR.Pvalue"   *"multi_HR"   //GSE65858  "FDR.KMpvalue" "uni_HR" 
hazards$log10P <- -log10(as.numeric(hazards$FDR.Pvalue)) # positive (for spot size)
hazards$plog10P <- log10(as.numeric(hazards$FDR.Pvalue)) # becoming negative (for X-axis)

# 5404 genes
hazards <- hazards[!duplicated(hazards$Gene.Symbol), ]

# to see how gene_labeled distributed in TCGA HNSCC world
gene_labeled_golden <- c("CAMK2N1", "CALML5", "FCGBP")
gene_labeled_silver <- c(
  "DKK1", #CAMK2N1, 
  "STC2", "PGK1", "SURF4", "USP10", "NDFIP1", "FOXA2", "STIP1", "DKC1", 
  "ZNF557", "ZNF266", "IL19", 
  "MYO1H", #FCGBP, 
  "LOC148709", "EVPLL", "PNMA5", "KIAA1683", "NPB"
)

# from pvalueTex without gene_labeled_golden, gene_labeled_silver
gene_labeled <- c(
  #  "DKK1", 
  #"CAMK2N1", 
  # "STC2", "PGK1", "SURF4", #"USP10", "NDFIP1",
  #"FOXA2", "STIP1", "DKC1", 
  #                  "ZNF557", "ZNF266", 
  #"IL19", "FCGBP",
  # "MYO1H", "LOC148709", "EVPLL", "PNMA5", "KIAA1683", 
  #"NPB",
  
  #    "ATP13A4",    
  #"PLAU",      
  #"ZNF662",     
  #"POMP",       
  #"DOT1L",     
  #"EPHX3" ,     #"IL34",       
  "FAM3D",     
  "MSMB", "LYPD2", "RBM11", # *** top KM P value in GSE65858 
  #"SERPINE1",   #"CALML5",    
  # "RASIP1",    # "FUT6",       
  #"BCAR3",     
  # "ST6GALNAC1", #"AQP1",   
  #"AIG1",        
  # "FAM151B"
  # "ABCB1",      
  #"SMPX"
  #"MASP1" ,   
  # "GRIA3"
  # CoxHR > 1.81:
  "ADAMTSL2",
  "BMP6" , "DUSP6" ,
  "FXYD5" ,   "GPX8",
  #"HERPUD2",
  "IFIT2" ,  
  #"KANK4",   
  "MYO1B" ,
  "NFKBID",   "PDIA5" ,   "PMEPA1",  "RAB43"  , 
  # "RET" ,  
  "SPP1" , 
  "TCP11L1",  "TIMM10"  ,
  "TNFAIP6" , # "TNPO1"  ,
  #"TPD52L2", 
  "TRAPPC10",
  "TTYH3" ,   "XYLT1" 
)

# **** or cv.mod.candidate from Cox-Lasso regression
# 16 genes
load(file="RLassoCox_TCGA_candidate62.Rda")
gene_labeled <- cv.mod.candidate$Gene.symbol
#
# *** TCGA HNSCC going for volcano plot by ggplot() ####
library(ggplot2)
library(ggrepel) # function "geom_label_repel"
cut_pvalue <- 0.05
ggplot(hazards, aes(y = multi_HR, x = plog10P)) +   # x = plog10P  #Creation of ggplot graph with pre-determined aesthetics
  geom_point(shape = 21, alpha = 0.8, aes(size = log10P, fill = multi_HR)) + 
  scale_fill_distiller(palette = "RdYlGn", trans="log", limits = c(min(hazards$multi_HR), max(hazards$multi_HR))) + # min(hazards$uni_HR)
  #Diverging  BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
  geom_vline(xintercept = log10(cut_pvalue), color="red", lty=2) + 
  geom_hline(yintercept = 0.6, lty = 2) +
  geom_hline(yintercept = 1.8, lty = 2) +
  guides(size=FALSE, fill=FALSE) + 
  guides(size=FALSE, fill=FALSE) + 
  ylab("Cox Regression Hazard Ratio") + 
  xlab("log10(FDR-adjusted P value) of KM estimate") + # -log10(P-value)
  annotate(geom="text", x=-5.0, y=2.3, #1.9, 
           label="HNSCC (TCGA)",
           color="black", size = 5) +
  ##  geom_label_repel(data=subset(hazards, log10P >= 1.3 & as.numeric(multi_HR) >= 1.5 | log10P >= 1.3 & as.numeric(multi_HR) <= 0.5),
  #  geom_label_repel(data=subset(hazards, log10P >= -log10(cut_pvalue)),
  
  scale_x_continuous(
    expand = expansion(mult = 0.1)) +
  
  # for gene_labeled_golden (CAMK2N1, IL19, and FCGBP)
  # hazards
  geom_label_repel(
    data=subset(hazards, Gene.Symbol %in% gene_labeled_golden),
    #    data          = subset(dat, wt > 3),
    aes(label=Gene.Symbol),
    color = "red",
    size = 6,
    nudge_x       = -5 - subset(hazards, 
                                Gene.Symbol %in% gene_labeled_golden)$plog10P,
    nudge_y       = 1.25 - subset(hazards, 
                                 Gene.Symbol %in% gene_labeled_golden)$multi_HR,
    segment.size  = 0.5,
    segment.color = "red",
    direction     = "both", #xy
    hjust         = 0,
    force_pull = 2 # pull label to it's data point
  ) +   ### ggplot() for gene_labeled_silver
  # hazards
  geom_label_repel(
    data=subset(hazards, Gene.Symbol %in% gene_labeled_silver),
    #    data          = subset(dat, wt > 3),
    aes(label=Gene.Symbol),
    color = "black",
    size = 3,
    nudge_x       = -3.2 - subset(hazards, 
                                Gene.Symbol %in% gene_labeled_silver)$plog10P,
    nudge_y       = 1.2 - subset(hazards, 
                                 Gene.Symbol %in% gene_labeled_silver)$multi_HR,
    segment.size  = 0.2,
    segment.color = "blue",
    direction     = "y",
    hjust         = 0,
    max.overlaps = Inf
  ) +

  # gene_labeled (other than CAMK2N1)
  # 
  # hazards
  geom_text_repel(data=subset(hazards, Gene.Symbol %in% gene_labeled),
                  aes(label=Gene.Symbol),
                  color = "grey50", #factor(golden)), #  == TRUE
                  #                   data          = subset(dat, wt < 3),
                  nudge_x       = 1 - subset(hazards, Gene.Symbol %in% gene_labeled)$plog10P,
                  #                   segment.color = "grey50",
                  direction     = "y",
                  hjust         = 1,
                  
                  # x  color = ifelse(hazards$Gene.Symbol %in% gene_labeled_golden, "red", "grey50"),
                  #                   scale_color_manual(values = cols),
                  # breaks = c("TRUE", "FALSE"),
                  # guide = "none"), # better than color=
                  size=3, box.padding = unit(0.4, "lines"), 
                  segment.alpha=0.8, segment.size = 0.1,
                  max.overlaps = Inf)  + 
  theme_bw() + 
  theme(text=element_text(family="sans"),
        axis.title=element_text(size=18))


##
ggsave("Rplot_TCGA_HNSCC_CoxHR_CAMK2N1_top3FDRKM.pdf") #, width = 12, height = 8, dpi = 84)



# [2021/06/15] Japan's AstraZeneca vaccine
# candidate_sample
hazard999 <- candidate_sample[complete.cases(candidate_sample),]
hazard999$fdr_KMPvalue <- p.adjust(hazard999$p_value, method="fdr")
# not one-hot variable
hazard999$dummy_rawPvalue <- ifelse(hazard999$p_value<0.001, "star3", "star1")
colnames(hazard999)[2] <- "Gene.symbol"

cox_HR_Pvalue <- data.frame("Gene.symbol"=matrix(whole_genome, ncol=1),
                            "multi_HR"=NaN, "multi_HRPvalue"=NaN)

#candidate_cox[[20500]][9, c(6,7)] # import Cox' multi_HR and it's Pvalue
for (i in 1:length(whole_genome)){
  tryCatch({
    cox_HR_Pvalue[i, c(2:3)] <- candidate_cox[[i]][9, c(6,7)]
  }, error = function(err) {
  # error handler picks up where error was generated
  #print(paste("MY_ERROR:  ", err))
  return()
  }) # END tryCatch
} # end of for loop i


hazard999 <- merge(hazard999, cox_HR_Pvalue, by="Gene.symbol")
View(hazard999[hazard999$multi_HRPvalue<0.05 & (hazard999$multi_HR>1.8 | hazard999$multi_HR<0.5), c(1,3,4,5)])
# n=64/106, they are raw Pvalue between 0.001 and 0.03

# plot of <0.001 (n=97) or 0.001-0.,03 (n=6332)
barplot(table(hazard999$dummy_rawPvalue))
#R4> colnames(hazard999) # from candidate_sample with KM+Cox
#[1] "Gene.symbol"     "number"          "p_value"        
#[4] "fdr_KMPvalue"    "dummy_rawPvalue" "multi_HR"       
#[7] "multi_HRPvalue"

require(dplyr)
library(tidyverse)
#ndf <- arrange(hazard999,number,p_value,multi_HR,multi_HRPvalue) %>%
num_breaks <- c(0.9, 1.1, 154, Inf) # max 166
# p_value, multi_HR
df999 <- hazard999 %>%
  group_by(number=cut(number, breaks=num_breaks))
hazard999$dummy_numberP <- df999$number

bar_num <- hazard999 %>%
  group_by(number=cut(number, breaks=num_breaks)) %>% # like mutate() 
#    group_by(bmi_cat = cut(mass/(height/100)^2, breaks=bmi_breaks)) %>%
  tally()
#1 (0.9,1.1]   354 => only one P value (median cut)
#2 (1.1,2.2]    211
#3 (2.2,3.3]    182
#4 (3.3,10.1]   819
#5 (10.1,Inf]  4863
#  summarise(num = n(), id = paste(Gene.symbol, collapse=","))
# number as 1-10: 1212; 10-166: 4863; NaN: 354
# total 6529 genes
barplot(bar_num$n, xlab="Number of Pvalues",
        names.arg=c("1", "2-154", "155-166"),
        main="Distribution of KM P values") #(hazard999$number)
plot(hazard999$dummy_numberP, -log10(hazard999$fdr_KMPvalue))
# number > 154 has significant KM pvalue

# to find cases_OS outside of range 180-220 (== cutoff not at median)
# pvalue plot as 倒V字形 "^"
# tar -xf HNSCC.survival.marginS.20500.tar.gz -T file20500_Rda_list.txt
# HNSCC_survivalAnalysis_marginS_*.Rda
for (i in c(1:length(whole_genome))){
  system(paste("echo 'HNSCC_survivalAnalysis_marginS_", whole_genome[i],".Rda'", 
               " >> file20500_Rda_list.txt" , sep=""))
}
# or sed -i '$HNSCC_survivalAnalysis_marginS_...' file20500_Rda_list.txt


# from here ->
# R4> which(whole_genome=="IL10RA") with M or ^ p_OS pvalues plots
# [1] 8333
OS_pvalueTex <- data.frame(matrix(NaN, nrow = 1, ncol = 3))
colnames(OS_pvalueTex) <- c("Gene.symbol", colnames(temp_OS_pvalue)[c(1,2)])

for (i in c(1:length(whole_genome))){
  fname <- paste("marginS/HNSCC_survivalAnalysis_marginS_", whole_genome[i],".Rda", sep="")
  if (file.exists(fname)) {
    #to  load(file=fname)
    # to get ojbects' name of load()
    temp.space <- new.env()
    saved_Rda <- load(file=fname, temp.space)
#    the.object <- get(saved_Rda, temp.space) # 1st ojbect
    if (nrow(temp.space[["OS_pvalue"]]) > 0){
      # save(list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue"), file=paste(gsub("_", "", marginTag), "/HNSCC_survivalAnalysis", marginTag, geneName, ".Rda", sep=""))
      # R2Excel export@TCGA_HNSCC_marginSFP.R 
      # OS_pvalue$p_OS >= 0.05 & cases_OS within 180-220
      temp_OS_pvalue <- temp.space[["OS_pvalue"]]
      temp_OS_pvalue <- temp_OS_pvalue %>% 
        mutate(pattern_v = (cases_OS >= 180) & (cases_OS <= 220))
      temp_OS_pvalueTex <-
          temp_OS_pvalue[(temp_OS_pvalue$p_OS >= 0.02) & (temp_OS_pvalue$pattern_v ==TRUE), c(1,2)]
      if (nrow(temp_OS_pvalueTex) > 0){ # might be more than 1 nrow

          onerow <- cbind(data.frame(matrix(c(whole_genome[i]), nrow=nrow(temp_OS_pvalueTex))),
                 temp_OS_pvalueTex)
          colnames(onerow)[1] <- colnames(OS_pvalueTex)[1]
          OS_pvalueTex <- rbind(OS_pvalueTex, onerow)
          # OS_pvalueTex appending
        print(tail(OS_pvalueTex, n=1))
      } # end if; nothing found
          
      
    } # else there is no OS_pvalue
    rm(temp.space)
  } # end of if file.exists
} # end of for loop i

View(OS_pvalueTex)
# 99 genes have one p_OS >= 0.05, in OS_pvalueTex005
OS_pvalueTex002 <- OS_pvalueTex
# 3641 genes have one p_OS >= 0.02, in OS_pvalueTex002
# special case: HNSCC_survivalAnalysis_marginS_NDFIP1.xlsx
# NDFIP1 (within our 20 candidatte biomarkers), median cut will let NDFIP1 to be kik-out prematurely
# show it's pvalue plot
View(OS_pvalueTex002[OS_pvalueTex002$Gene.symbol=="NDFIP1",])
OS_pvalueTex_NDFIP1 <- OS_pvalueTex002[OS_pvalueTex002$Gene.symbol=="NDFIP1",]
plot(OS_pvalueTex_NDFIP1$cases_OS, OS_pvalueTex_NDFIP1$p_OS)
# => supplementary figure as Rplot_pvaluePlot_NDFIP1.pdf
###


# CALML5 as new candidate ####
# [2021/07/02]
save(consensus_9, top_TCGA, top_gse65858, file="consensus_CAMK2N1_CALML5_FCGBP.Rda")


#[2021/06/27]
# ***TCGA median cut KM plot ####
# how many genes with P value < 0.05?
# from pvalueTex
# > fit KM survival curves ####
# cancer type should be defined at TCGA_cohort <- "LUAD" "HNSC"
#x GEO dataset: GSE65858, dim 256 * 7620
#x load(file="gset_OS_mRNA2.Rda") # as gset_OS_mRNA2
#
# load TCGA HNSCC dataset
# xload(file="HNSCC.clinical.RNAseq_tobacco.Fire.Rda") #clinical data only
load(file="HNSCC.clinical.RNAseq.Fire.Rda") # 20500 RNA-Seq dataset
#oscc <- clean6_oscc_tobacco
#oscc0 <- merge(oscc, clean6_oscc_tobacco, by "")

# "status" "time"  n=521
survData <- clean6_oscc_tobacco[, c(10,12)] # "Unique.ID","OS_IND" "OS_months"
colnames(survData) <- c("status", "time")
survData$status <- ifelse((survData$status==1), TRUE, FALSE)
rownames(survData) <- clean6_oscc_tobacco[, c(1)] # ID

# Gene.symbol is colname: "z.score_A1BG"...p=20239 genes
mRNA_matrix <- clean6_oscc_tobacco[, c(16:20254)]
rownames(mRNA_matrix) <- clean6_oscc_tobacco[, c(1)] # ID

## connecting pvalueTex code to GEO code:####
# *** merge survData, mRNA_matrix into gset_OS_mRNA2, np521, p=20254
#df <- as.data.frame(matrix(NaN, ncol = 1, nrow = nrow(gset_OS_mRNA2))) # Don't skip 1st row 
df <- as.data.frame(row.names(survData))
colnames(df)[1] <- "Gene.symbol"
df <- cbind(df, survData, mRNA_matrix)
# OS time (days), OS event
colnames(df)[c(3,2)] <- c("OS_time", "OS_event")
library(stringr) # extract last word
# # "z.score_" removal by str_replace
#df$OS_time <- str_extract(gset_OS_mRNA2$OS_time,"\\w+$")
gsymbol <- str_replace(colnames(df)[4:20242], "z.score_", "")
colnames(df)[4:20242] <- gsymbol
#df$OS_event <- str_extract(gset_OS_mRNA2$OS_event,"\\w+$") #first word
# if # extract first word: "\\w+" in regular expression

df$OS_time <- as.numeric(df$OS_time) # months
df$OS_event <- ifelse(df$OS_event == "TRUE", 1, 0) # OS event TRUE as 1 (dead)
# 1: 88 vs 0:168
#gset_OS_mRNA2[, c("OS_time", "OS_event")] <- df[, c("OS_time", "OS_event")]
gset_OS_mRNA2 <- df

# probes number: coln122-coln7+1 =116 (100 genes)
#library('taRifx')
#gset_OS_mRNA3 <- remove.factors(gset_OS_mRNA2)
# unlist or lapply if legnth > 1 to coerce
#gset_OS_mRNA3[, c(coln7:coln122)] <- lapply(gset_OS_mRNA3[, c(coln7:coln122)], as.numeric)
#gset_OS_mRNA3 <- as.numeric(gset_OS_mRNA3[, c(coln7:coln122)])
#gset_OS_mRNA2 <- gset_OS_mRNA3
save(gset_OS_mRNA2, file="TCGA_KM_gset_OS_mRNA2.Rda") # TCGA pvalueTex [2021/06/27]

# data cleaning####
load(file="TCGA_KM_gset_OS_mRNA2.Rda")
# as gset_OS_mRNA2
#removal missing event status or time-to-event (2 NaN  cases at no 165 416)
table(complete.cases(gset_OS_mRNA2))
summary(gset_OS_mRNA2[,1:3])
which(is.na(gset_OS_mRNA2$OS_time))
gset_OS_mRNA2 <- gset_OS_mRNA2[complete.cases(gset_OS_mRNA2[,1:3]), ]
# n=519 now
# check missing of gene expression
# (159 missing)
col_nan <- which(is.na(gset_OS_mRNA2[1,]))
gset_OS_mRNA2 <- subset(gset_OS_mRNA2, select = -col_nan)
# 20083= 20242-159
# dim 519 * 20083 now
# knockout TCGA-BA-A6DF (at 32)
gset_OS_mRNA2 <- gset_OS_mRNA2[-32,]
# # dim 519 * 20083 now
# still NaN at somewhere (n=225) => check NaN expression over 30% of each gene => skip (during for i loop)
# If you need NA count Column wise – sapply(z, function(x) sum(is.na 3.6k(x)))
#
# x
# knock-out at 165 416: TCGA-CQ-A4CA, TCGA-H7-A6C4
#mRNA_matrix[c(165, 416), 1] <- NaN
# and at 32: TCGA-BA-A6DF
#survData[c(32), 2] <- NaN
#survData <- survData[complete.cases(survData$time), ]
# missing at 32 x (imputation this missing data)
#mRNA_matrix <- mRNA_matrix[complete.cases(mRNA_matrix$"1"), ]

# run from here => go ####
# *************
# removal of all generated pdf
library(survival)
# 
mysurv <- Surv(gset_OS_mRNA2$OS_time, gset_OS_mRNA2$OS_event == 1) #1==dead, time by months

# check MASP1, warnings() on more than one probe to a gene
#which("MASP1"== colnames(gset_OS_mRNA2))
# [1] 1503 1580 1928
#
# for loop (20242-4+1=20239 genes) save KM pvalues
fdr100_KMsurvival <- data.frame("Gene Symbol"=NaN, "FDR Pvalue"=NaN, 
                                "KM Pvalue"=NaN, "cases_High"=NaN, "cases_Low"=NaN)
# a empty 1 x 5 data.frame
#i <- 4 # for TCGA
#coln7
coln4 <- 4 # exp of gene ZNF589 -> A1BG
# error at 16114: data set has no non-missing observations
coln4 <- 16114 +1
coln4 <- 18326+1 # skip "TRYX3": all are -0.04203
coln122 <- ncol(gset_OS_mRNA2) # 2298 now for FDR2000; 20242 -> 20083 for TRCGA clean6_oscc_tobacco of pvalueTex
n_cases <- nrow(gset_OS_mRNA2)
#
for (i in c(coln4:coln122)) {
  # expression cutoff at median value -> grouping by dichomolized it
  # cutoff <- sapply(gset_OS_mRNA2[, i], median)
  #  still NaN at somewhere (n=226) => check NaN expression over 30% of each gene => skip
  if (sum(is.na(gset_OS_mRNA2[, i]))/n_cases < 0.3) {
  cutoff <- median(gset_OS_mRNA2[, i], na.rm = TRUE)
  # dichotomize
  #oscc$ageDx[gset_OS_mRNA2[, i] < cutoff] <- 0 # underexpression
  #oscc$ageDx[gset_OS_mRNA2[, i] >= cutoff] <- 1 # overexpression than cutoff
  #
  # *** one group issue: survdiff.fit
  # A1CF 13% cases >1, while 87% cases = -0.30949; so median is -0.30949; 
#  % all cases belong to >= median
# solved by assign: [, i] > cutoff
  gset_OS_overexpression <- ifelse((gset_OS_mRNA2[, i] > cutoff), 1, 0)
  # imputation NaN as "< cutoff" group
  # *** check "all are the same"
  if (sum(gset_OS_overexpression) == 0)  {
    print(paste("[xx ", i, "/", coln122, "] skip due to 'all are the same' of", colnames(gset_OS_mRNA2)[i]))
    next # skip to next i;
    # not break 
  }
  # Test for difference (log-rank test) between groups (by PMM1_median 0 vs 1)
  #tryCatch(surv_OS1 <- survdiff(mysurv ~ as.vector(gset_OS_mRNA2[, gset_OS_mRNA2M_pos], mode="numeric"), data=gset_OS_mRNA2), error = function(e) return(NA)) # PMM1 high or low
  # or OS.km ; as.vector(unlist(gset_OS_mRNA2[, osccCleanNAM_pos]), mode="numeric"), data=gset_OS_mRNA2
  # survfit(, type= "kaplan-meier", conf.type = "log-log") 
  # for fit, plot; 
  surv_OS1 <- survfit(mysurv ~ gset_OS_overexpression, type= "kaplan-meier", conf.type = "log-log") # log-rank test
  
  # survdiff for two group with P value
  # pchisq gives the distribution function
  surv_OSp <- survdiff(mysurv ~ gset_OS_overexpression) # log-rank test
  # N Observed Expected (O-E)^2/E (O-E)^2/V
  # broom::glance(surv_OS1)$p.value
  p_OS1 <- format(pchisq(surv_OSp$chisq, length(surv_OSp$n)-1, lower.tail = FALSE), digits=3)
  # (p_OS1 == p_OS0) is TRUE
  #cases_OS1 <- surv_OSp$n[1]
  # #Value of survdiff => a list with components: help("survdiff")
  # n => the number of subjects in each group.
  # obs
  # the weighted observed number of events in each group. If there are strata, this will be a matrix with one column per stratum.
  # exp
  # the weighted expected number of events in each group. If there are strata, this will be a matrix with one column per stratum.
  # chisq
  # the chisquare statistic for a test of equality.
  # var
  # the variance matrix of the test.
  # strata
  # optionally, the number of subjects contained in each stratum.
  
  # [survfit] - Kaplan-Meier curve: P-Value ####
  # confidence intervals as log hazard or log(-log(survival))
  #OS.km <- survfit(mysurv ~ as.vector(unlist(gset_OS_mRNA2[, gset_OS_mRNA2M_pos]), mode="numeric"), data=gset_OS_mRNA2, type= "kaplan-meier", conf.type = "log-log")
  # 365.25 days in a year for xscale => 3650 days for 10 years
  # maximal followup to 6.8 years => 2483.7 days
  # 12 months per year for 5 years => 60 months, 1826.25 days
  # xscale=12, xmax=60
  
  pdf_fn <- paste("rplot_TCGA_", colnames(gset_OS_mRNA2)[i], ".pdf", sep="")
  if (file.exists(pdf_fn)) {
    pdf_fn <- 
      paste("rplot_TCGA_", colnames(gset_OS_mRNA2)[i], "_", i, ".pdf", sep="")
  }
  pdf(pdf_fn) # gene symbol
  # save PDF, size 746 x 431 pixel(?)
  # xmax by months; if xscale by days 365.25
  plot(surv_OS1, lty=1, xscale=12, xmax=200, col=c("blue","red"), 
       sub=paste("Kaplan-Meier P Value =", p_OS1), 
       main=paste("OS in TCGA", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", colnames(gset_OS_mRNA2)[i]), ylab="Percent Survival", xlab="Years")
  legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:1, col=c("blue","red"))
  dev.off()
  #plot(OS.km, lty=1, col=c("blue","red"), sub="p= 0.816", main="OS in gset_OS_mRNA2(n=505)/gene level", ylab="Percent Survival", xlab="Days")
  #legend("topright", legend=c('Low', 'High'), lty=1:2, col=c("blue","red"))
  # summary(OS.km, times = seq(0, 3000, 100))
  # 
  # save/store P value
  # [1] "Gene.Symbol" "FDR.Pvalue"  "KM.Pvalue"   "cases_High"  "cases_Low" 
  fdrpvalue <- 2 # false value (we calculate FDR p value after for loop i) #fdrpvalue <- candidate_fdr100$p_value[which(candidate_fdr100$Gene.symbol==colnames(gset_OS_mRNA2)[i])]
#  if (length(fdrpvalue) > 1) {fdrpvalue <- fdrpvalue[1]} # prevent duplicated genes (FDR Pvalue)
  fdr100_KMsurvival[i-which(colnames(gset_OS_mRNA2) == "OS_event"), ] <- 
    c(colnames(gset_OS_mRNA2)[i], 
      fdrpvalue, 
      p_OS1, surv_OS1$n[2], surv_OS1$n[1])
  print(paste("[", i, "/", coln122, "] Generating median-cut KM plot of", colnames(gset_OS_mRNA2)[i]))
  
  #"(", which(whole_genome==geneName), ")"))
} # end of if 30% NaN
} # end of for loop

View(fdr100_KMsurvival)
# p=20080
# removal of 1st NaN row
fdr100_KMsurvival <- fdr100_KMsurvival[-1,]
# 3117 genes has unadjusted KM pvalue < 0.05
# FDR calculation
fdr100_KMsurvival$FDR.Pvalue <- p.adjust(fdr100_KMsurvival$KM.Pvalue, method = "fdr")
# or
fdr100_KMsurvival$FDR.Pvalue <- p.adjust(fdr100_KMsurvival$KM.Pvalue[fdr100_KMsurvival$KM.Pvalue < 0.05], method = "fdr")
# 12/100 or 152/2000 discovery rate by GSE65858
# 530/7614
plot(fdr100_KMsurvival$FDR.Pvalue[fdr100_KMsurvival$KM.Pvalue < 0.05], fdr100_KMsurvival$KM.Pvalue[fdr100_KMsurvival$KM.Pvalue < 0.05], log = "y", xlab="KM FDR P value", ylab="KM unadjusted P value", main="TCGA HNSCC Kaplan-Meier estimate with cutoff at median expression")
abline(h=0.05)
abline(v=0.05)

save(fdr100_KMsurvival, file="TCGA_HNSCC_KMsurvival_medianCut20080.Rda")
# 3118 genes
View(fdr100_KMsurvival[fdr100_KMsurvival$KM.Pvalue<0.05 ,])
# 209 genes
View(fdr100_KMsurvival[fdr100_KMsurvival$FDR.Pvalue<0.05 ,])

# ****** IL19 必須要靠 pvalueTex 才能發現
> gene_labeled_golden[!gene_labeled_golden %in% fdr100_KMsurvival$Gene.Symbol[fdr100_KMsurvival$FDR.Pvalue<0.05]]
#[1] "IL19"; median cut at 258 vs 260
# KM pvalue = 0.0038
# FDR Pvalue = 0.1154315
# when [optimal cutoff] by pvalueTex
# KM pvalue = 3.73e-7
# FDR Pvalue = 6.54e-6
# % scp  tex@35.201.169.0:~/R/HNSCC_Tex_survival/hnscc_github/rplot_GSE65858.tar.gz ./Downloads
# end of major revision and answer (1-1 to 2-6)




### (spared code) ####
###### Download platform data from GEO and get sample (phenotype) information ## 
GPL = "GPL10558"
data.platform = getGEO(GPL)
data.index = match(GPL, sapply(gset, annotation))
data.p = pData(gset[[data.index]])

data.expr = exprs(gset[[data.index]])
common = intersect(colnames(data.expr), rownames(data.p))
m1 = match(common, colnames(data.expr))
m2 = match(common, rownames(data.p))
data.expr = data.expr[,m1]
data.p = data.p[m2,]

## generate boxplot of expression profiles ##
title = "samples" 
s.num = 1:ncol(data.expr)
n = ncol(data.expr)
if (n > 30) {
  s.num = sample(1:n, 30)
  title = "selected samples"
}
title = paste0(GSE, "/", GPL, " ", title)

fixed.df <- as.data.frame(x=data.expr[,s.num], stringsAsFactors = FALSE)

x1 <- reshape2::melt(fixed.df, na.rm = TRUE, id.vars = NULL,
                     variable.name = "variable", value.name = "value")

exp.prof.plot <- ggplot(x1, aes(variable, value)) +
  geom_boxplot(outlier.colour = "green") +
  labs(title = title, y = "log2 expression", x = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(exp.prof.plot)
#########%###

# > fit linear model ####
fit <- lmFit(gset, design)  # package: limma

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "FDR-adj P value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)



#############################################################2###
# > General expression data analysis####
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE65858", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE65858", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE65858")


### shinyGEO ####
### for Kaplan-Meier survivial analsys
### Availability and Implementation: Web application and source code are available from
# download from https://github.com/gdancik/shinyGEO
# x https://gdancik.shinyapps.io/shinyGEO/
runApp()
# R code for survival
# https://github.com/gdancik/shinyGEO/blob/master/server/server-survival.R
# 
shlibrary(shiny)

ui <- fluidPage(
  
)

server <- function(input, output, session) {
  
}

shinyApp(ui, server)

