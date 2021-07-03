# [Results]  ####
#= Analysis of output .Rda of _marginS_ on GCE i-4 or _marginFree_ on GCE i-4Free
# or i-4Plus
# >*** [choice ONE]: _marginFree_ or _marginS_ loading from .Rda
#marginTag <- "_marginS_" #at ./marginS
#marginTag <- "_marginFree_" #at ./marginFree
#marginTag <- "_marginPlus_" #at ./marginPlus
TCGA_cohort <- "HNSCC" # start-over again [2020/03/08]
raw <- 
 readline("_margin(S)_, _marginFree(F)_ or _margin(P)lus_ -- process run: ")
select_margin <- function(x) {
  switch(x,
         s = "_marginS_",
         f = "_marginFree_",
         p = "_marginPlus_",
         S = "_marginS_",
         #         F = "_marginFree_",
         P = "_marginPlus_",
         stop("Unknown input")
  )
}
marginTag <- select_margin(raw)
print(marginTag)
rm(raw)
# [2019/07/03-2019/07/22] they are stored at ./xlsx => moved to ./marginS
# [2019/07/03-2019/07/22] they are stored at ./xlsx => moved to ./marginFree
# [2019/07/03-2019/07/22] they are stored at ./xlsx => moved to ./marginPlus
# [2019/07/30-2019/08/09] they are stored at ./xlsx => moved to ./marginS/tobacco

# set path to the source of .Rda from cutoff finding
path_cohort <- getwd() #  "/home/tex/R/HNSCC_Tex_survival/hnscc_github"
path_ZSWIM2 <- file.path(path_cohort, gsub("_", "", marginTag), "") # e.x. marginS/

# ZSWIM2_archive1000_20180408_0042_0933.Rda; 9 hours for 1,000 genes to be scanned
# 3 hours for 1,000 genes to be scanned under GCP Rstudio server
# STRAP to SLC8A1
#which(whole_genome==geneName)
#
#_marginS_
# x merge ZSWIM2a ZSWIM2b and ZSWIM2c => ZSWIM2abc_archive.Rda [2018/06/19]
# para_seg <- c(20499, 13528, 6833, 1) # ZZZ3-> PLCE1, PNMT-> GARNL3, GAR1-> A1BG
#load(file=file.path(path_cohort, "ZSWIM2_archive.Rda")) # as ZSWIM2
#ZSWIM2 <- ZSWIM2abc #(done, merged)
#save(data=ZSWIM2, file=file.path(path_ZSWIM2, "ZSWIM2_archive.Rda")) # or ZSWIM2abc_archive 
#x OR
# _marginFree_
#load(file="ZSWIM2_free_archive.Rda") # as well as ZSWIM2



#(at the first time) 
#>Dissection of ZSWIM2 ####
#1) KM plot is not at best cutoff, why??
#2) ZSWIM2 should be saved by appending.
#3)
load(file=file.path(path_cohort, "ZSWIM2_archive.Rda")) # x merged abc
ZSWIM2$X3 <- addNA(ZSWIM2$X2, ifany=F)
ZSWIM2_2 <- table(ZSWIM2$X3) # try to tryCatch (1, 2)
# 0     1     2     3     4     5     6     7     8     9 
# 682    64   187    26    23    21     9    18    13    14 
# 10    11    12    13    14    15    16    17    18    19 
# 14    10    12    10     5     6     5     6     2     5 
# 20    21    22    23    24    25    26    27    28    29 
# 3     5     1     4     2     3     1     2     5     2 
# 30    31    34    35    36    37    40    41    44    49 
# 3     3     3     2     2     2     1     1     1     1 
# 51    52    57    60    66  <NA> 
#   1     2     1     1     1 19314 
ZSWIM2_2 <- data.frame(ZSWIM2_2)
#ZSWIM2_2$Freq[46] <- ZSWIM2_2$Freq[46] - (17406+(LUAD_n-18602)) #number of NA:11, from TRIL(18602) to "STX10"(17406)
#print(run16h <- (sum(ZSWIM2_2$Freq)-(LUAD_n-(18602-17406))))
# running 18602-17406=1196 genes in 16 hours
plot((ZSWIM2_2))
sum(ZSWIM2_2$Freq)
#[1] 20500  # scan completely
# why error_01: whole_genome[ZSWIM2$X2 == 1]
# 


# Only reture(0) is ok for subsequent analyses
pielabels <- sprintf("%s", ZSWIM2_2$Freq[ZSWIM2_2$Freq>0])
pie_colors <- c("yellow", "green", "violet", 
              "orange", "blue", "pink", "cyan","red") 
text_pie <- paste("Workable genes \n in", pie_colors[1], "(", round(ZSWIM2_2$Freq[1]/sum(ZSWIM2_2$Freq)*100), "% of whole genome)") # _marginS_ about 45.9% usable RNAseq data
# n=9416
pie(ZSWIM2_2$Freq[ZSWIM2_2$Freq>0],
    labels=NA,
    clockwise=F,
    col=pie_colors,
    border="white",
    radius=0.7,
    cex=0.8,
    main=text_pie)
legend("bottomright", legend=pielabels, bty="n",
       fill=pie_colors)
# brewer.pal(7,"Set1") # from RColorBrewer::brewer.pal()
#browsers<-read.table("browsers.txt",header=TRUE)
#browsers<-browsers[order(browsers[,2]),]
#pielabels <- sprintf("%s = %3.1f%s", browsers[,1],
#                     100*browsers[,2]/sum(browsers[,2]), "%")
#pie(ZSWIM2_2$Freq, labels=ZSWIM2_2$Var1, main=text_pie)

#
# deal with return(1), return(2) and return(3) errors, try to solve it
# _marginS_
n_percent_Bonferroni <- ZSWIM2_2$Freq[1]/sum(ZSWIM2_2$Freq) # 0.4593171 -> 0.4581951
error01_sample <- which(ZSWIM2$X3==1) #  32.2% error01: "ZSWIM2 skip from contingencyTCGA()"): one group issue in Melt()
# => *** re-run
error02_sample <- which(ZSWIM2$X3==2) # 14.2% error02: There has only one group in survdiff.
error03_sample <- which(ZSWIM2$X3==3) # 6.85% error03: There has only one group in survdiff on KM cutofFinder_func.R.
#=> Test for difference (log-rank test is multiple Chi-square alone survival time) with P-value; or Test for different KM curve between two groups (seperated by PMM1_median 0 vs 1)
error05_sample <- which(ZSWIM2$X3==5) # 0.78%(159 genes) error05: for why?
# done #





# #debug {CoxPH table 3, univariate 全部是同一個值 from gender to tobacco} only RNAseq 是正確的 [2019/08/14]
# >A summary table of P-value, Z-score standardization from all HNSCC_survival*.Rda ####
# 不必每次重 run (saved at HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda)
# go ahead post2, a simple way for p.adjut() =>
# _marginFree_ or _marginS_ or _marginPlus_ from .Rda
# number of OS_pvalue: from cutoff finder, the number (frequency) of OS P-values, which KM P-value < 0.05, in each gene;


# get_Rda_pvalue <- function(geneName) {
#   load(file=paste("HNSCC_survivalAnalysis_marginS_", geneName, ".Rda", sep=""))
#   # load list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue")
#   # a example: geneName <- "TRIP13"
#   #OS_pvalue$p_OS[which.min(OS_pvalue$p_OS)]
#   
#   if (nrow(OS_pvalue) > 0)  { # we hit this gene with P-value < 0.05 in KM plot
#     return(min(OS_pvalue$p_OS))} else {return(NA)}
# }

# candidate_sample(KM) and candidate_cox (Cox list or Cox table), two parts
# candidate_sample: Kaplan–Meier plots for genes significantly associated with survival.
# Table1: list genes, which it's KM plot with P-value < 0.05, and ranking by P-value (ascending) and z-score (decending)

aa <- LUAD_n; bb<- 1
#aa <- ZSWIM2_2$Freq[1]  # n=9416
# however, $ ls HNSCC_sur*.Rda | wc -l
# 16017

#{ declare empty data.frame and list
#*
candidate_sample <- data.frame(matrix(data = NA, nrow = aa, ncol = 3)) # for num, gene ID and it's P-value
colnames(candidate_sample) <- c("number", "gene_id", "p_value") # KM OS P-value only (there is no DFS in HNSCC); 
#x"number" position in ZSWIM2
#xcandidate_sample$number <- data.frame(which(ZSWIM2$X3==0)) # gene "number"
candidate_sample$gene_id <- whole_genome[bb:aa] # retrieving gene name from whole_genome; return() at X2
#

#* a list for all cox survival datas
# http://www.cookbook-r.com/Manipulating_data/Converting_between_data_frames_and_contingency_tables/
# column name might be created by variable: e.x. df.sex[,"per"] <- df.sex$count1/sum(df.sex$count1); # df.sex$per
candidate_cox <- replicate(aa, list()) # empty list(), which doesn't need colnames
#data.table(matrix(data = NA, nrow = aa, ncol = 3)) # for num, gene ID and it's P-value
# colnames(candidate_cox) <- c( "Features", "HR",       "P_value_uni",  "sig",      "Features", "HR",      
#                               "P_value_multi",  "sig",      "Features", "P_value_KM",  "sig" )
# #}


setwd(paste(path_ZSWIM2, "Rda", sep="")) # set for a while (for ip loop) at ./marginS/Rda
# * all .Rda 要先解壓縮 HNSCC.survival.marginS.20500.tar.gz
# retrieved then saved as HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda
for (ip in (bb:aa)) { # 3 hours for each run
  geneName <- candidate_sample$gene_id[ip]
  #print(paste("At", path_ZSWIM2, "=> (", ip, ")", geneName), sep="")
  #candidate_sample$p_value <- lapply(unlist(candidate_sample$gene_id[bb:aa]), get_Rda_pvalue) # retrieving P-value by gene name from .Rda; return() at X3
  #[choice ONE]: _marginFree_ or _marginS_ loading from .Rda
  #  _marginS_
  # aka: "HNSCC_survivalAnalysis_marginS_DKK1.Rda"
  # list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue")
  load_filename <- file.path(paste(TCGA_cohort, "_survivalAnalysis", marginTag, geneName, ".Rda", sep=""))
  #OR (automatically defined)
  #_marginFree_
  #load_filename <- paste("HNSCC_survivalAnalysis_marginS_", geneName, ".Rda", sep="")
  #
  #
  if (load_filename %in% dir()) {
    
    load(file=load_filename)
    #      if (is.na(tryCatch(load(file=load_filename), error = function(e) return(NA)))) {} # load file with error free :-)
    #  load(file=paste("HNSCC_survivalAnalysis_marginS_", geneName, ".Rda", sep=""))
    # load list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue")
    # a example: geneName <- "TRIP13"
    print(paste(marginTag, ", there is ", nrow(OS_pvalue), " OS P-values (uncorrected)"))
    # create a temp variable: OS_pvalueS
    OS_pvalueS <- as.data.frame(matrix(NA, nrow=nrow(OS_pvalue), ncol=3))
    colnames(OS_pvalueS) <- c("p_value", "p_value_adj_bonferroni", "p_value_adj_FDR")
    OS_pvalueS$p_value <- OS_pvalue$p_OS
    OS_pvalueS$p_value_adj_bonferroni <- p.adjust(OS_pvalue$p_OS, method='bonferroni')
    OS_pvalueS$p_value_adj_FDR <- p.adjust(OS_pvalue$p_OS, method='fdr')
    if ((nrow(OS_pvalue) > 0) & any(OS_pvalueS$p_value_adj_FDR < 0.05)) {
      
      candidate_sample$p_value[ip] <- min(OS_pvalueS$p_value_adj_FDR) # most sigificant P-value => 要改為 min(FDR P-value)
      # for 20500 genes comparison, candidate_sample$p_value 仍須要再一次 Bonferroni correction
      candidate_sample$number[ip] <- nrow(OS_pvalueS[OS_pvalueS$p_value_adj_FDR<0.05, ])
      #number of OS_pvalue: from cutoff finder, the number (frequency) of OS P-values, which KM P-value < 0.05, in each gene;
    }
    #candidate_sample$p_value[ip] <- get_Rda_pvalue(candidate_sample$gene_id[ip]) # retrieving P-value by gene name from .Rda; return() at X3
    #  print(paste(ip, geneName, ": ", candidate_sample$p_value[ip], sep=" "))
    
    #***
    # Significant features from table 2, table 3 and table 4 of output .Rda files ##
    #{
    # [uni_cox_pvalue and uni_HR]
    # [multi_cox_pvalue and multi_HR]
    # [exp_pvalue]
    # from tableChi1 (Table 2), tableOS1 (Table 3) /[tableRFS1]
    
    # focusing on tableOS1, Cox proportiopnal hazard model: (both univariate and multivariate) P-value <= 0.05,
    # and sorted by HR > 1 vs HR < 1
    # P-value notes:
    # Significant codes:  0 ???***??? 0.001 ???**??? 0.01 ???*??? 0.05 ???.??? 0.1 ??? ??? 1
    # if <0.001 => mark as "***"
    # osHR$X4[osHR$X4<0.001] <- "***"
    
    # uni_cox_pvalue
    uni_cox_pvalue <- tableOS1[c(FALSE, TRUE), c(1,2,5)] # get significant P-value, hazard ratio at column 2
    uni_cox_pvalue$`P-value` <- as.character(uni_cox_pvalue$`P-value`)
    uni_cox_pvalue$`P-value`[(uni_cox_pvalue$`P-value` == "***")] <- "0.001"
    uni_cox_pvalue$`P-value` <- as.numeric(uni_cox_pvalue$`P-value`)
    uni_cox_pvalue$sig <- uni_cox_pvalue$`P-value` <= 0.05 # "sig" marking for significant
    
    # multi_cox_pvalue
    multi_cox_pvalue <- tableOS1[c(FALSE, TRUE), c(1,6,9)] # get significant P-value, hazard ratio at column 6
    multi_cox_pvalue$`P-value` <- as.character(multi_cox_pvalue$`P-value`)
    multi_cox_pvalue$`P-value`[(multi_cox_pvalue$`P-value` == "***")] <- "0.001"
    multi_cox_pvalue$`P-value` <- as.numeric(multi_cox_pvalue$`P-value`)
    multi_cox_pvalue$sig <- multi_cox_pvalue$`P-value` <= 0.05 # "sig" marking for significant
    
    # *** we don't have RFS in HNSCC cohort :-)
    
    # and tableChi1: P-value <= 0.05 # exp_pvalue
    # "sig" marking the significant "features": check "odd" position
    exp_pvalue <- tableChi1[c(TRUE, FALSE), c(1,7)] # get significant P-value from column 7, name from column 1; remark at column 8
    exp_pvalue$`P-value` <- as.numeric(as.character(exp_pvalue$`P-value`))
    exp_pvalue$sig <- exp_pvalue$`P-value` <= 0.05 # "sig" marking for significant
    
    # # merging them as one by common row names
    # > colnames(exp_pvalue)
    exp1 <- data.frame("X1"= geneName, "X2" = NA,  "X3"=NA) # append one row, with this geneName
    colnames(exp1) <- colnames(exp_pvalue) # precisely matched
    candidate_cox_ip <- cbind(uni_cox_pvalue, multi_cox_pvalue, rbind(exp_pvalue, exp1)) # colname is duplicated, bind 3 tables together
    colnames(candidate_cox_ip) <- c( "uni_Features", "uni_HR",       "uni_P_value",  "uni_sig",      "multi_Features", "multi_HR",
                                     "multi_P_value",  "multi_sig",      "KM_Features", "KM_P_value",  "KM_sig") # KM_sig is Remark at table 2(tableChi1)
    
    # it must be [[]] for assign the content !!! https://stat.ethz.ch/R-manual/R-devel/library/base/html/Extract.html.
    candidate_cox[[ip]] <- candidate_cox_ip #http://cran.r-project.org/doc/manuals/R-lang.html#Indexing
    print(paste(ip, geneName, ": cox features saved; ncol = ", length(candidate_cox[[ip]]))) # a vector in [] is ok for indexing
    # length == 11 => that is correct !
    #}
  }
} # end of ip for loop
# #[2019/06/26] _marginS_ finished at 19:48 (about 3 hours)
# #[2020/03/12] _marginS_ finished at 01:00 (about 3 hours) with Bonferroni and FDR P-value correction
# n=6429

#_marginS_ or _marginFree_ by SFree; saving on ./run04_marginS_, files => 17030
save(candidate_sample, candidate_cox, n_percent_Bonferroni, 
     file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep=""))) #ok; with KM_sig and Remark, and Cox HR
setwd(path_cohort) 
#_marginS_ by marginTag
# saved file="HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda" above
# R4> candidate_sample[4:5,] # yes with FDR correction # n=6429
# number gene_id    p_value
# 4     14   A2LD1 0.04930000
# 5     99   A2ML1 0.04648171
# old version preserved as HNSCC_OS_marginS_pvalueKM_candidate_cox_bak_cutoffnoFDR.Rda
# R4> candidate_sample[4:5,] # noFDR # n=6429
# number gene_id p_value
# 4     14   A2LD1 0.03210
# 5     99   A2ML1 0.00552


# devtools::install_github(c('jeroenooms/jsonlite', 'rstudio/shiny', 'ramnathv/htmlwidgets', 'timelyportfolio/listviewer'))
# library(listviewer)
# jsonedit( candidate_cox )
# # finished (both) processes

# archive .xlsx and transfer .Rda from GCE i-5, i-6 to i-4 (3 in 1) ####
#$ find HNSCC_survivalAnalysis_marginS_*.* -print > survival.txt
#$ tar -czf  HNSCC.survival.marginS.20500.tar.gz -T survival.txt --remove-files

# macOS$ scp tex@35.201.224.219:~/R/HNSCC_Tex_survival/hnscc_github/marginFree/HNSCC_OS_marginFree_pvalueKM_candidate_cox.Rda ./
# $ scp /Users/apple/R/HNSCC_OS_marginFree_pvalueKM_candidate_cox.Rda   tex@35.201.169.0:~/R/HNSCC_Tex_survival/hnscc_github/marginFree/
#
# macOS$ scp tex@35.234.23.175:~/R/HNSCC_Tex_survival/hnscc_github/marginPlus/HNSCC_OS_marginPlus_pvalueKM_candidate_cox.Rda ./
# $ scp /Users/apple/R/HNSCC_OS_marginPlus_pvalueKM_candidate_cox.Rda   tex@35.201.169.0:~/R/HNSCC_Tex_survival/hnscc_github/marginPlus/ 
$ la -al *.Rda




## Post1 cut by P-value and/or Z-score ####
##  table1 (KM): candidate_sample(KM) ##
# #sort by p_value (ascending) and cyl (descending)
# 楊老師: Bonferroni correction (adjustment) sets the significance cut-off at α/n. of KM P-value; Cox P-value <=0.05 as well.
# 王老師: 珍貴的結果 strong evidence
# https://www.stat.berkeley.edu/~mgoldman/Section0402.pdf
# newdata <- mtcars[order(mpg, -cyl),]
# >*** [choice ONE]: _marginFree_ or _marginS_ loading from .Rda
#marginTag <- "_marginS_" #at ./marginS
#marginTag <- "_marginFree_" #at ./marginFree
#marginTag <- "_marginPlus_" #at ./marginPlus 
raw <- 
  readline("_margin(S)_, _marginFree(F)_ or _margin(P)lus_ -- process run: ")
select_margin <- function(x) {
  switch(x,
         s = "_marginS_",
         f = "_marginFree_",
         p = "_marginPlus_",
         S = "_marginS_",
#         F = "_marginFree_",
         P = "_marginPlus_",
         stop("Unknown input")
  )
 }
marginTag <- select_margin(raw)
print(marginTag)
rm(raw) 
# notice: //
path_ZSWIM2 <- file.path(path_cohort, gsub("_", "", marginTag), "") # e.x. marginS/
load(file=paste(path_ZSWIM2, TCGA_cohort, "_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep="")) 
# as candidate_sample with candidate_cox and n_percent_Bonferroni at HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda 


# marginS, marginPlue or marginFree: using a "common" HNSCC_OS_marginS_pvalue_sorted
# run 3 times, 就有三組 marginS, marginPlus and marginFree;
attach(candidate_sample) # n=20500
HNSCC_OS_marginS_pvalue_sorted <- candidate_sample[order(p_value, -number),] # sorting by order(ascending)
detach(candidate_sample)


# ***其實直接 jump to p.adjust() post2 run 454
# x try 人工 Bonferroni correction (adjustment)
# -> sets the significance cut-off at α/n 
# -> p.adjust(p, method = p.adjust.methods="bonferroni", n = length(p)) # bonferroni, or fdr
# n= number of comparisons, must be at least length(p); only set this (to non-default) when you know what you are doing!
# {
attach(HNSCC_OS_marginS_pvalue_sorted)
alpha_HNSCC <- 0.05
Bonferroni_cutoff <- alpha_HNSCC / (LUAD_n * n_percent_Bonferroni)
# => 4.786521e-06 (LUAD run04); => 5.31011e-06 (2019 run02) => 5.323113e-06 [2020/03/12] with 100cut FDR correction
# } 人工 Bonferroni end

# 以下是 plotting (it is not final figure); calculation go 463 ####
# #number/OS_pvalue: from cutoff finder, the number (frequency) of OS P-values, which KM P-value < 0.05, in this gene;
# higher probability to be "found" significantly on a random selection of "cutoff" value in the tranditional manner.
#tiff("Rplot10_Freq_Pvalue.tiff", units="cm", width=5, height=5, res=300)
plot(p_value, number, type="p", ylab="Frequency", xlab="P-value", main="P-value plot of KM survival analyses", log="x", cex=0.3) # log scale x or y
abline(h=100, lty=2, col="blue")
abline(v=Bonferroni_cutoff, lty=2, col="red") # 5.31011e-06
legend("topright", legend=c(paste("Frequency at 100"), paste("Bonferroni ", signif(Bonferroni_cutoff, 2))), lty=2:2, col=c("blue","red"), cex=0.7) # box and font size
#dev.off()
# KM P-value with Bonferroni correction
# ok n=28 in _marginS_ of HNSCC_OS_marginS_pvalueBonferroni_sorted; (1/3)
HNSCC_OS_marginS_pvalueBonferroni_sorted <- HNSCC_OS_marginS_pvalue_sorted[which(p_value<=Bonferroni_cutoff & !is.na(p_value)), 1:3]
# ok n=14 in _marginFree_; (2/3) 
HNSCC_OS_marginFree_pvalueBonferroni_sorted <- HNSCC_OS_marginS_pvalueBonferroni_sorted
# ok n=6 in _marginPlus_; (3/3); return to 258
HNSCC_OS_marginPlus_pvalueBonferroni_sorted <- HNSCC_OS_marginS_pvalueBonferroni_sorted

# afer 3/3 then save 3 SFP
save(HNSCC_OS_marginS_pvalueBonferroni_sorted, 
     HNSCC_OS_marginFree_pvalueBonferroni_sorted, 
     HNSCC_OS_marginPlus_pvalueBonferroni_sorted, 
     file=paste(path_cohort, "/", TCGA_cohort, "_OS", "_marginSFP_pvalueBonferroni_KM_candidate_cox.Rda", sep="")) 
# a.k.a. HNSCC_OS_marginSFP_pvalueBonferroni_KM_candidate_cox.Rda
# go to 413 then try FDR 442
# 
### if no Bonferroni: try z-cut, z_score (Tex發明的)
non_Bonferroni_cutoff <- alpha_HNSCC / 1
plot(p_value, number, type="p", ylab="Frequency", xlab="P-value", main="P-value plot of KM survival analyses", log="x", cex=0.3) # log scale x or y
abline(h=150, lty=2, col="blue")
abline(v=non_Bonferroni_cutoff, lty=2, col="red") # ***as alpha_HNSCC instead of 5.31011e-06
legend("bottomleft", legend=c(paste("Frequency at 150"), paste("P-value at ", non_Bonferroni_cutoff)), lty=2:2, col=c("blue","red"), cex=0.7) # box and font size
# KM P-value <= 0.05 (alpha_HNSCC)
HNSCC_OS_marginS_pvalue005_sorted <- HNSCC_OS_marginS_pvalue_sorted[which(p_value<=non_Bonferroni_cutoff & !is.na(p_value)), 1:3]
detach(HNSCC_OS_marginS_pvalue_sorted)  # keeping drawing
# no correction: n=6624 in _marginS_; n=? in _marginFree_; and n=? in _marginPlus_
# of HNSCC_OS_marginS_pvalue005_sorted


#  plot(OS.km, lty=1, col=c("blue","red"), sub=paste("p-value =", p_OS[j]), main=paste("OS in OSCC(n=", surv_OS$n[1]+surv_OS$n[2],")/", geneName, "cutoff=", i), ylab="Percent Survival", xlab="Days")
#  legend("topright", legend=c(paste("low(",surv_OS$n[1], ")"), paste("high(",surv_OS$n[2], ")")), lty=1:2, col=c("blue","red"))
library(stats)
library(scales)
library(minpack.lm)
# n=6601
attach(HNSCC_OS_marginS_pvalue005_sorted)
# reg <- lm(number ~ p_value, data = HNSCC_OS_marginS_pvalue005_sorted)
# abline(reg, col="blue")
# n=6601 in HNSCC_OS_marginS_pvalue005_sorted
HNSCC_OS_marginS_pvalue005_sorted$z_score <- scale(number, center=T, scale = T) # z_score z_cut were born
# "number" frequency is standardized as z-score (-1.04 to +2.30)
# ("Z" because the normal distribution is also known as the "Z distribution").
# => scale(): scale=TRUE, center=TRUE then scaling is done by dividing the (centered) columns of x by their standard deviations

# cut Z-score at
zcut <- 0.8 # n=1475
boxplot(HNSCC_OS_marginS_pvalue005_sorted$z_score)
abline(h=zcut, lty=2, col="blue") # mean=0, Max.: 2.3034
#DO NOT use Feature Scaling - normalization by (here: scaling to [0, 1] range); Rescaling data to have values between 0 and 1. This is usually called feature scaling. (min-max scaling)
# it will "condensate" the spreading of "number".
#x HNSCC_OS_marginS_pvalue005_sorted$number_01 <- scales:::rescale(z_score, to = c(0, 1))
#x rescaled to range minnew to maxnew (aka. 0 to 1 for binomial glm)
#x number_01 (0~1) is not Z-score, it is Feature Scaling 須要 正名之


#tiff("Rplot10_Zscore_Pvalue.tiff", units="cm", width=5, height=5, res=300)
plot(p_value, z_score, type="p", ylab="Z-score", xlab="P-value", main="P-value plot of KM survival analysis", log="x") # log scale x or y
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
abline(h=zcut, lty=2, col="blue")
abline(v=alpha_HNSCC, lty=2, col="red") # ***as alpha_HNSCC instead of 5.31011e-06

# run a logistic regression model (categorical 0 vs 1)
#g <- glm(number_01 ~ p_value, family=poisson, data=HNSCC_OS_marginS_pvalue005_sorted)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model

# x# run a non-linear regression: non-linear least squares approach (function nls in R) ; A nice feature of non-linear regression in an applied context is that the estimated parameters have a clear interpretation (Vmax in a Michaelis-Menten model is the maximum rate) which would be harder to get using linear models on transformed data for example.
# # the Levenberg-Marquardt algorithm for nonlinear regression
# # "eyeball" plot to set approximate starting values
# # https://stats.stackexchange.com/questions/160552/why-is-nls-giving-me-singular-gradient-matrix-at-initial-parameter-estimates
# # nls "singular gradient matrix at initial parameter estimates" error, using nlsLM instead
# a_start <- 0.9 #param a is the y value when x=0
# b_start <- 2*log(2)/a_start #b is the decay rate
# m <- nlsLM(number_01 ~ p_value, data=HNSCC_OS_marginS_pvalue005_sorted, start=list(a=a_start, b=b_start))
# curve(predict(m, data.frame(p_value = x), type="response", col="green"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
# #lines(p_value, predict(m), lty=2, col="green", lwd=3)

legend("bottomleft", legend=c(paste("Z-score at ", zcut), paste("P-value at ", alpha_HNSCC)), lty=2:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm
#detach(HNSCC_OS_marginS_pvalue005_sorted)

# then...
# a candidate table, with "number of P-value" under cutoff finding (becoming z-score: bigger, much more significant cutoff "sites")
# => a kind of local minimal or global minimal of curve fitting (Levenberg-Marquardt Optimization)?
# subsetting the candidated genes table, according P-value (Bonferroni_cutoff) and Z-score
# ***after Bonferroni correction => n=26 in _marginS_
# ***after Bonferroni correction => n=33 in _marginFree_
#attach(HNSCC_OS_marginS_pvalue005_sorted)
HNSCC_OS_marginS_pvalue005_zcut <- HNSCC_OS_marginS_pvalue005_sorted[which(p_value <= alpha_HNSCC & z_score >= zcut), -1]
# discard "number"
HNSCC_OS_marginS_pvalue005_zcut$p_value <- signif(HNSCC_OS_marginS_pvalue005_zcut$p_value, 3)
HNSCC_OS_marginS_pvalue005_zcut$z_score <- signif(HNSCC_OS_marginS_pvalue005_zcut$z_score, 3)
detach(HNSCC_OS_marginS_pvalue005_sorted)
# > colnames(HNSCC_OS_marginS_pvalue005_zcut)
# x[1] "gene_id"   "p_value"   "z_score"   "number_01"
# [1] "gene_id" "p_value" "z_score"
# _marginS_ (ok) as marginTag

# *** there is 3 collections:
# save n=28, cut by P-value <= Bonferroni_cutoff (Bonferroni correction)
save(HNSCC_OS_marginS_pvalueBonferroni_sorted, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalueBonferroni_sorted.Rda", sep="")))
#
# #R4> HNSCC_OS_marginS_pvalueBonferroni_sorted$gene_id
# [1] "ZNF557"    "DKK1"      "ZNF266"    "CAMK2N1"  
# [5] "IL19"      "MYO1H"     "STC2"      "PGK1"     
# [9] "SURF4"     "FCGBP"     "LOC148709" "USP10"    
# [13] "EVPLL"     "PNMA5"     "NDFIP1"    "FOXA2"    
# [17] "GRAP"      "KIAA1683"  "ZNF846"    "ZNF20"    
# [21] "CELSR3"    "SLC26A9"   "FAM3D"     "GPR15"    
# [25] "KLRA1"     "NPB"       "TCP11"     "STIP1" 

# venn of 3: marginS, marginFree and marginPlus:
venn_marginSFP <- list(HNSCC_OS_marginFree_pvalueBonferroni_sorted$gene_id, HNSCC_OS_marginS_pvalueBonferroni_sorted$gene_id, 
                       HNSCC_OS_marginPlus_pvalueBonferroni_sorted$gene_id)
names_marginSFP <- c(paste("margin[-]"), paste("margin[+/-]"), paste("margin[+]"))
library(venn)
# https://cran.r-project.org/web/packages/venn/venn.pdf
#tiff("venn_marginSFP.tiff", units="cm", width=5, height=5, res=300) 
# # saving as .tiff (by tiff())
venn(venn_marginSFP, snames=names_marginSFP,
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "red, deeppink, pink", opacity = 0.6, size = 15, cexil = 2.0, cexsn = 0.8, borders = FALSE)
#      predefined colors if "style"
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# meta-language 1 0 or -
title <- c(paste(TCGA_cohort, "survival analyses: candidate genes"), paste("(KM P-Value <=", signif(Bonferroni_cutoff, 2), ")")) #, collapse = "\n")
#coords <- unlist(getCentroid(getZones(venn_HR2p5, snames="uni_CoxHR>=2p5, multi_CoxHR>=2p5")))
# coords[1], coords[2], 
text(500,950, labels = title[1], cex = 1.30) # (0,0) on bottom_left corner
text(500,100, labels = title[2], cex = 1.20) 
# n=47
#dev.off() # saving instead of showing


# retrieving the intersection genes by gplots::venn
library(gplots)
tmp_cand <- venn(venn_marginSFP, names=names_marginSFP, show.plot=F) #library(gplots); the group count matrix alone
isect_marginSFP <- attr(tmp_cand, "intersections")
#isect_marginSFP
detach(package:gplots)


## > Post2 KM adding candidate_cox ####
## # (ok)重新改寫 [Bonferroni] and [FDR] P-value correction
# with FDR or Bonferroni correction ##
# *** best [2019/11/07] pValue_adj <- p.adjust(survOutput_km$pValue, method="fdr") #, n = nrow(survOutput_km))
# the possibility of type I error: alpha (single test)
# multiple testing problem
# FDR = E[F/S] <= q-value; http://www.omicshare.com/forum/thread-173-1-1.html
# >*** [choice ONE]: _marginFree_ or _marginS_ loading from .Rda
#marginTag <- "_marginS_" #at ./marginS
#marginTag <- "_marginFree_" #at ./marginFree
#marginTag <- "_marginPlus_" #at ./marginPlus 
raw <- 
  readline("_margin(S)_, _marginFree(F)_ or _margin(P)lus_ -- process run: ")
select_margin <- function(x) {
  switch(x,
         s = "_marginS_",
         f = "_marginFree_",
         p = "_marginPlus_",
         S = "_marginS_",
         #         F = "_marginFree_",
         P = "_marginPlus_",
         stop("Unknown input")
  )
}
marginTag <- select_margin(raw)
print(marginTag)
rm(raw) 
# notice: //
path_ZSWIM2 <- file.path(path_cohort, gsub("_", "", marginTag), "") # e.x. marginS/
load(file=paste(path_ZSWIM2, TCGA_cohort, "_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep="")) 
# as candidate_sample with candidate_cox and n_percent_Bonferroni at HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda 
# after 100cut FDR correction and selection


# marginS, marginPlue or marginFree: using a "common" HNSCC_OS_marginS_pvalue_sorted
# run 3 times, 就有三組 marginS, marginPlus and marginFree;
attach(candidate_sample) # n=20500
HNSCC_OS_marginS_pvalue_sorted <- candidate_sample[order(p_value, -number),] # sorting by order(ascending)
detach(candidate_sample)


#install.packages("BiocManager")
library("BiocManager")
# #source("https://bioconductor.org/biocLite.R")
# #biocLite("GDCRNATools") 
# #biocLite("BiocParallel") 
# install.packages("remotes")
# library(remotes)
# devtools::install_github("jeroen/jsonlite")
# 
# # https://github.com/Jialab-UCR/Jialab-UCR.github.io/blob/master/GDCRNATools.workflow.R
# remotes::install_github("Jialab-UCR/GDCRNATools")
# library(GDCRNATools) # Pathview, citation("pathview") within R
# #
# a simple way
# na.omit -> only n=6624 genes have P-value available; removal of NA
HNSCC_OS_marginS_pvalue_sorted_noNA <- HNSCC_OS_marginS_pvalue_sorted[complete.cases(HNSCC_OS_marginS_pvalue_sorted), ]
#x keeping z_score
# HNSCC_OS_marginS_pvalue_sorted_noNA$z_score <- scale(HNSCC_OS_marginS_pvalue_sorted_noNA$number, center=T, scale = T)

# 1) try bonferroni by p.adjust()
#p_value_adj_bonferroni <- p.adjust(HNSCC_OS_marginS_pvalue_sorted_noNA$p_value, method="bonferroni", n = nrow(HNSCC_OS_marginS_pvalue_sorted_noNA))
p_value_adj_bonferroni <- p.adjust(HNSCC_OS_marginS_pvalue_sorted_noNA$p_value, method="bonferroni") # n=default 不必寫

# 2) estimate FDR from p-values by p.adjust()
# https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
# by p.adjust(); n=6615 of FDR<0.05; n=9 FDR=0.05
p_value_adj_FDR <- p.adjust(HNSCC_OS_marginS_pvalue_sorted_noNA$p_value, method="fdr", n = nrow(HNSCC_OS_marginS_pvalue_sorted_noNA))
# trying FDR<0.02, n=8
# p_value_adj[p_value_adj<=0.02]; n=8
# p_value_adj[p_value_adj<=0.01]; n=0

HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS <- cbind(HNSCC_OS_marginS_pvalue_sorted_noNA, p_value_adj_bonferroni, p_value_adj_FDR)
# HNSCC_OS_marginS_pvalue_sorted_noNA <- cbind(HNSCC_OS_marginS_pvalue_sorted_noNA, p_value_adj)
#HNSCC_OS_marginS_pvalue_sorted_noNA_FDR <- HNSCC_OS_marginS_pvalue_sorted_noNA
save(HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue_sorted_noNA_p_adjustS.Rda", sep="")))
# [2020/02/16] p_value_adj is FDR and bonferroni;
# "HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS.Rda" saved at "/home/tex/R/HNSCC_Tex_survival/hnscc_github/marginS/"
# n=6624 => n=6429


# x[2019/09/05]
#HNSCC_OS_marginS_pvalue005_sorted
# save n=6624, cut by P-value <= alpha_HNSCC
#x save(HNSCC_OS_marginS_pvalue005_sorted, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue005_sorted.Rda", sep="")))
# save n=1437, cut by P-value <= alpha_HNSCC AND Z-score >= zcut(e.q. 0.8)
#x save(HNSCC_OS_marginS_pvalue005_zcut, Bonferroni_cutoff, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue005_zcut.Rda", sep=""))) 
# *** rename HNSCC_OS_marginS_pvalue1e_6_zscore0_6 => HNSCC_OS_marginS_pvalue005_zcut

# 
#x or
#x x) by fdrtool() with chart; n=6615 FDR<0.05
## https://www.rdocumentation.org/packages/fdrtool/versions/1.2.15/topics/fdrtool
#xinstall.packages("fdrtool")
#xHNSCC_fdr = fdrtool(HNSCC_OS_marginS_pvalue_sorted_noNA$p_value, statistic="pvalue", cutoff.method=c("fndr"))
#cutoff.method=c("fndr", "pct0", "locfdr")
# Step 1... determine cutoff point
# Step 2... estimate parameters of null distribution and eta0
# Step 3... compute p-values and estimate empirical PDF/CDF
# Step 4... compute q-values and local fdr
# Step 5... prepare for plotting
#HNSCC_fdr$qval # estimated Fdr values => all are zero
#HNSCC_fdr$lfdr # estimated local fdr => all are zero
#---


 
#library("BiocParallel")
#register(MulticoreParam(2)) # 2 cores





##xx it has [part I][part II][part III]
# "sig" marking for significant P-value (<=0.05)
# [uni_cox_pvalue, uni_HR, uni_sig]
# [multi_cox_pvalue, multi_HR, multi_sig]
# [exp_pvalue]
# to generate HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR.Rda; n=6601 or 1476
# from 
# x# 1) HNSCC_OS_marginS_pvalue005_sorted # n=6601; P-value < 0.05 (without correction)
# load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, 
#                                        "pvalue005_sorted.Rda", sep="")))
# x# 2) HNSCC_OS_marginS_pvalue005_zcut, n=1476;
# load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, 
#                                        "pvalue005_zcut.Rda", sep="")))
# x# or 3) # *** try using FDR: HNSCC_OS_marginS_pvalue_sorted_noNA_FDR.Rda; FDR
# load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue_sorted_noNA_FDR.Rda", sep="")))


load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue_sorted_noNA_p_adjustS.Rda", sep="")))
# as HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS
#
load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep="")))
# as candidate_sample, candidate_cox, n_percent_Bonferroni (after 100cut FDR)

# [Bonferroni][FDR] new variable: KM + Cox, n=6429; P-value < 0.05 (Bonferroni correction)
# # KM P-value <= 0.05 (alpha_HNSCC): HNSCC_OS_marginS_pvalue005_sorted
HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR <- HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS
#dataframe[,"newName"] <- NA # add more named columns (NA)
HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[,c("uni_HR", "uni_P_value", "uni_sig","multi_HR", "multi_P_value", "multi_sig")] <- NA
#HNSCC_OS_marginS_pvalue005_sorted[,c(colnames(candidate_cox[[1]][8, c(2:4, 6:8)]))] <- NA
#colnames(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR) <- c(...."uni_HR", "uni_P_value", "uni_sig",
#                                                "multi_HR", "multi_P_value", "multi_sig")

#attach(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
#ipp <- 0
# *** row and col position 非常重要
for (ip in 1:nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)) {
  pos_gene <- which(whole_genome==HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id[ip]) # HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id
  # by geneid ?
  if (length(candidate_cox[[pos_gene]])==11) {
    # reference: candidate_cox_ip, which listing as candidate_cox
    HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, 6:11]  <- candidate_cox[[pos_gene]][9, c(2:4, 6:8)] 
                                                                # taking 9 RNAseg: sig marked and P-values, HRs
                                                                #*** 8 is tobacco; 
        # we don't have number0_1 any more:  colnames "number"  "gene_id" "p_value" "z_score"
    #   ipp <- ipp+1
    print(paste(ip, " out of ", nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), " (", whole_genome[pos_gene], ")", sep=""))
    print(paste(c("..adding...", HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, c(2, 6, 9)]), sep=";"))
  }
  # on : [2018-05-18 07:26:42] [error] handle_read_frame error: websocketpp.transport:7 (End of File)
  #"exp_pvalue" is a correlation of gene expression vs features (TNM....): 
  # <- candidate_cox[which(gene_id[ip]==whole_genome)][1:7, c(11)] # correlation
  # colnames <-  c("KM_Features", "KM_P_value",  "KM_sig")
  
}
# HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS with [Bonferroni][FDR] and variable: KM + Cox
save(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue005KM_sorted_pvalueCox_HR.Rda", sep="")))
# as HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR.Rda
# 2020/02/24 it is tobacco :-); 2020/03/08 ok RNA-seq; go "Pickup all..."


# # xx
# #x [part I] new variable: KM + Cox, n=6601; P-value < 0.05 (without correction)
# #xHNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR <- HNSCC_OS_marginS_pvalue005_sorted
# #dataframe[,"newName"] <- NA # add more named columns (NA)
# #xHNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[,c("uni_HR", "uni_P_value", "uni_sig","multi_HR", "multi_P_value", "multi_sig")] <- NA
# #HNSCC_OS_marginS_pvalue005_sorted[,c(colnames(candidate_cox[[1]][8, c(2:4, 6:8)]))] <- NA
# #colnames(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR) <- c(...."uni_HR", "uni_P_value", "uni_sig",
# #                                                "multi_HR", "multi_P_value", "multi_sig")
# 
# #attach(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# #ipp <- 0
# for (ip in 1:nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)) {
#   pos_gene <- which(whole_genome==HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id[ip]) # HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id
#   # by geneid ?
#   if (length(candidate_cox[[pos_gene]])==11) {
#     # reference: candidate_cox_ip, which listing as candidate_cox
#     HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, 5:10]  <- candidate_cox[[pos_gene]][8, c(2:4, 6:8)] # taking RNAseg: sig marked and P-values, HRs
#     # we don't have number0_1 any more:  colnames "number"  "gene_id" "p_value" "z_score"
#     #   ipp <- ipp+1
#     print(paste(ip, " out of ", nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), " (", whole_genome[pos_gene], ")", sep=""))
#     print(paste(c("..adding...", HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[ip, c(2, 5)]), sep=";"))
#   }
#   # on : [2018-05-18 07:26:42] [error] handle_read_frame error: websocketpp.transport:7 (End of File)
#   #"exp_pvalue" is a correlation of gene expression vs features (TNM....): 
#   # <- candidate_cox[which(gene_id[ip]==whole_genome)][1:7, c(11)] # correlation
#   # colnames <-  c("KM_Features", "KM_P_value",  "KM_sig")
#   
# }
# #detach(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# # c( "uni_Features", "uni_HR",       "uni_P_value",  "uni_sig",      "multi_Features", "multi_HR",
# #"multi_P_value",  "multi_sig",      "KM_Features", "KM_P_value",  "KM_sig")
# # add colname as HRs, "uni_cox"/"multi_cox", sig ... for HRs, P-values, sig of RNAseq(z-score)
# 
# 
# # [part II] new variable: KM + Cox, n=1475; under zcut
# HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR <- HNSCC_OS_marginS_pvalue005_zcut
# # add more columns
# HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR[,c("uni_HR", "uni_P_value", "uni_sig","multi_HR", "multi_P_value", "multi_sig")] <- NA
# for (ip in 1:nrow(HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR)) {
#   pos_gene <- which(whole_genome==HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR$gene_id[ip]) # HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR$gene_id
#   if (length(candidate_cox[[pos_gene]])==11) {
#     # reference: candidate_cox_ip, which listing as candidate_cox
#     HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR[ip, 4:9]  <- candidate_cox[[pos_gene]][8, c(2:4, 6:8)] # taking RNAseg: sig marked and P-values, HRs
#     # we don't have number0_1 any more:   colnames "gene_id" "p_value" "z_score"
#     #   ipp <- ipp+1
#     print(paste(ip, " out of ", nrow(HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR), " (", whole_genome[pos_gene], ")", sep=""))
#     print(paste(c("..adding...", HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR[ip, c(1, 4)]), sep=";"))
#   }
# }
# ##
# # save table1 + table2 (ok) [2019/07/01], n= 6601 or 1475
# save(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue005KM_sorted_pvalueCox_HR.Rda", sep="")))
# save(HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue005_zcutKM_sorted_pvalueCox_HR.Rda", sep="")))
# #save(list(HNSCC_OS_marginS_pvalue005_zcut, ???))
# x[i], or might be x[i:j]
# x[i, j]
# x[[i]]; x[[expr]]; it can NOT be x[[i:j]]
# x[[i, j]]
# x$a
# x$"a"

# # [part III] under FDR q-value
# # # HNSCC_OS_marginS_pvalue_sorted_noNA_FDR.Rda
# HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR <- HNSCC_OS_marginS_pvalue_sorted_noNA_FDR
# # add more columns
# HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR[,c("uni_HR", "uni_P_value", "uni_sig","multi_HR", "multi_P_value", "multi_sig")] <- NA
# for (ip in 1:nrow(HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR)) {
#   pos_gene <- which(whole_genome==HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR$gene_id[ip]) # HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR$gene_id
#   if (length(candidate_cox[[pos_gene]])==11) {
#     # reference: candidate_cox_ip, which listing as candidate_cox
#     HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR[ip, 5:10]  <- candidate_cox[[pos_gene]][8, c(2:4, 6:8)] # taking RNAseg: sig marked and P-values, HRs
#     # we don't have number0_1 any more:   colnames "gene_id" "p_value" "z_score"
#     #   ipp <- ipp+1
#     print(paste(ip, " out of ", nrow(HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR), " (", whole_genome[pos_gene], ")", sep=""))
#     print(paste(c("..adding...", HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR[ip, c(1, 5)]), sep=";"))
#   }
# }
# #
# save(HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR.Rda", sep="")))
# ## =6624



#> Pickup all significant genes list -> HNSCC_OS_marginS_THREE_pvalue005 ####
# x FDR < 0.05; HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR, n=6624
# x Z-score > 0.8: HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR, n=1475

### => HNSCC_OS_marginS_pvalue_sorted_noNA_p_adjustS with [Bonferroni][FDR] and variable: KM + Cox HR
# as HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR
# _marginS_ on [2019/07/02]
load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue005KM_sorted_pvalueCox_HR.Rda", sep="")))

# # > colnames(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR)
# [Bonferroni] and [FDR] P-value corrected (no z_cut)
# [1] "number"                
# [2] "gene_id"               
# [3] "p_value"               
# [4] "p_value_adj_bonferroni"
# [5] "p_value_adj_FDR"       
# [6] "uni_HR"                
# [7] "uni_P_value"           
# [8] "uni_sig"               
# [9] "multi_HR"              
# [10] "multi_P_value"         
# [11] "multi_sig"              
# R4> 
#xHNSCC_OS_marginS_THREE_pvalue005_FDR <- subset(HNSCC_OS_marginS_pvalue_sorted_noNA_FDRKM_sorted_pvalueCox_HR, 
#                                                select=c(gene_id, p_value, p_value_adj, uni_HR, uni_P_value, multi_HR, multi_P_value))
                                                # p_value_adj is FDR
# n=6429
HNSCC_OS_marginS_THREE_pvalue005 <- subset(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, 
                                           select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
save(HNSCC_OS_marginS_THREE_pvalue005, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005.Rda", sep=""))) 

#xHNSCC_OS_marginS_THREE_pvalue005_1475 <- subset(HNSCC_OS_marginS_pvalue005_zcutKM_sorted_pvalueCox_HR, 
#                                           select=c(gene_id, p_value, z_score, uni_HR, uni_P_value, multi_HR, multi_P_value))
#... <- subset(..., (uni_P_value <= 0.05) & (multi_P_value <= 0.05), 
#save(HNSCC_OS_marginS_THREE_pvalue005_FDR, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005_FDR.Rda", sep=""))) 

#save(HNSCC_OS_marginS_THREE_pvalue005_6601, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005_6601.Rda", sep=""))) 
#save(HNSCC_OS_marginS_THREE_pvalue005_1475, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005_1475.Rda", sep=""))) 
# as HNSCC_OS_marginS_THREE_pvalue005 (6601 or 1475)


# >pubmed.mineR ####
# # no zcut: HNSCC_OS_marginS_THREE_pvalue005_6601 <- HNSCC_OS_marginS_THREE_pvalue005, n=6601
# load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005_6601.Rda", sep="")))
# # zcut, Z-score > 0.8: HNSCC_OS_marginS_THREE_pvalue005_1475 <- HNSCC_OS_marginS_pvalue005_zcut, n=1475
# load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005_1475.Rda", sep="")))
# # if FDR < 0.05
# load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005_FDR.Rda", sep="")))


# start to working...# *** 
load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005.Rda", sep="")))
# as HNSCC_OS_marginS_THREE_pvalue005
# [GDC data portal]
# https://www.rdocumentation.org/packages/pubmed.mineR/versions/1.0.1
# https://docs.gdc.cancer.gov/Data_Portal/Users_Guide/Exploration/
# Excluding COSMIC Cancer Gene Census (CGC, such as TRIM24, ...) from GDC data portal: d/l Cancer Gene Census (575 genes)
# [2019/07/17] downloaded, n=575; https://cancer.sanger.ac.uk/census/? 再也找不到了
# [2020/02/27] Census_allThu Feb 27 09_29_04 2020.tsv: total n=723; n=7+12=19 in tumour types: head and neck cancer, HNSCC, oral SCC, oral squamous cell carcinoma]
# load genes.2019-07-17.json
#install.packages("jsonlite")
library(jsonlite)
#x cgc_gene575 <- fromJSON("~/R/genes.2019-07-17.json", simplifyDataFrame=TRUE)
cgc_gene723 <- read.table(file = "~/R/Census_allThu Feb 27 09_29_04 2020.tsv", sep = '\t', header = TRUE)
# subsetting by omitting row
# xcc_Index2019 <- which(HNSCC_OS_marginS_THREE_pvalue005$gene_id %in% cgc_gene575$symbol)
cgc_Index2020 <- which(HNSCC_OS_marginS_THREE_pvalue005$gene_id %in% cgc_gene723$Gene.Symbol)
#venn(list(cgc_gene575$symbol[cc_Index2019], cgc_gene723$Gene.Symbol[cc_Index2020]), snames=c("CGC_2019", "CGC_2020"), ilabels = T, counts = T, ellipse = FALSE, zcolor = "red, deeppink, pink", opacity = 0.6, size = 15, cexil = 2.0, cexsn = 0.8, borders = FALSE)
View(HNSCC_OS_marginS_THREE_pvalue005$gene_id[cgc_Index2020]) # n=213 or 52
HNSCC_OS_marginS_THREE_pvalue005_noCancerGene <- HNSCC_OS_marginS_THREE_pvalue005[-cgc_Index2020, ] # removal
# now n=6368 [2020/03/08]

# There is updated HNSCC related gene from publication ####
# *** go to embase

# x [PubMed] search keywoard:"((((marker[Abstract] OR biomarker[Abstract])) AND (HNSCC OR HNSC OR "oral"))) AND (prognostic[Title] OR prognosis[Title])"
#install.packages("pubmed.mineR") #  two formats (text and XML) from PubMed
# scp /Users/apple/Downloads/pmc_result_HNSCC_genes1486_MEDLINE.txt tex@35.201.169.0:~/R/
# *** export MEDLINE format (with abstract)
library(pubmed.mineR) # http://ramuigib.blogspot.com
#install.packages(c("rentrez", "XML"))
library(rentrez)
library(XML)
# gene_atomization(), official_fn(gene_atomization, abs, filename, terms), 
#abs <- xmlreadabs("~/R/pmc_result_HNSCC_genes1486.xml") #an S4 object of class Abstracts containing journals, abstracts and PMID
#abs <- read.table("~/R/pmc_result_HNSCC_genes1486.xml", sep="|") #an S4 object of class Abstracts containing journals, abstracts and PMID
abs <- readabs("~/R/pmc_result_HNSCC_genes1486_MEDLINE.txt") # Abstracts <- .txt; the slot operator @ 
# # 5.5Mb; objects of class "Abstracts" have 3 slots, Journal, Abstract, and PMID
# readabs(), xmlgene_atomizations(m)
# x pub_hnscc_genes <- xmlgene_atomizations(abs)
#pub_hnscc_genes <- gene_atomization(abs) # n=615; however, "FAU" "T" "GC" "HR" are not gene symbol
## (there is misleading data)
#pub_Index <- which(pub_hnscc_genes[ , 1] %in% whole_genome) # curation1; n=611
#pub_hnscc_genes <- pub_hnscc_genes[pub_Index, ] # n=611
# https://www.rdocumentation.org/packages/pubmed.mineR/versions/1.0.1
# https://omictools.com/pubtator-tool PubTator
#x curation2 <- mapply(Give_Sentences, pub_hnscc_genes[ , 2], abs@Abstract) # or searchabsT
#x curation2 <- mapply(searchabsT, abs, include=pub_hnscc_genes[ , 2]) # or searchabsT
#  getabs(abs, pub_hnscc_genes[1, 2], FALSE)

pubmed_hnscc_PMID <- entrez_search(db="pubmed", 
                                   term="((((((((marker[Abstract] OR biomarker[Abstract])) AND (prognostic[Title] OR prognosis[Title]))))) AND (((Upper Aerodigestive) OR HNSCC OR (oral cancer)))))) NOT ((salivary OR bladder OR breast OR bone OR kidney OR leukoplakia OR dysplasia))", retmax=2000, use_history=TRUE)
# n=1325
# x=> refinement with ((D006258[MeSH Major Topic]) OR D006258[MeSH Terms]) <= Head and Neck Neoplasms
#x pubmed_hnscc_PMID <- entrez_search(db="pubmed", term="((((marker[Abstract] OR biomarker[Abstract])) AND ((D006258[MeSH Major Topic]) OR D006258[MeSH Terms]) AND (prognostic[Title] OR prognosis[Title])", retmax=5000)
# MESH is not universal available :-)

# paper_links <- entrez_link(dbfrom="pubmed", id=25500142, cmd="llinks")
# paper_links$linkouts; linkout_urls(paper_links)
# entrez_fetch(); entrez_summary()
# error: download abstracts (HTTP failure 414, the request is too large. For large requests, try using "web history" as described in the rentrez tutorial)
# "web history" { , use_history=TRUE => web_history object
#upload_ncbi <- entrez_post(db="pubmed", id = pubmed_hnscc_PMID$ids)
#upload_ncbi

pubmed_hnscc_abs <- entrez_fetch(db="pubmed", web_history=pubmed_hnscc_PMID$web_history,
                       rettype="xml", retmax=2000, parsed = T)
#R4> str(pubmed_hnscc_abs)
#Classes 'XMLInternalDocument', 'XMLAbstractDocument' <externalptr> 
# *.txt generated
# /Users/apple/Downloads/dataout.txt => gene without abstract available
# /Users/apple/Downloads/newabs.txt => abstracts
# /Users/apple/Downloads/result.txt => ditto with pubmed_hnscc_PMID_list **
# /Users/apple/Downloads/table.txt => genes list
# }

# x pubmed_hnscc_abs <- entrez_fetch(db = "pubmed", id = pubmed_hnscc_PMID$ids,
#                rettype = "xml", parsed = T)

# >text annotation uisng online PubTator (rentrez::pubtator_function)
#pub_hnscc_PMID <- lapply(abs@PMID, pubtator_function) # Gene, Chemical, Disease and PMID
# ** ditto with result.txt from entrez_fetch
pubmed_hnscc_PMID_list <- lapply(pubmed_hnscc_PMID$ids, pubtator_function) # Gene, Chemical, Disease and PMID
# => n=1325; list of Genes, Diseases, ... and PMID
# e.x. 
# R4> pubmed_hnscc_PMID_list[[555]]$Diseases
#$Diseases
#[1] "cancer>MESH:D009369"                  
#[2] "hematopoiesis>MESH:C536227"           
#[3] "inflammation>MESH:D007249"            
#[4] "tumor>MESH:D009369"                   
#[5] "squamous cell carcinoma>MESH:D002294" 
#[6] "squamous cell carcinomas>MESH:D002294"
#[7] "tumors>MESH:D009369" 
# e.x. pubmed_hnscc_PMID[[555]]$Genes is TATA-binding protein>6908 => Entrez Gene: 6908
#pub_hnscc_table <- pubtator_result_list_to_table(pub_hnscc_PMID) # to exclude HNSCC "disease"
# check genes with HNSCC OR HNSC OR "oral cancer" related
pubmed_hnscc_df <- as.data.frame(matrix(NA, nrow=length(pubmed_hnscc_PMID_list), ncol=3))
colnames(pubmed_hnscc_df) <- c("PMID", "Diseases", "Genes")
for (pubtator_i in 1: length(pubmed_hnscc_PMID_list)) {
  if (pubmed_hnscc_PMID_list[[pubtator_i]] == " No Data ") {next}
  else {pubmed_hnscc_df[pubtator_i, 1] <- pubmed_hnscc_PMID_list[[pubtator_i]]$PMID[1]}

  tryCatch({
    pubmed_hnscc_df[pubtator_i, 2] <- pubmed_hnscc_PMID_list[[pubtator_i]]$Diseases[1]
    pubmed_hnscc_df[pubtator_i, 3] <- pubmed_hnscc_PMID_list[[pubtator_i]]$Genes[1]
    }, error = function(err) {
      # error handler picks up where error was generated
      #print(paste("MY_ERROR:  ", err))
      return()
    }) # END tryCatch
}

pubmed_hnscc_abstractFull <- read.csv(file.path(path_ZSWIM2, "newabs.txt"), header=T) 
# /Users/apple/Downloads/newabs.txt => abstracts
# 29453876: ameloblastoma, Genes of FGF2, Bcl-2
# excluding (salivary OR bladder OR breast OR bone OR kidney OR leukoplakia OR dysplasia) 

#[embase]
# >embase search then exported .csv ####
# https://dev.elsevier.com/documentation/EmbaseAPI.wadl
# https://www.rdocumentation.org/packages/rscopus/versions/0.6.3/topics/scopus_search
#https://cran.r-project.org/web/packages/rscopus/rscopus.pdf
#install.packages("rscopus")
library(rscopus)
library(pubmed.mineR) # http://ramuigib.blogspot.com
#library(RISmed) # https://amunategui.github.io/pubmed-query/; 2017 no more read.ris()
# https://cran.r-project.org/web/packages/RefManageR/RefManageR.pdf

#embase_retrieval(id, identifier = c("lui", "pii", "doi", "embase",
#                                    "pubmed_id", "medline"), http_end = NULL, ...)
# Type of search. See https://dev.elsevier.com/api_docs.html
# [2019/07/20] my elsevier api_key = "9925b63f77af6f11b74b38a3d6a80f41" from web site
# [2020/02/29] manual Quick search in embase: 
# (prognosis:ti OR prognostic:ti OR 'tumor marker') AND cancer:ti AND 'head and neck squamous cell carcinoma' OR ('SLC2A4' AND hnscc) OR ('thymosin beta 4' AND hnscc) OR (hnscc AND 'hsiao m.':au AND 'chang w.-m':au)
# download records-503.csv -> records-432.csv
res_df <- read.csv(file=file.choose(), header=T) # with abstract, PUI as embase
em_hnscc_PMID_list <- lapply(res_df$Medline.PMID, pubtator_function) # list of Gene, Chemical, Disease and PMID
# => n=371 => 134; 503
# e.x. The neuropeptide genes SST, TAC1, HCRT, NPY, and GAL are powerful epigenetic biomarkers in head and neck cancer: a site-specific analysis.
#(PMID:29682090 PMCID:PMC5896056) => The somatostatin (SST)
# parsing and pickup 1st: "PMID", "Diseases", "Genes"
pubmed_hnscc_df <- as.data.frame(matrix(NA, nrow=length(em_hnscc_PMID_list), ncol=3))
colnames(pubmed_hnscc_df) <- c("PMID", "Diseases", "Genes")
for (pubtator_i in 1: length(em_hnscc_PMID_list)) {
  print(pubtator_i)
  if (em_hnscc_PMID_list[[pubtator_i]] == " No Data ") {next}
  else {pubmed_hnscc_df[pubtator_i, 1] <- em_hnscc_PMID_list[[pubtator_i]]$PMID[1]}
  
  tryCatch({
    pubmed_hnscc_df[pubtator_i, 2] <- em_hnscc_PMID_list[[pubtator_i]]$Diseases[1]
    pubmed_hnscc_df[pubtator_i, 3] <- em_hnscc_PMID_list[[pubtator_i]]$Genes[1]
  }, error = function(err) {
    # error handler picks up where error was generated
    #print(paste("MY_ERROR:  ", err))
    return()
  }) # END tryCatch
}
# [HNSCC gene signature]
#Cancer GeneticsWeb: check it out
#  http://www.cancerindex.org/geneweb/HRPT2.htm

# find published HNSCC genes from Embase
# pubmed_hnscc_df data cleaning: sorting by "Disease", manual curation
pubmed_hnscc_df_sorted <- pubmed_hnscc_df[order(pubmed_hnscc_df$Diseases, na.last=NA), ] # NA removal, n=259
pubmed_hnscc_df_sorted <- pubmed_hnscc_df_sorted[complete.cases(pubmed_hnscc_df_sorted$Genes), ] # n=169 => 64
# n=203; including DDX58; GLUT4; HNSCC; Metastasis; TRIM24, PMID: 28061796 
# manually add GLUT4>6517 and DDX58>23586 from 28061796
pubmed_hnscc_df_sorted <- rbind(pubmed_hnscc_df_sorted, pubmed_hnscc_df_sorted[pubmed_hnscc_df_sorted$PMID==28061796, ]) # n=204
pubmed_hnscc_df_sorted <- rbind(pubmed_hnscc_df_sorted, pubmed_hnscc_df_sorted[pubmed_hnscc_df_sorted$PMID==28061796, ]) # n=206
pubmed_hnscc_df_sorted <- pubmed_hnscc_df_sorted[-(nrow(pubmed_hnscc_df_sorted)), ] # n=205
pubmed_hnscc_df_sorted[(nrow(pubmed_hnscc_df_sorted)), 3] <- "GLUT4>6517"
pubmed_hnscc_df_sorted[(nrow(pubmed_hnscc_df_sorted)-1), 3] <- "DDX58>23586"

# loading all genes mentioned in each article
pubmed_hnscc_genes <- as.data.frame(matrix(NA, nrow=nrow(pubmed_hnscc_df_sorted)*4, ncol=5))
colnames(pubmed_hnscc_genes) <- c("No", "PMID", "Diseases", "GeneSymbol", "EntrezGene")
No_i <- 1 # a serial number

for (curation_i in 1: nrow(pubmed_hnscc_df_sorted)) {
  print(curation_i)
  No_j <- 1 # index of genes in each article
  # em_hnscc_PMID_list is 1:371 from res_df$Medline.PMID; there is duplicated article (PMID:24565202)
  res_index <- which(res_df$Medline.PMID==pubmed_hnscc_df_sorted$PMID[curation_i]) # always choice res_index[1], first PMID
  for (No_j in 1:length(em_hnscc_PMID_list[[res_index[1]]]$Genes)) {
    pubmed_hnscc_genes[No_i + No_j - 1, ] <- c((No_i + No_j - 1), pubmed_hnscc_df_sorted$PMID[curation_i], 
                                pubmed_hnscc_df_sorted$Diseases[curation_i], 
                                unlist(strsplit(em_hnscc_PMID_list[[res_index[1]]]$Genes[No_j], ">")))
  }
  No_i <- No_i + No_j
# EntrezGene 4313 (MMP2), 4318 (MMP9); 83639 (TEX101) at PMID 23481573
# Cancer Biomark. 2012-2013;12(3):141-8. doi: 10.3233/CBM-130302.
#  [Overexpression of TEX101, a potential novel cancer marker, in head and neck squamous cell carcinoma.]
} # from 1 to 976 entries in [pubmed_hnscc_genes]
# there is GLUT4, RUNX2 and TMSB4X => repeat searching including 手工加入(ok)
#  ('glucose transporter 4' OR 'thymosin beta 4') AND ('hsiao m.' or 'Chang W.-M'):au (yes)
#  'mRNA expression level'
#  if 'Chang W.-M' AND 'hsiao m.':au => 
# resulting pubmed_hnscc_genes, 256 entries
# \Sci Rep. 2017[Parathyroid Hormone-Like Hormone is a
# Poor Prognosis Marker of Head and Neck Cancer and Promotes Cell Growth via RUNX2 Regulation.]
# n=820
pubmed_hnscc_genes <- pubmed_hnscc_genes[complete.cases(pubmed_hnscc_genes$GeneSymbol), ]
# n=660
# removal non-HNSCC article: PMID 23582651, 26293675, 25487446, 26745068, 16092551, 30125310; and TNM SCC (it is not gene name)
pubmed_hnscc_genes <- pubmed_hnscc_genes[!(pubmed_hnscc_genes$PMID %in% c(23582651, 26293675, 25487446, 26745068, 16092551, 30125310, 23672832, 11329776, 25639759, 30636185, 20629076, 25608735, 23818362, 24405882, 24190760)), ]
# removal one item HNSCC-ALDH1; and rest of TNM
pubmed_hnscc_genes <- pubmed_hnscc_genes[!(pubmed_hnscc_genes$GeneSymbol %in% c("HNSCC-ALDH1", "TNM")), ]
# n=635
save(pubmed_hnscc_genes, file=file.path(path_cohort, "pubmed_hnscc_genes.Rda"))

# then checking duplicated by EntrezGene; distinct() [dplyr package] 
library(tidyverse)
pubmed_hnscc_genes371 <- pubmed_hnscc_genes %>% dplyr::distinct(EntrezGene, .keep_all = TRUE)
# *** There is n= 371 unique HNSCC related genes, Tex curated.
save(pubmed_hnscc_genes371, file=file.path(path_cohort, "pubmed_hnscc_genes371.Rda"))
load(file=file.path(path_cohort, "pubmed_hnscc_genes371.Rda")) # as pubmed_hnscc_genes371
# TEX101: https://www.ncbi.nlm.nih.gov/pubmed/?term=23481573%5Buid%5D
 
# pubmed_hnscc_genes [2020/02/29 <- <- <- ] x自此更名?(n=371)
#X pubmed_hnscc_genes <- pubmed_hnscc_genes371 # for supplementary .csv file
write.csv(pubmed_hnscc_genes371[, -1], file=file.path(path_cohort, "pubmed_hnscc_genes371.csv"), quote = F,  row.names = F)

# n=635, with article frequency
pubmed_articles <- table(pubmed_hnscc_genes$EntrezGene) # EGFR, TP53 == p53 :-), EntrezGene 7157
pubmed_articles_5 <- pubmed_articles[pubmed_articles>5]
##pubmed_articles_5[[4]] <- pubmed_articles_5[[4]] + pubmed_articles_5[[9]]
##pubmed_articles_5[[8]] <- pubmed_articles_5[[8]] + pubmed_articles_5[[11]]
#x pubmed_articles_5 <- data.frame(cbind(pubmed_articles_5))[c(1:8, 10), ]
##pubmed_articles_5 <- as.data.frame(pubmed_articles_5)[c(1:8, 10), ]
pubmed_articles_5 <- as.data.frame(pubmed_articles_5)
colnames(pubmed_articles_5) <- c("EntrezGene", "Freq")

Entrez_Symbol <- subset(pubmed_hnscc_genes371, select=c(GeneSymbol, EntrezGene), (pubmed_hnscc_genes371$EntrezGene %in% pubmed_articles_5$EntrezGene))
  # pubmed_hnscc_genes371$GeneSymbol[pubmed_hnscc_genes371$EntrezGene %in% pubmed_articles_5$EntrezGene])
pubmed_articles_5 <- merge(pubmed_articles_5, Entrez_Symbol, by = "EntrezGene")
pubmed_articles_5 <- pubmed_articles_5[order(-pubmed_articles_5$Freq), ] #sorting by order (descending)
#pubmed_articles_5vector <- as.vector(unlist(pubmed_articles_5$Freq))
pubmed_articles_5vector <- pubmed_articles_5[["Freq"]]
# barplot needs a vector or table (freq)
# tiff(file=paste(TCGA_cohort, "_genes_embase5_barplot.tiff"), units="cm", width=5, height=5, res=300)
par(mar = c(7, 4, 2, 2) + 0.2) #add room for the rotated labels
end_point = 0.5 + nrow(pubmed_articles_5) + nrow(pubmed_articles_5)-1 #this is the line which does the trick (together with barplot "space = 1" parameter)
barplot(pubmed_articles_5vector, col="grey50",
        xlab="", horiz=F, space=1, names.arg="", cex.names=0.5,
        ylim=c(0, 5+max(pubmed_articles_5vector)))
# ***barplot has bug (vector checking is not well)
# name of bar label: names.arg (taken from colname of vector)
title(main=paste(TCGA_cohort, "genes published articles, until ", Sys.Date(), ""), 
      ylab="Frequence", )
#rotate 60 degrees, srt=60
#bar_label_index <- which(pubmed_hnscc_genes$EntrezGene==rownames(pubmed_articles_5))
text(seq(1.5, end_point, by=2), par("usr")[3]-0.25, 
     srt = 60, adj= 1, xpd = TRUE,
     labels = paste(pubmed_articles_5$GeneSymbol), cex=0.8)

#dev.off()
# R4> em_Index # n=75, removal those published genes
em_Index <- which(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$gene_id %in% pubmed_hnscc_genes371$GeneSymbol)
# [1]   26   67  206  284  427  484  530  609  809  971 1020 1067 1137 1178 1193
# [16] 1365 1404 1423 1501 1518 1523 1698 1819 1881 2036 2104 2197 2212 2360 2380
# [31] 2490 2696 2706 2759 2808 2926 2953 2992 3287 3392 3550 3700 3796 3848 3854
# [46] 3938 3985 4076 4095 4125 4213 4246 4352 4360 4604 4667 4671 4677 4731 4742
# [61] 4849 4900 4948 5034 5095 5550 5582 5657 5659 5699 5708 6018 6094 6153 6263
HNSCC_OS_marginS_THREE_pvalue005_noCancerGene <- HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[-em_Index, ] # removal
# n=6293 genes
# (ok) Excluded the HNSCC cancer driver genes, list retrieved from Embase
# ?BioXpress*.csv # data from https://hive.biochemistry.gwu.edu/cgi-bin/prd/bioxpress/servlet.cgi
# <done> [2019/07/20] anniversary for discovery of GLUT4 in HNSCC cell line [2014/07/20]
save(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, file=file.path(path_cohort, "HNSCC_OS_marginS_THREE_pvalue005_noCancerGene.Rda"))
load(file=file.path(path_cohort, "HNSCC_OS_marginS_THREE_pvalue005_noCancerGene.Rda"))





#> uni/multi_HR plot, bad_FC ####
# there is no "filtering", plot is just for demo
#x HNSCC_OS_marginS_THREE_pvalue005 <- HNSCC_OS_marginS_THREE_pvalue005_noCancerGene

bad_FC <- 1.5 # cutoff for fold change
good_FC <- 0.5 # there is also good guy
# ** choice proper adjusted P-value: "p_value" "p_value_adj_bonferroni" "p_value_adj_FDR" 
raw <- 
  readline("raw P-value(p), Bonferroni P-value(b) or FDR P-value(f) -- process run: ")
select_p <- function(x) {
  switch(x,
         p = "uncorrected P-value",
         b = "Bonferroni P-value",
         f = "FDR P-value",
         stop("Unknown input")
  )
}
select_pos <- function(x) { # col of 3 P-values
  switch(x,
         p = 1,
         b = 2,
         f = 3,
         stop("Unknown input")
  )
}
pvalue3 <- select_p(raw)
pvalue_pos <- select_pos(raw)
print(c(pvalue3, pvalue_pos))
rm(raw)


# univariate HR by Cox
#tiff("Rplot07_cox_uniHR.tiff", units="cm", width=5, height=5, res=300) 
# x# saving as .tiff (by tiff()) => x y axis is too small in .tiff
# adds a loess regression smoother to a scatter plot by smoothScatter
#smoothScatter(p_value, uni_HR, nrpoints = 10000, cex=0.6, xlab="P-value (KM survival)", ylab="Cox's Harzard Ratios (univariate)", bandwidth=c(7, 7), nbin=20)
# ** choice proper adjusted P-value: "p_value" "p_value_adj_bonferroni" "p_value_adj_FDR" 
# *** taking uni_P-value < 0.05, then plot; n=6204
HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_uni <- HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[
  HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$uni_P_value < 0.05, ] # a temp variable

plot(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_uni[ , pvalue_pos+1], 
     HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_uni$uni_HR, type="p", main=paste(TCGA_cohort, "Cox's Harzard Ratios (univariate)\n", "versus", pvalue3), 
     ylab="Cox's HR", xlab="P-value (KM survival)", log="x", cex=0.3, , cex.main=1.2, cex.lab=1.3, cex.axis=1)
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05))) #1=below, 2=left, 3=above and 4=right
x0_p <- min(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_uni[ , pvalue_pos+1])
x1_p <- max(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_uni[ , pvalue_pos+1])
y0_p <- min(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_uni$uni_HR)
y1_p <- max(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_uni$uni_HR)
segments(x0=x0_p, y0=bad_FC, x1=x1_p, y1=bad_FC, lty=2, col="red") # no abline()
segments(x0=alpha_HNSCC, y0=y0_p, x1=alpha_HNSCC, y1=y1_p, lty=2, col="blue")
segments(x0=x0_p, y0=good_FC, x1=x1_p, y1=good_FC, lty=2, col="green")
#segments(x0=0.5e-7, y0=150, x1=1e-2, y1=150, lty=2, col="blue")
#segments(x0=Bonferroni_cutoff, y0=0, x1=Bonferroni_cutoff, y1=170, lty=2, col="red") # 5.31011e-06
# run a logistic regression model (categorical 0 vs 1)
#g <- glm(uni_HR ~ p_value, family=poisson, data=HNSCC_OS_marginS_THREE_pvalue005)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
legend("topright", legend=c(paste("HR at", bad_FC), paste("HR at", good_FC), paste("P-value at", alpha_HNSCC)), lty=2:2, col=c("red", "green", "blue"), cex=1) # box and font size
# text(500,900, labels=paste(TCGA_cohort, "KM P-value plot"), cex = 0.60) # (0,0) on bottom_left corner
# figure legend: logistic regression, LR, by Generalized linear model, glm
#dev.off()


# multivariate HR by Cox
# **plot multi_HR, HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, under Bonferroni
# tiff("Rplot07_cox_multiHR.tiff")
# ** choice proper adjusted P-value: "p_value" "p_value_adj_bonferroni" "p_value_adj_FDR" 
#x pvalue3 <- "p_value"
# # *** keeping uni_P-value < 0.05, then plot (n=4942)
HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_multi <- HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[
  HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$multi_P_value < 0.05, ]  # a temp variable

plot(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_multi[ , pvalue_pos+1], 
     HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_multi$multi_HR, type="p", main=paste(TCGA_cohort, "Cox's Harzard Ratios (multivariate)\n", "versus", pvalue3), 
     ylab="Cox's HR", xlab="P-value (KM survival)", log="x", cex=0.3, , cex.main=1.2, cex.lab=1.3, cex.axis=1)
#axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05))) #1=below, 2=left, 3=above and 4=right
x0_p <- min(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_multi[ , pvalue_pos+1])
x1_p <- max(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_multi[ , pvalue_pos+1])
y0_p <- min(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_multi$multi_HR)
y1_p <- max(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene_multi$multi_HR)
segments(x0=x0_p, y0=bad_FC, x1=x1_p, y1=bad_FC, lty=2, col="red") # no abline()
segments(x0=alpha_HNSCC, y0=y0_p, x1=alpha_HNSCC, y1=y1_p, lty=2, col="blue")
segments(x0=x0_p, y0=good_FC, x1=x1_p, y1=good_FC, lty=2, col="green")
#segments(x0=0.5e-7, y0=150, x1=1e-2, y1=150, lty=2, col="blue")
#segments(x0=Bonferroni_cutoff, y0=0, x1=Bonferroni_cutoff, y1=170, lty=2, col="red") # 5.31011e-06
# run a logistic regression model (categorical 0 vs 1)
#g <- glm(uni_HR ~ p_value, family=poisson, data=HNSCC_OS_marginS_THREE_pvalue005)
# (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
legend("topright", legend=c(paste("HR at", bad_FC), paste("HR at", good_FC), paste("P-value at", alpha_HNSCC)), lty=2:2, col=c("red", "green", "blue"), cex=1) # box and font size
# figure legend: logistic regression, LR, by Generalized linear model, glm


# Figure 2 ####
# HNSCC_OS_marginS_THREE_pvalue005_noCancerGene => for manuscript Figure 2
# By uncorrected P-value below 0.05, we selected 967 genes which HR is greater than 1.5 or less than 0.5 
# (see Figure \ref{fig:figure2a} univariate,
HNSCC_OS_marginS_uni1p5_Figure2 <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, (p_value < alpha_HNSCC) & (uni_P_value < 0.05) & (uni_HR >= bad_FC), 
                                          select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=942
HNSCC_OS_marginS_uni0p5_Figure2 <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, (p_value < alpha_HNSCC) & (uni_P_value < 0.05) & (uni_HR <= good_FC), 
                                          select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=25
# ... and Figure \ref{fig:figure2b} multivariate, respectively).
HNSCC_OS_marginS_multi1p5_Figure2 <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, (p_value < alpha_HNSCC) & (uni_P_value < 0.05) & (uni_HR >= bad_FC), 
                                          select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=942
HNSCC_OS_marginS_multi0p5_Figure2 <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, (p_value < alpha_HNSCC) & (uni_P_value < 0.05) & (uni_HR <= good_FC), 
                                          select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=25

# x Venn_marginS_: cross matching HR>1.8 of Uni+Multi, HR<0.5 of Uni+Multi
#library(venn) #https://cran.r-project.org/web/packages/venn/venn.pdf
# venn(x, snames = "", ilabels = FALSE, counts = FALSE, ellipse = FALSE,
#      zcolor = "bw", opacity = 0.3, size = 15, cexil = 0.6, cexsn = 0.85,
#      borders = TRUE, ...)
#library(gplots)
# To get the list of gene present in each Venn compartment we can use the gplots package
#library(gplots) # capture the list of genes from venn

# x *** how about Z-score? included in HNSCC_OS_marginS_THREE_pvalue005_1475
# x use Z-score at the end of our analysis




#{> [pickup1] (bad guy genes)#### 
# from HNSCC_OS_marginS_THREE_pvalue005_noCancerGene # ** uni_p multi_p 尚未篩選
# Broader gene candidate (first 100): Cox HR (>1.5 or >=2.5), bad_FC fold change 
# (uni_P_value <= 0.05) & (multi_P_value <= 0.05)
# Bonferroni_cutoff = 5.31011e-06 is enoughly restricted in this cohort
bad_FC <- bad_FC # 1.5; 1.6 or 1.8 # it was decided at [# P-values plot uni_HR] section
# #HNSCC_OS_marginS_uni_CoxHR2p5 <- subset(HNSCC_OS_marginS_THREE_pvalue005, (p_value <= Bonferroni_cutoff) & (uni_HR >= 2.5), 
# #                                       select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# # xmulti_HR >=1.5 # & (multi_P_value <= 0.05) 
# HNSCC_OS_marginS_unimulti_CoxHR2p5_uncorrected <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene,  (p_value < alpha_HNSCC) & (multi_P_value < 0.05) & (multi_HR >= bad_FC) & (uni_P_value < 0.05) & (uni_HR >= bad_FC), 
#                                                          select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# #x n=679; uncorrected P-value < 0.05

# n=6293
# univariate p_value_adj_bonferroni  Bonferroni P-value < 0.05 才是王道
HNSCC_OS_marginS_uni_CoxHR2p5_Bonf <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, (p_value_adj_bonferroni < alpha_HNSCC) & (uni_P_value < 0.05) & (uni_HR >= bad_FC), 
                                        select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=10 同下
# [1] "DKK1"    "CAMK2N1" "STC2"    "PGK1"    "SURF4"   "USP10"  
# [7] "NDFIP1"  "FOXA2"   "STIP1"   "DKC1"

# multivariate Bonferroni P-value < 0.05 才是王道
HNSCC_OS_marginS_multi_CoxHR2p5_Bonf <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene,  (p_value_adj_bonferroni < alpha_HNSCC) & (multi_P_value < 0.05) & (multi_HR >= bad_FC), 
                                          select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# # n=10; 同上
# [1] "DKK1"    "CAMK2N1" "STC2"    "PGK1"    "SURF4"   "USP10"  
# [7] "NDFIP1"  "FOXA2"   "STIP1"   "DKC1"


# FDR uni/multi
HNSCC_OS_marginS_uni_CoxHR2p5_FDR <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, (p_value_adj_FDR < alpha_HNSCC) & (uni_P_value < 0.05) & (uni_HR >= bad_FC), 
                                        select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
HNSCC_OS_marginS_multi_CoxHR2p5_FDR <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene,  (p_value_adj_FDR <= alpha_HNSCC) & (multi_P_value < 0.05) & (multi_HR >= bad_FC), 
                                          select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=942 and n=845; FDR P-value < 0.05; 跟沒有校正的，基因一樣 順序不一樣
# comparison...
# # R4> HNSCC_OS_marginS_multi_CoxHR2p5_FDR$gene_id %in% HNSCC_OS_marginS_uni_CoxHR2p5_FDR$gene_id
#    Mode   FALSE    TRUE 
# logical     166     679 交集 n=679
library(venn)
venn::venn(list(HNSCC_OS_marginS_uni_CoxHR2p5_FDR, HNSCC_OS_marginS_multi_CoxHR2p5_FDR),
     snames=c("univariate", "multivariate"),
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "red, orange", opacity = 0.4, size = 15, cexil = 0.7, cexsn = 0.3, borders = TRUE)
#      predefined colors if "style"
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# meta-language 1 0 or -
title <- c(paste(TCGA_cohort, "Cox's Hazard ratio, FC >= ", bad_FC, ""), paste("( while KM P-Value (FDR) <", signif(alpha_HNSCC, 3), ")")) #, collapse = "\n")
#coords <- unlist(getCentroid(getZones(venn_HR2p5, snames="uni_CoxHR>=2p5, multi_CoxHR>=2p5")))
# coords[1], coords[2], 
text(500,900, labels = title[1], cex = 1.5) # (0,0) on bottom_left corner
text(500,855, labels = title[2], cex = 1) 
library(gplots)
tmp_bad <- venn(list(HNSCC_OS_marginS_uni_CoxHR2p5_FDR, HNSCC_OS_marginS_multi_CoxHR2p5_FDR),
                names=c("univariate", "multivariate"), show.plot=F) #library(gplots); the group count matrix alone
isect_HR2p5_FDR <- attr(tmp_bad, "intersections")

# # R4> HNSCC_OS_marginS_multi_CoxHR2p5_FDR$gene_id %in% HNSCC_OS_marginS_multi_CoxHR2p5_uncorrected$gene_id
# Mode    TRUE 
# logical    1106 

# => 判斷 ***take genes by Bonferroni P-value correction: 
# HNSCC_OS_marginS_uni_CoxHR2p5_Bonf
# HNSCC_OS_marginS_multi_CoxHR2p5_Bonf
#save(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue005_noCancerGene.Rda", sep=""))) # in HNSCC_OS_marginS_THREE_pvalue005_noCancerGene (final results)
#isect_HR2p5 <- HNSCC_OS_marginS_multi_CoxHR2p5_Bonf
save(HNSCC_OS_marginS_uni_CoxHR2p5_Bonf, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "uni_CoxHR2p5_Bonf.Rda", sep=""))) # as HNSCC_OS_marginS_uni_CoxHR2p5_Bonf (final bad guy candidates)
save(HNSCC_OS_marginS_multi_CoxHR2p5_Bonf, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "multi_CoxHR2p5_Bonf.Rda", sep=""))) # as HNSCC_OS_marginS_multi_CoxHR2p5_Bonf (final bad guy candidates) with *** final Bonferroni
# go to pickup2 (good guy) then Export r2excel 2020



#FDR: venn1 diagram of HR>=1.5 of Uni & Multi ###
#for list of genes by grouping; 
venn_HR2p5 <- list(HNSCC_OS_marginS_uni_CoxHR2p5$gene_id, HNSCC_OS_marginS_multi_CoxHR2p5$gene_id)
# no cutoff by Bonferroni_cutoff
names_HR2p5 <- c(paste("Cox's HR(univariate) >=", bad_FC), paste("Cox's HR(multivariate) >=", bad_FC))
library(gplots)
tmp_bad <- venn(venn_HR2p5, names=names_HR2p5, show.plot=F) #library(gplots); the group count matrix alone
isect_HR2p5 <- attr(tmp_bad, "intersections")
#isect_HR2p5
detach(package:gplots)


#x  https://www.statmethods.net/advgraphs/parameters.html
library(venn)
# https://cran.r-project.org/web/packages/venn/venn.pdf
tiff("Rplot09_venn_HR1.5.tiff", units="cm", width=5, height=5, res=300) 
# # saving as .tiff (by tiff())
venn::venn(venn_HR2p5, snames=names_HR2p5,
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "red, deeppink", opacity = 0.6, size = 15, cexil = 0.7, cexsn = 0.3, borders = TRUE)
#      predefined colors if "style"
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# meta-language 1 0 or -
title <- c(paste(TCGA_cohort, "survival analysis"), paste("(KM P-Value <=", signif(alpha_HNSCC, 3), ")")) #, collapse = "\n")
#coords <- unlist(getCentroid(getZones(venn_HR2p5, snames="uni_CoxHR>=2p5, multi_CoxHR>=2p5")))
# coords[1], coords[2], 
text(500,900, labels = title[1], cex = 0.60) # (0,0) on bottom_left corner
text(500,855, labels = title[2], cex = 0.40) 
# n=47
dev.off() # saving instead of showing

#https://stackoverflow.com/questions/43324180/adding-legend-to-venn-diagram
#legend("top", legend=c("B:multi_Cox HR >=2.5", "A:uni_Cox HR >=2.5"), lty=1:2, col=c("blue","red"), cex=0.7) # box and font size
# figure legend: 




# {> [pickup2]  (good guy genes) ####
# from HNSCC_OS_marginS_THREE_pvalue005_noCancerGene
#* Cox HR <0.5 # 
good_FC <- good_FC  # 0.5
#HNSCC_OS_marginS_uni_CoxHR0p5 <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, (p_value <= Bonferroni_cutoff) & (uni_P_value <= 0.05) & (multi_P_value <= 0.05) & (uni_HR <0.0), 
#                                       select=c(number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=0 while (uni_P_value <= 0.05) & (multi_P_value <= 0.05) & (uni_HR <0.0)
# 
#HNSCC_OS_marginS_uni_CoxHR0p5 <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, (p_value <= alpha_HNSCC) & (uni_P_value <= 0.05) & (uni_HR <good_FC), 
#                                        select=c(gene_id, p_value, z_score, uni_HR, uni_P_value, multi_HR, multi_P_value))
# uni multi candidate 相近:
# univariate p_value_adj_bonferroni  Bonferroni P-value < 0.05 才是王道
HNSCC_OS_marginS_uni_CoxHR0p5_Bonf <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene, (p_value_adj_bonferroni < alpha_HNSCC) & (uni_P_value < 0.05) & (uni_HR <= good_FC), 
                                             select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# n=12 不同下
# [1] "ZNF557"    "ZNF266"    "IL19"      "MYO1H"     "FCGBP"    
# [6] "LOC148709" "EVPLL"     "PNMA5"     "KIAA1683"  "SLC26A9"  
# [11] "NPB"       "FRMD4A"  

# multivariate Bonferroni P-value < 0.05 才是王道
HNSCC_OS_marginS_multi_CoxHR0p5_Bonf <- subset(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene,  (p_value_adj_bonferroni < alpha_HNSCC) & (multi_P_value < 0.05) & (multi_HR <= good_FC), 
                                               select=c(gene_id, p_value, p_value_adj_bonferroni, p_value_adj_FDR, uni_HR, uni_P_value, multi_HR, multi_P_value))
# # n=13; 不同上
# [1] "ZNF557"    "ZNF266"    "IL19"      "MYO1H"     "FCGBP"    
# [6] "LOC148709" "EVPLL"     "PNMA5"     "KIAA1683"  "GPR15"    
# [11] "NPB"       "CALML5"    "PPFIA3" 
save(HNSCC_OS_marginS_uni_CoxHR0p5_Bonf, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "uni_CoxHR0p5_Bonf.Rda", sep=""))) # as HNSCC_OS_marginS_uni_CoxHR0p5_Bonf (final good guy candidates)
save(HNSCC_OS_marginS_multi_CoxHR0p5_Bonf, file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "multi_CoxHR0p5_Bonf.Rda", sep=""))) # as HNSCC_OS_marginS_multi_CoxHR0p5_Bonf (final good guy candidates) with *** final Bonferroni
# go to Export r2excel 2020 (line 1798)


# FDR uni/multi

# FDR: venn2 diagram of HR < 0.5 of Uni & Multi ###
#for list of genes by grouping; library(gplots)
venn_HR0p5 <- list(HNSCC_OS_marginS_uni_CoxHR0p5$gene_id, HNSCC_OS_marginS_multi_CoxHR0p5$gene_id)
names_HR0p5 <- c(paste("Cox's HR(univariate) <", good_FC), paste("Cox's HR(multivariate) <", good_FC))
library(gplots)
tmp_good <- venn(venn_HR0p5, names=names_HR0p5, show.plot=F) #library(gplots); the group count matrix alone
isect_HR0p5 <- attr(tmp_good, "intersections")
detach(package:gplots)

library(venn)
# venn(venn_HR0p5, snames=names_HR0p5,
#      ilabels = T, counts = T, ellipse = FALSE,  zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 0.6, cexsn = 0.85, borders = TRUE)
#      predefined colors if "style"
tiff("Rplot09_venn_HR0.5.tiff", units="cm", width=5, height=5, res=300) 
# # saving as .tiff (by tiff())
venn(venn_HR0p5, snames=names_HR0p5,
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 0.7, cexsn = 0.3, borders = TRUE)
#      predefined colors if "style"
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# meta-language 1 0 or -
title <- c(paste(TCGA_cohort, "survival analysis"), paste("(KM P-Value <=", signif(alpha_HNSCC, 3), ")")) #, collapse = "\n")
#coords <- unlist(getCentroid(getZones(venn_HR2p5, snames="uni_CoxHR>=2p5, multi_CoxHR>=2p5")))
# coords[1], coords[2], 
text(500,900, labels = title[1], cex = 0.60) # (0,0) on bottom_left corner
text(500,855, labels = title[2], cex = 0.40) 
# n=43
dev.off() # saving instead of showing

# signif(Bonferroni_cutoff, 3)


#x #library(VennDiagram)
# VENN.LIST <- list(HNSCC_OS_marginS_uni_CoxHR0p5$gene_id, HNSCC_OS_marginS_multi_CoxHR0p5$gene_id)
# #venn.plot <- venn.diagram(VENN.LIST , NULL, fill=c("darkmagenta", "darkblue"), alpha=c(0.5,0.5), cex = 2, cat.fontface=4, category.names=c("A", "B"), main="HNSCC: Cox HR <0.5")
# # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
# #grid.draw(venn.plot)
# # We can summarize the contents of each venn compartment, as follows:
# # in 1) ConditionA only, 2) ConditionB only, 3) ConditionA & ConditionB
# lapply(isect_HR0p5, head)
#}

# # Embase HNSCC genes was removed; saved genes in left, right and intersection
# save(isect_HR2p5, isect_HR0p5, 
#      file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "candidates_Venn.Rda", sep="")))
#} venn the end




##x Export r2excel and .Rda 2019 ####
# 直接 jump to Export 2020
# _marginS_ [2020/03/04] z_cut 要刪去
# Embase HNSCC genes was removed; no venn is necessary
load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "multi_CoxHR2p5_Bonf.Rda", sep=""))) # in HNSCC_OS_marginS_multi_CoxHR2p5_Bonf (final candidates)
# as HNSCC_OS_marginS_multi_CoxHR2p5_Bonf
isect_HR2p5 <- HNSCC_OS_marginS_multi_CoxHR2p5_Bonf # n=9
save(isect_HR2p5, 
     file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "candidates_Venn.Rda", sep="")))
# sink() for .xlsx export as well :-) https://stackoverflow.com/questions/34038041/how-to-merge-multiple-data-frame-into-one-table-and-export-to-excel
# HNSCC_OS_marginS_candidates_Venn.xlsx: show up Bonferroni_cutoff and 5e-6 (expression style).
# *** loading 素材們 
#load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue005_zcut.Rda", sep=""))) # in HNSCC_OS_marginS_pvalue005_zcut, Z-score >= zcut (e.x. 0.8)
#load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep=""))) # in candidate_cox (a list), candidate_samplie, candidate_cox, n_percent_Bonferroni
load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "candidates_Venn.Rda", sep=""))) # isect_HR2p5 <- HNSCC_OS_marginS_multi_CoxHR2p5_Bonf; [x isect_HR0p5]

library("xlsx") # as well as rJava: ok
# # installing curl directly then r2excel
# install.packages("devtools")
# library(devtools)
# install_github("kassambara/r2excel", force=T) # the both ways are OK; reinstall it
# # or
# install.packages(c("cli", "gh", "usethis", "githubinstall"))
# library("githubinstall")
# githubinstall("r2excel")
library("r2excel")
# ***going to R2Excel (2020) 



#xxx Create an Excel workbook. Both .xls and .xlsx file formats can be used.
filenamex <- paste(TCGA_cohort, "_OS", marginTag, "candidates_Venn", ".xlsx", sep = "") 
# "HNSCC_OS_marginS_candidates_Venn.xlsx"
wb <- createWorkbook(type="xlsx")

# Create a sheet in that workbook
sheet <- xlsx::createSheet(wb, sheetName = paste("Survival_candidates"))
# [add data row by row, start from column 2]
#+++++++++++++++++++++++++++++++
## Add paragraph : Author
library("tis") # by Brian Salzer
# today(), arg must be ti, tis, ts, tif, or tifName
select_title <- function(x) {
  switch(x,
         "_marginS_" = "The cohort with surgical margins status (0 or 1)",
         "_marginFree_" = "The cohort with surgical margins free from tumor (0)",
         "_marginPlus_" = "The cohort with surgical margins involving tumor (1)",
         stop("Unknown input")
  )
}
title_candidates_Venn_xlsx <- select_title(marginTag)
author <- paste("Reported by Tex Li-Hsing Chi. \n",
                "tex@gate.sinica.edu.tw \n", title_candidates_Venn_xlsx, "\n", Sys.Date(), sep="")
cat(author)
xlsx.addParagraph(wb, sheet, value=author, isItalic=TRUE, colSpan=5, 
                  rowSpan=4, fontColor="darkgray", fontSize=24)
xlsx.addLineBreak(sheet, 3)
# header
# Z-score:a normalized number (frequency) of OS P-values, which KM P-value < 0.05, in each gene
xlsx.addHeader(wb, sheet, value=paste("Table 1. The top 30 candiate genes expressed in ", TCGA_cohort,
                                      " (ranking by KM P-value, selected by Z-score >= ", signif(zcut, 2), ") ", "\n" , sep=""),
               level=5, color="black", underline=0) # nrow(HNSCC_OS_marginS_pvalue005_zcut), n=6601
#xlsx.addHeader(wb, sheet, value=paste("Cutoff at ", round(cutoff1, 3), " (", percent(surv_OS1$n[1]/(surv_OS1$n[1]+surv_OS1$n[2])), ")", sep = ""),
#               level=5, color="red", underline=0) # total n is taken from surv_OS1$n

xlsx.addLineBreak(sheet, 1) # add one blank line

#xlsx.addTable(wb, sheet, data = t(data.frame(c(paste(geneName, "expression"), "", paste(geneName, "expression"), "", "(Optimised)"))), fontSize=12, startCol=4,
#              fontColor="darkblue", row.names = F, col.names = F) #, colSpan=1, rowSpan=1)

# x a candidate genes table
# 要修
candidates_alpha_pvalue <- HNSCC_OS_marginS_multi_CoxHR2p5_Bonf # after Bonferroni correcion
# > colnames(HNSCC_OS_marginS_pvalue005_zcut)
# [1] "gene_id" "p_value" "z_score" rename to "Gene_id", "P_value", "Z_score" (why?)
colnames(candidates_alpha_pvalue)[1:3] <- c("Gene_id", "P_value", "Bonferroni P_value")

# # # no more Bonferroni Bonferroni_cutoff; [, c(1,2,4)]
# attach(HNSCC_OS_marginS_pvalue005_zcut)
# candidates_alpha_pvalue <- HNSCC_OS_marginS_pvalue005_zcut[which(P_value<=alpha_HNSCC), ]
# detach(HNSCC_OS_marginS_pvalue005_zcut)

attach(candidates_alpha_pvalue) # removal of as.factor (P_value)
# in scientific notation: formatC of [, c(2)]
candidates_alpha_pvalue$P_value <- formatC(as.numeric(as.character(P_value)), format = "e", digits = 2)
candidates_alpha_pvalue$`Bonferroni P_value` <- formatC(as.numeric(as.character(`Bonferroni P_value`)), format = "e", digits = 2)
#candidates_alpha_pvalue$Z_score <- signif(Z_score, digits=4)
candidates_alpha_pvalue <- candidates_alpha_pvalue[order(`Bonferroni P_value`), ] #sorting by order(ascending)
detach(candidates_alpha_pvalue)

# export top 30 genes [1:30,], out of 1702 genes
xlsx.addTable(wb, sheet, data = candidates_alpha_pvalue[, ], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = TRUE)

xlsx.addLineBreak(sheet, 5)  # add two blank lines

# Export (no Z-score)/P-Value plot by plotFunction() ####
xlsx.addLineBreak(sheet, 1)

#x Define function
# xlsx.addPlot.candidates<-function( wb, sheet, startRow=NULL,startCol=2,
#                                    width=480, height=480,... )
# { # plot of Z-score vs P-value digram from "the summary"
#   
#   # png(filename = "plot.png", width = width, height = height,...)
#   # #{
#   # # plot(OS.km, lty=1, xscale=12, xmax=60, col=c("blue","red"), sub=paste("P-value =", p_OS1), main=paste("OS in TCGA ", TCGA_cohort, "(n=", surv_OS1$n[1]+surv_OS1$n[2],")/", geneName), ylab="Percent Survival", xlab="Years")
#   # # legend("topright", legend=c(paste("low(",surv_OS1$n[1], ")"), paste("high(",surv_OS1$n[2], ")")), lty=1:2, col=c("blue","red"))
#   # load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalue005_sorted.Rda", sep=""))) 
#   # # is HNSCC_OS_marginS_pvalue005_sorted
#   # # "swap data" for following code running
#   # 
#   # attach(HNSCC_OS_marginS_pvalue005_sorted) # Z-score
#   # plot(p_value, z_score, type="p", ylab="Z-score", xlab="P-value (KM survival analysis)", log="x") # log scale x or y
#   # #axis(side=3, at=c(1e-7, 1e-6, 1e-5, 0.01, 0.05)) #1=below, 2=left, 3=above and 4=right
#   # abline(h=0.6, lty=2, col="blue")
#   # abline(v=alpha_HNSCC, lty=2, col="red") #
#   # # run a logistic regression model (categorical 0 vs 1)
#   # #g <- glm(number_01 ~ p_value, family=poisson, data=HNSCC_OS_marginS_pvalue005_sorted)
#   # # (in this case, generalized linear model with log link)(link = "log"), poisson distribution
#   # #curve(predict(g, data.frame(p_value = x), type="response"), add=TRUE, col="blue") # draws a curve based on prediction from logistic regression model
#   # legend("bottomleft", legend=c(paste("Z-score at 0.6"), paste("P-value at", alpha_HNSCC)), lty=2:2, col=c("blue","red"), cex=0.9) # box and font size
#   # # figure legend: logistic regression, LR, by Generalized linear model, glm
#   # detach(HNSCC_OS_marginS_pvalue005_sorted)
#   # #}
#   # dev.off()
#   # 
#   #Append plot to the sheet
#   if(is.null(startRow)){
#     rows<- getRows(sheet) #list of row object
#     startRow=length(rows)+1
#   } 
#   # append to the file created previously
#   addPicture("Rplot10_Zscore_Pvalue.png", sheet=sheet,  startRow = startRow, startColumn = startCol) 
#   xlsx.addLineBreak(sheet, round(width/20)+1)
#   #res<-file.remove("plot.png") # Z-score vs P-value plot
#   #file.rename("plot.png", "Zscore_Pvalue_scatter.png")
# } # Define function
# 
# # calling and convert tiff to png
# system("convert Rplot10_Zscore_Pvalue.tiff Rplot10_Zscore_Pvalue.png")
# currentRow <- xlsx.addPlot.candidates(wb, sheet) # startRow + a plot
# #print(paste("The z-score summary plot:", filenamex,"was exported successfully."))


# # Export Cox uni/multi HR vs P-value plot ####
# Define function
xlsx.addPlot.candidates<-function( wb, sheet, file="plot.png", startRow=NULL,startCol=2,
                                   width=480, height=480,... )
{ # plot of Cox uni/multi HR vs P-value plot
  #Append plot to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  # append to the file created previously, add 
  # {
  # Cox's Harzard Ratios (univariate) plot [Rplot07_cox_uniHR.tiff]
  # Cox's Harzard Ratios (multivariate) plot [Rplot07_cox_multiHR.tiff]
  # }
  addPicture(file, sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
} # Define function

# calling
# converting .tiff to .png for export (from Version: ImageMagick 6.9.7-4)
#system("convert Rplot07_cox_uniHR.tiff Rplot07_cox_uniHR.png", intern = FALSE)
# Cox's Harzard Ratios (univariate) plot [Rplot07_cox_uniHR.tiff]
#currentRow <- xlsx.addPlot.candidates(wb, sheet, file="Rplot07_cox_uniHR.png") # startRow + a plot
# Cox's Harzard Ratios (multivariate) plot [Rplot07_cox_multiHR.tiff] => cox_multi_Bonf_KMPvalue.tiff
# cox_multi_Bonf_KMPvalue.tiff is Bonferroni KM P-value plot in multivariate analysis
system("convert cox_multi_Bonf_KMPvalue.tiff cox_multi_Bonf_KMPvalue.png", intern = FALSE) # convert is a duplicate
currentRow <- xlsx.addPlot.candidates(wb, sheet, file="cox_multi_Bonf_KMPvalue.png") # startRow + a plot




#x Export Venn1 (bad guy genes) by plotFunction() ####
xlsx.addLineBreak(sheet, 1)

# Define function
xlsx.addPlot.venn<-function( wb, sheet, startRow=NULL,startCol=2,
                             width=480, height=480,... )
{ # plot of venn digram: Rplot09_venn_HR1.5.tiff"
  
  # png(filename = "plot1.png", width = width, height = height,...)
  # #{
  # venn(venn_HR2p5, snames=names_HR2p5,
  #      ilabels = T, counts = T, ellipse = FALSE,  zcolor = "red, deeppink", opacity = 0.6, size = 15, cexil = 3, cexsn = 0.85, borders = TRUE)
  # #      predefined colors if "style"; cexil = 0.6 (default text size of counts)
  # title <- paste(c(TCGA_cohort, " survival analysis", "KM P-Value <= ", signif(alpha_HNSCC, 3)), sep = "", collapse = "\n")
  # text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner
  # 
  # #}
  # dev.off()
  
  
  #Append Rplot09_venn_HR1.5.tiff to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 

    # Add the file created previously
  system("convert -resize 2048X1536 Rplot09_venn_HR1.5.tiff Rplot09_venn_HR1.5.png", intern = FALSE)
  addPicture("Rplot09_venn_HR1.5.png", sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
  #res<-file.remove("plot1.png")
  
} # end of Define function

# calling
#detach(package:gplots)
#Append saved "Rplot09_venn_HR1.5.tiff" to the sheet
#print(paste("The Venn diagram (bad guy) and the candidate genes on", filenamex, ", which was exported successfully."))
# converting .tiff to .png for export

library(venn)
xlsx.addPlot.venn(wb, sheet)



# (FDR: bad guy table) prognostic features of those genes in TCGA HNSCC cohort ####
# store KM_sig remrark as a Byte; converted as base10 in list_KM_sigBin
load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep=""))) # as candidate_cox (a list), since 2019/07/03; with KM_sig, Remark
# candidate_sample, candidate_cox, n_percent_Bonferroni
# isect_HR2p5 <- HNSCC_OS_marginS_multi_CoxHR2p5_Bonf # n=9
load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "candidates_Venn.Rda", sep=""))) # as isect_HR2p5

# # x load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005.Rda", sep="")))
# #load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "THREE_pvalue005_1475.Rda", sep="")))
# # Z-score > zcut e.x. 0.8: HNSCC_OS_marginS_pvalue005_zcut.Rda, n=1475
# #HNSCC_OS_marginS_THREE_pvalue005 <- HNSCC_OS_marginS_THREE_pvalue005_1475 # *** try it
# # as as HNSCC_OS_marginS_THREE_pvalue005, Bonferroni_cutoff
# 



# #.. bad guy gene candidate
# # => see "remark=1"; "which" or NOT "which", the order is wrong.
# # isect_HR2p5 is a list, which comes from line 2909
# # R4> isect_HR2p5
# # $`uni_Cox HR >=1:multi_Cox HR >=1`
# # [1] "CAMK2N1" "USP10"   "PGK1"    "SURF4"   "EFNB2"   "EIF2AK1" "STC2"   

# R4> isect_HR2p5 (FDR with venn)
# $`uni_Cox HR >=1:multi_Cox HR >=1`
# [1] "CAMK2N1" "USP10"   "PGK1"    "SURF4"   "EFNB2"   "EIF2AK1" "STC2"   
isect_HR2p5 <- isect_HR2p5_FDR
geneid_bad_uni_HR2p5 <- c(isect_HR2p5[[3]], isect_HR2p5[[1]]) # geneid_bad_uni_HR2p5 <- c(isect_HR2p5$`A:B`, isect_HR2p5$`A`) # from venn_HR2p5: A + A:B = group uni_
# whole_genome[which(whole_genome %in% geneid_bad_uni_HR2p5)]
#  [1] "CD1A"    "CENPN"   "COX7A2L" "ERCC6"   "HIGD1B"  "HMGN5"   "NEK6"    "NUF2"   
# [9] "PA2G4P4" "PCTP"    "PLOD2"   "RRAS2"   "SFXN1" 
candidate_bad_uni_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_uni_HR2p5), HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[which(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$gene_id %in% geneid_bad_uni_HR2p5), 2:8])
#as.binary(candidate_cox[which(whole_genome %in% geneid_bad_uni_HR2p5)]$KM_sig, n=7, logic=T) # store KM_sig remrark as a Byte
# can NOT use "whole_genome" any more
#x list_KM_sig <- candidate_cox[which(whole_genome %in% geneid_bad_uni_HR2p5)] 
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_bad_uni_HR2p5))] # a list (rownames has whole_genome's right position)
# km <- function(x) {print(list_KM_sig[[x]][8,9])}; lapply(c(1:10), km)
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_bad_uni_HR2p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8]) # 8 column, without tobacco feature
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as.factor [0 or 1]
}
candidate_bad_uni_HR2p5 <- cbind(candidate_bad_uni_HR2p5, list_KM_sigBin)
# ex. "NEK6" has bad impact on clinical T when expression is high

# => see "remark=1"                                                # number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value])
# # isect_HR2p5 is a list
geneid_bad_multi_HR2p5 <- c(isect_HR2p5[[3]], isect_HR2p5[[2]]) # from venn_HR2p5
candidate_bad_multi_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_multi_HR2p5), HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$gene_id %in% geneid_bad_multi_HR2p5), 2:7])
# as.binary
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_bad_multi_HR2p5))] # a list
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_bad_multi_HR2p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_bad_multi_HR2p5 <- cbind(candidate_bad_multi_HR2p5, list_KM_sigBin)
# ex. "DUSP5" has bad impact on M stage

# => do NOT need to see "remark=1"
geneid_bad_unimulti_HR2p5 <- isect_HR2p5[[3]] #$`A:B`
candidate_bad_unimulti_HR2p5 <- cbind(data.frame(gene_id=geneid_bad_unimulti_HR2p5), HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$gene_id %in% geneid_bad_unimulti_HR2p5), 2:7])


# export 3 tables
xlsx.addLineBreak(sheet, 10) # 
# header of candidate_bad_unimulti_HR2p5
xlsx.addHeader(wb, sheet, value=paste("Table 2A. The ", nrow(candidate_bad_unimulti_HR2p5), " candiate genes (uni_ & multi_CoxHR > ", bad_FC,") in ", TCGA_cohort, " (bad guy)", sep=""),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_bad_unimulti_HR2p5[,2] <- formatC(candidate_bad_unimulti_HR2p5[,2], format = "e", digits = 2) #as.numeric(as.character)
candidate_bad_unimulti_HR2p5[,3] <- formatC(candidate_bad_unimulti_HR2p5[,3], format = "f", digits = 2) #as.numeric
xlsx.addTable(wb, sheet, data = candidate_bad_unimulti_HR2p5[, c(1:7)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_bad_uni_HR2p5
xlsx.addHeader(wb, sheet, value=paste("Table 2B. The candiate genes (uni_CoxHR >", bad_FC,") in ", TCGA_cohort, " (bad guy)", sep=""),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_bad_uni_HR2p5[,2] <- formatC(candidate_bad_uni_HR2p5[,2], format = "e", digits = 2) #as.numeric(as.character
candidate_bad_uni_HR2p5[,3] <- formatC(candidate_bad_uni_HR2p5[,3], format = "f", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_bad_uni_HR2p5[, c(1:7, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_bad_multi_HR2p5
xlsx.addHeader(wb, sheet, value=paste("Table 2C. The candiate genes (multi_CoxHR >", bad_FC,") in ", TCGA_cohort, " (bad guy)", sep=""),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_bad_multi_HR2p5[,2] <- formatC(candidate_bad_multi_HR2p5[,2], format = "e", digits = 2) #as.numeric(as.character
candidate_bad_multi_HR2p5[,3] <- formatC(candidate_bad_multi_HR2p5[,3], format = "f", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_bad_multi_HR2p5[, c(1:7, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines





#x Export Venn2 by plotFunction() ####
xlsx.addLineBreak(sheet, 1)

# Define function
xlsx.addPlot.venn<-function( wb, sheet, startRow=NULL,startCol=2,
                             width=480, height=480,... )
{ # plot of venn digram"
  
  
  # png(filename = "plot2.png", width = width, height = height,...)
  # #{
  # venn(venn_HR0p5, snames=names_HR0p5,
  #      ilabels = T, counts = T, ellipse = FALSE, zcolor = "forestgreen, darkolivegreen1", opacity = 0.6, size = 15, cexil = 3, cexsn = 0.85, borders = TRUE)
  # #      predefined colors if "style"
  # title <- paste(c("HNSCC survival analysis", "KM P-Value <= ", signif(Bonferroni_cutoff, 3)), sep = "", collapse = "\n")
  # text(500,900, labels = title, cex = 0.85) # (0,0) on bottom_left corner
  # 
  # #}
  # dev.off()
  
  # append Rplot09_venn_HR0.5.tiff
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  }
  system("convert -resize 2048X1536 Rplot09_venn_HR0.5.tiff Rplot09_venn_HR0.5.png", intern = FALSE)
  addPicture("Rplot09_venn_HR0.5.png", sheet=sheet,  startRow = startRow , startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
  #res<-file.remove("plot2.png")
} # end of Define function

# calling
#detach(package:gplots)
# converting .tiff to .png for export
library(venn)
xlsx.addPlot.venn(wb, sheet)
#print(paste("The Venn diagram2 and the candidate genes (good guy): ", filenamex, " were exported successfully."))


#x (good guy table) prognostic features of those genes on lists in TCGA HNSCC cohort ####
# => see "remark=1" 
geneid_good_uni_HR0p5 <- c(isect_HR0p5[[3]], isect_HR0p5[[1]]) #$`A:B`
candidate_good_uni_HR0p5 <- cbind(data.frame(gene_id=geneid_good_uni_HR0p5), HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$gene_id %in% geneid_good_uni_HR0p5), 2:7])
# number, gene_id, p_value, uni_HR, uni_P_value, multi_HR, multi_P_value])
# the list (order) of gene_id has error ???
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_good_uni_HR0p5))] # a list
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_good_uni_HR0p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_good_uni_HR0p5 <- cbind(candidate_good_uni_HR0p5, list_KM_sigBin)
# ex."TXNDC11" higher expression has better prognosis impact as well as smaller T and less N and lower stage

# => see "remark=1" 
geneid_good_multi_HR0p5 <- c(isect_HR0p5[[3]], isect_HR0p5[[2]])
candidate_good_multi_HR0p5 <- cbind(data.frame(gene_id=geneid_good_multi_HR0p5), HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$gene_id %in% geneid_good_multi_HR0p5), 2:7])
# as.binary
list_KM_sig <- candidate_cox[as.numeric(rownames(candidate_good_multi_HR0p5))] # a list
list_KM_sigBin <- data.frame(matrix(NA, nrow = length(geneid_good_multi_HR0p5), ncol = nrow(candidate_cox[[1]])-1))
colnames(list_KM_sigBin) <- c(rownames(candidate_cox[[1]])[-8])
for (ikm in 1:length(list_KM_sig)) 
{
  #list_KM_sigBin[ikm, 2:nrow(candidate_cox[[1]])] <- 
  list_KM_sigBin[ikm,] <- as.numeric(list_KM_sig[[ikm]]$KM_sig[-8]) # or as .factor
}
candidate_good_multi_HR0p5 <- cbind(candidate_good_multi_HR0p5, list_KM_sigBin)
# ex."PLD3" has less N


# => do NOT need to see "remark=1"
geneid_good_unimulti_HR0p5 <- isect_HR0p5[[3]] #$`A:B`
candidate_good_unimulti_HR0p5 <- cbind(data.frame(gene_id=geneid_good_unimulti_HR0p5), HNSCC_OS_marginS_THREE_pvalue005_noCancerGene[(HNSCC_OS_marginS_THREE_pvalue005_noCancerGene$gene_id %in% geneid_good_unimulti_HR0p5), 2:7])
# *** using table 2 to interpretate the impact on survival by DEG of this candidate gene list 
# c("CRHR2", "MYLIP", "NUP210L", "ZKSCAN4")
# => see "remark=1" on HNSCC_survivalAnalysis_marginS_CRHR2.xlsx => higher CRHR2 expression is associated with less LN metastaese
# => see "remark=1" on HNSCC_survivalAnalysis_marginS_CRHR2.xlsx => higher MYLIP expression is associated with smaller tumor size (T)
# => see "remark=1" on HNSCC_survivalAnalysis_marginS_CRHR2.xlsx => higher NUP210L expression is associated with less LN metastaese
# => see "remark=1" on HNSCC_survivalAnalysis_marginS_CRHR2.xlsx => higher ZKSCAN4 expression is associated with smaller tumor size (T)

# export tables
xlsx.addLineBreak(sheet, 10) # 
# header of candidate_good_unimulti_HR0p5
xlsx.addHeader(wb, sheet, value=paste("Table 3A. The ", nrow(candidate_good_unimulti_HR0p5), " consensus candiate genes (uni_ & multi_CoxHR < ", good_FC, ") in ", TCGA_cohort, " (good guy)", sep=""),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_good_unimulti_HR0p5[,2] <- formatC(candidate_good_unimulti_HR0p5[,2], format = "e", digits = 2) #as.numeric(as.character
candidate_good_unimulti_HR0p5[,3] <- formatC(candidate_good_unimulti_HR0p5[,3], format = "f", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_good_unimulti_HR0p5[, c(1:7)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_good_uni_HR0p5
xlsx.addHeader(wb, sheet, value=paste("Table 3B. The candiate genes (uni_CoxHR < ", good_FC, ") in ", TCGA_cohort, " (good guy)", sep=""),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_good_uni_HR0p5[,2] <- formatC(candidate_good_uni_HR0p5[,2], format = "e", digits = 2) #as.numeric(as.character
candidate_good_uni_HR0p5[,3] <- formatC(candidate_good_uni_HR0p5[,3], format = "f", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_good_uni_HR0p5[, c(1:7, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines

xlsx.addLineBreak(sheet, 3)
# header of candidate_good_multi_HR0p5
xlsx.addHeader(wb, sheet, value=paste("Table 3C. The candiate genes (multi_CoxHR < ", good_FC, ") in ", TCGA_cohort, " (good guy)", sep=""),
               level=5, color="black", underline=0)
xlsx.addLineBreak(sheet, 1) # add one blank line
candidate_good_multi_HR0p5[,2] <- formatC(candidate_good_multi_HR0p5[,2], format = "e", digits = 2) #as.numeric(as.character
candidate_good_multi_HR0p5[,3] <- formatC(candidate_good_multi_HR0p5[,3], format = "f", digits = 2) #as.numeric(as.character
xlsx.addTable(wb, sheet, data = candidate_good_multi_HR0p5[, c(1:7, 9:13)], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = F)
xlsx.addLineBreak(sheet, 2)  # add two blank lines



#-------------#
# Write wb to disk ####
##..[2019/07/22] _marginS_ done
# Finalize the workbook to an Excel file and write the file to disk.
#setwd(path_ZSWIM2) # at /marginS/
xlsx::saveWorkbook(wb, filenamex) # file name only, no path
#setwd(path_cohort)
#xlsx.openFile(filenamex) # open file to review

# x the END of R2Excel (with z-cut) ###





#{>*** R2Excel (2020) ####
#2020/03/08 Create an Excel workbook .xlsx file formats can be used.
# 2020/03/12 updated 100cut FDR and wholegenome Bonferroni correction
## export to HNSCC_OS_marginS_candidates_Bonferroni.xlsx
library("r2excel")
filenamex <- paste(TCGA_cohort, "_OS", marginTag, "candidates_Bonferroni", ".xlsx", sep = "") 
# no more "HNSCC_OS_marginS_candidates_Venn.xlsx"
wb <- createWorkbook(type="xlsx")

# Create a sheet in that workbook
sheet <- xlsx::createSheet(wb, sheetName = paste("Survival_candidates"))
# [add data row by row, start from column 2]
#+++++++++++++++++++++++++++++++
## Add paragraph : Author
library("tis") # by Brian Salzer
# today(), arg must be ti, tis, ts, tif, or tifName
select_title <- function(x) {
  switch(x,
         "_marginS_" = "Cohort with surgical margins status (0 or 1)",
         "_marginFree_" = "Cohort with surgical margins free from tumor (0)",
         "_marginPlus_" = "Cohort with surgical margins involving tumor (1)",
         stop("Unknown input")
  )
}
title_candidates_Venn_xlsx <- select_title(marginTag)
author <- paste("Reported by Tex Li-Hsing Chi. \n",
                "tex@gate.sinica.edu.tw \n", title_candidates_Venn_xlsx, "\n", Sys.Date(), sep="")
cat(author)
xlsx.addParagraph(wb, sheet, value=author, isItalic=TRUE, colSpan=5, 
                  rowSpan=4, fontColor="darkgray", fontSize=24)
xlsx.addLineBreak(sheet, 3)


# # Export Cox uni/multi HR vs P-value plot ###
# Define function
xlsx.addPlot.candidates<-function( wb, sheet, file="plot.png", startRow=NULL,startCol=2,
                                   width=480, height=480,... )
{ # plot of Cox uni/multi HR vs P-value plot
  #Append plot to the sheet
  if(is.null(startRow)){
    rows<- getRows(sheet) #list of row object
    startRow=length(rows)+1
  } 
  # append to the file created previously, add 
  # {
  # Cox's Harzard Ratios (univariate) plot [Rplot07_cox_uniHR.tiff]
  # Cox's Harzard Ratios (multivariate) plot [Rplot07_cox_multiHR.tiff]
  # }
  addPicture(file, sheet=sheet,  startRow = startRow, startColumn = startCol) 
  xlsx.addLineBreak(sheet, round(width/20)+1)
} # Define function
# calling
# converting .tiff to .png for export (from Version: ImageMagick 6.9.7-4)
# Cox's Harzard Ratios (multivariate) plot [Rplot07_cox_multiHR.tiff] => cox_multi_Bonf_KMPvalue.tiff
# cox_multi_Bonf_KMPvalue.tiff is Bonferroni KM P-value plot in multivariate analysis
system("convert cox_multi_Bonf_KMPvalue.tiff cox_multi_Bonf_KMPvalue.png", intern = FALSE) # convert is a duplicate
currentRow <- xlsx.addPlot.candidates(wb, sheet, file="cox_multi_Bonf_KMPvalue.png") # startRow + a plot


# header
# a candidate genes table
# HNSCC_OS_marginS_uni_CoxHR2p5_Bonf; bad guy uni/multi 相同
# HNSCC_OS_marginS_multi_CoxHR2p5_Bonf
# [HNSCC_OS_marginS_uni_CoxHR0p5_Bonf; good guy 取交集]
# [HNSCC_OS_marginS_multi_CoxHR0p5_Bonf]
candidates_bad_guy <- HNSCC_OS_marginS_multi_CoxHR2p5_Bonf # after Bonferroni correction
# no more Z-score:a normalized number (frequency) of OS P-values, which KM P-value < 0.05, in each gene
xlsx.addHeader(wb, sheet, value=paste("Table 1. The ", nrow(candidates_bad_guy), " candiate genes overexpressed with poor prognosis in ", TCGA_cohort, "",
                                      " (ranked by Bonferroni corrected Kaplan Meier P-value) ", "\n" , sep=""),
               level=5, color="black", underline=0) # nrow(HNSCC_OS_marginS_pvalue005_zcut), n=6601
#xlsx.addHeader(wb, sheet, value=paste("Cutoff at ", round(cutoff1, 3), " (", percent(surv_OS1$n[1]/(surv_OS1$n[1]+surv_OS1$n[2])), ")", sep = ""),
#               level=5, color="red", underline=0) # total n is taken from surv_OS1$n

#xlsx.addLineBreak(sheet, 1) # add one blank line

# > colnames(HNSCC_OS_marginS_pvalue005_zcut)
# [1] "gene_id" "p_value" "z_score" rename to "Gene_id", "P_value", "Z_score" (why?)
colnames(candidates_bad_guy)[1:4] <- c("Gene_id", "P_value", "Bonferroni P_value", "FDR P_value")
# # # Bonferroni Bonferroni_cutoff; [, c(1,2,4)]
attach(candidates_bad_guy) # removal of as.factor (P_value)
# in scientific notation: formatC of [, c(2)]
candidates_bad_guy$P_value <- formatC(as.numeric(as.character(P_value)), format = "e", digits = 2)
candidates_bad_guy$`Bonferroni P_value` <- formatC(as.numeric(as.character(`Bonferroni P_value`)), format = "f", digits = 3)
candidates_bad_guy$`FDR P_value` <- formatC(as.numeric(as.character(`FDR P_value`)), format = "f", digits = 4)
#candidates_bad_guy$Z_score <- signif(Z_score, digits=4)
candidates_bad_guy <- candidates_bad_guy[order(`Bonferroni P_value`), ] #sorting by order(ascending)
detach(candidates_bad_guy)
# export up to top 10 genes [1:30,], don't show FDR P-value
xlsx.addTable(wb, sheet, data = candidates_bad_guy[, -4], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = TRUE)

xlsx.addLineBreak(sheet, 1)
xlsx.addParagraph(wb, sheet, value=paste("Selection criteria:", "\n", "Bonferroni KM P-value < 0.05", "\n", "Cox's univariate & multivariate HR >=", bad_FC), isItalic=FALSE, colSpan=4, 
                  rowSpan=4, fontColor="darkblue", fontSize=14)
xlsx.addLineBreak(sheet, 5)  # add two blank lines


# ...
# HNSCC_OS_marginS_uni_CoxHR0p5_Bonf; good guy 取交集 inner join
# HNSCC_OS_marginS_multi_CoxHR0p5_Bonf
#candidates_good_guy <- merge(HNSCC_OS_marginS_uni_CoxHR0p5_Bonf, HNSCC_OS_marginS_multi_CoxHR0p5_Bonf, by="gene_id", all=FALSE)
candidates_good_guy <- merge(HNSCC_OS_marginS_uni_CoxHR0p5_Bonf, HNSCC_OS_marginS_multi_CoxHR0p5_Bonf) # after Bonferroni correction
# no more Z-score:a normalized number (frequency) of OS P-values, which KM P-value < 0.05, in each gene
xlsx.addHeader(wb, sheet, value=paste("Table 2. The ", nrow(candidates_good_guy), " candiate genes overexpressed with better prognosis in ", TCGA_cohort,
                                      " (ranked by Bonferroni corrected Kaplan Meier P-value) ", "\n" , sep=""),
               level=5, color="black", underline=0) # nrow(HNSCC_OS_marginS_pvalue005_zcut), n=6601
#xlsx.addHeader(wb, sheet, value=paste("Cutoff at ", round(cutoff1, 3), " (", percent(surv_OS1$n[1]/(surv_OS1$n[1]+surv_OS1$n[2])), ")", sep = ""),
#               level=5, color="red", underline=0) # total n is taken from surv_OS1$n

#xlsx.addLineBreak(sheet, 1) # add one blank line

# > colnames(HNSCC_OS_marginS_pvalue005_zcut)
# [1] "gene_id" "p_value" "z_score" rename to "Gene_id", "P_value", "Z_score" (why?)
colnames(candidates_good_guy)[1:4] <- c("Gene_id", "P_value", "Bonferroni P_value", "FDR P_value")
# # # Bonferroni Bonferroni_cutoff; [, c(1,2,4)]
attach(candidates_good_guy) # removal of as.factor (P_value)
# in scientific notation: formatC of [, c(2)]
candidates_good_guy$P_value <- formatC(as.numeric(as.character(P_value)), format = "e", digits = 2)
candidates_good_guy$`Bonferroni P_value` <- formatC(as.numeric(as.character(`Bonferroni P_value`)), format = "f", digits = 3)
candidates_good_guy$`FDR P_value` <- formatC(as.numeric(as.character(`FDR P_value`)), format = "f", digits = 4)
#candidates_good_guy$Z_score <- signif(Z_score, digits=4)
candidates_good_guy <- candidates_good_guy[order(`Bonferroni P_value`), ] #sorting by order(ascending)
detach(candidates_good_guy)
# export up to top 10 genes [1:30,], don't show FDR P-value
xlsx.addTable(wb, sheet, data = candidates_good_guy[, -4], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = TRUE)

xlsx.addLineBreak(sheet, 1)
xlsx.addParagraph(wb, sheet, value=paste("Selection criteria:", "\n", "Bonferroni KM P-value < 0.05", "\n", "Cox's univariate & multivariate HR <=", good_FC), isItalic=FALSE, colSpan=4, 
                  rowSpan=4, fontColor="darkblue", fontSize=14)
xlsx.addLineBreak(sheet, 5)  # add two blank lines


### 最棒 P-value plot of KM survival analyses (PvalueplotKM_20genes tiff or png 600*500 ** with Bonferroni correction 決定) ####
### 這張圖最棒 Figure 3
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
plot(HNSCC_OS_marginS_pvalue_sorted$p_value, 
     HNSCC_OS_marginS_pvalue_sorted$number, type="p", xaxt="none", yaxt="none",
     axes=FALSE, xlim=c(1e-8, 0.05),
     ylab="Frequency", xlab="Adjusted KM P-value", log="x", cex=0.3) # log scale x or y; main="P-value plot of KM survival analyses", 
# summary(complete.cases(HNSCC_OS_marginS_pvalue_sorted)) == 6429
mtext(side=3, line=0.7, "P-value plot of KM survival analyses", font=3, cex=1.6)
# X-axis 為 0.05 # https://rpubs.com/riazakhan94/297778
# x0 3e-7, 8e-8
# *** no more axis:
#axis(1, seq(0, 0.05, 7e-6), font=2) # X axis, max(HNSCC_OS_marginS_pvalue_sorted$p_value) # 5e-6
#axis(2, seq(0, 190, 50), font=2, las=2, pos=-7e-6) # Y axis, max(HNSCC_OS_marginS_pvalue_sorted$number)
#***p value plot abline 改為 segments(x0, y0, x1, y1, ....)

#abline(h=150, lty=2, col="blue")
#abline(v=Bonferroni_cutoff, lty=2, col="red") # 5.31011e-06
segments(x0=8e-8, y0=100, x1=4e-2, y1=100, lty=2, col="blue") # h line
segments(x0=alpha_HNSCC/nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), y0=0, 
         x1=alpha_HNSCC/nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), y1=170, 
         lty=2, col="red") # Bonferroni_cutoff= 7.8e-06 or 選 0.05/6624; v line
legend("topright", legend=c(paste("Frequency at 100"), paste("Cutoff at 0.05")), lty=2:2, col=c("blue","red"), cex=1.1) # box and font size
       #                                                      signif(alpha_HNSCC/nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), 2), 
       #                                                      "*"))
       
#mtext(side=1, line=5, paste("* P-value cutoff by Bonferroni method (", alpha_HNSCC, "/", nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), 
#                            "=", signif(alpha_HNSCC/nrow(HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR), 2), ")"), font=2, cex=1, adj=0) # (1=bottom, 2=left, 3=top, 4=right)
# n=6429
#legend("topright", legend=c(paste("Frequency at 100"), paste("Cutoff at ", signif(Bonferroni_cutoff, 2), "*")), lty=2:2, col=c("blue","red"), cex=1.1) # box and font size
#mtext(side=1, line=5, paste("* P-value cutoff by Bonferroni method (", alpha_HNSCC, "/", (LUAD_n * n_percent_Bonferroni), "=", signif(Bonferroni_cutoff, 2), ")"), font=2, cex=1, adj=0) # (1=bottom, 2=left, 3=top, 4=right)

# points 10 candidate genes on this P-value plot, 標出 10 bad and 10 good 
# **on this P-value plot (PvalueplotKM_10genes.tiff) => X-axis cutoff 要改為 0.05
# #https://bookdown.org/ndphillips/YaRrr/low-level-plotting-functions.html
#load(file=file.path(path_ZSWIM3, paste(TCGA_cohort, "_OS", marginTag, "multi_CoxHR2p5_Bonf.Rda", sep=""))) # as HNSCC_OS_marginS_multi_CoxHR2p5_Bonf (final candidates)
path_ZSWIM3 <- "marginS"
load(file=file.path(path_ZSWIM3, paste(TCGA_cohort, "_OS", marginTag, "pvalue005KM_sorted_pvalueCox_HR.Rda", sep=""))) # as HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR, with "number"
# good guy candidates_good_guy, n=10
points_x <- candidates_good_guy$P_value # using uncorrected P-value
points_y <- HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[
  HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id %in% candidates_good_guy$Gene_id, 1] # number 
points(x = points_x, y = points_y,
       pch = 16,
       col = "green" ) # transparent("coral2", trans.val = .8)
# text(x = points_x[c(1,3,5,7,9)], y = points_y[c(1,3,5,7,9)],
#      labels = candidates_good_guy$Gene_id[c(1,3,5,7,9)],
#      cex = 0.5, adj = 0,
#      pos = 1)            # Put labels below the points
# text(x = points_x[c(2, 4, 6, 8)], y = points_y[c(2, 4, 6, 8)],
#      labels = candidates_good_guy$Gene_id[c(2, 4, 6, 8)],
#      cex = 0.5, adj = 0,
#      pos = 3)            # Put labels below the points

# bad guy candidates_bad_guy, n=10
points_x <- candidates_bad_guy$P_value # using uncorrected P-value
points_y <- HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR[
  HNSCC_OS_marginS_pvalue005KM_sorted_pvalueCox_HR$gene_id %in% candidates_bad_guy$Gene_id, 1] # number 
points(x = points_x, y = points_y,
       pch = 16,
       col = "red" ) # transparent("coral2", trans.val = .8)
text(x = points_x[c(1,3,5,7,9)], y = points_y[c(1,3,5,7,9)],
     labels = candidates_bad_guy$Gene_id[c(1,3,5,7,9)],
     cex = 0.3, adj = 0,
     pos = 1)            # Put labels below the points
text(x = points_x[c(2, 4, 6, 8, 10)], y = points_y[c(2, 4, 6, 8, 10)],
     labels = candidates_bad_guy$Gene_id[c(2, 4, 6, 8, 10)],
     cex = 0.3, adj = 0,
     pos = 3)            # Put labels above the points

#detach(HNSCC_OS_marginS_pvalue_sorted_noNA_bonferroni)
# save as PDF instead of .tiff for best resolution
# export to .xlsx
system("convert PvalueplotKM_20genes.tiff PvalueplotKM_20genes.png", intern = FALSE) # convert is a duplicate
currentRow <- xlsx.addPlot.candidates(wb, sheet, file="PvalueplotKM_20genes.png") # startRow + a plot


## Write wb to disk ###
##..[2019/07/22] _marginS_ done
# Finalize the workbook to an Excel file and write the file to disk.
xlsx::saveWorkbook(wb, filenamex) # file name only, no path
#xlsx.openFile(filenamex) # open file to review
#system(paste("xdg-open ", filenamex), intern=FALSE) # not in non-graphic Ubuntu console
# the END of R2Excel (2020)###
#}


#{>>** R2Excel (2020/09) ####
#2020/09/08 Create an Excel workbook/ export to HNSCC_OS_marginS_candidates_Bonferroni_Table3.xlsx
library("r2excel")
marginTag <- "_marginS_"
filenamex <- paste(TCGA_cohort, "_OS", marginTag, "candidates_Bonferroni_Table3", ".xlsx", sep = "") 
# no more "HNSCC_OS_marginS_candidates_Venn.xlsx"
wb <- createWorkbook(type="xlsx")

# Create a sheet in that workbook
sheet <- xlsx::createSheet(wb, sheetName = paste("Survival_candidates"))
# [add data row by row, start from column 2]
# > colnames(HNSCC_OS_marginS_pvalue005_zcut)
# [1] "gene_id" "p_value" "z_score" rename to "Gene_id", "P_value", "Z_score" (why?)
# Gene ID, Gene Description, P-value, Adjusted P-value, HR*, 95% CI, HR*, 95% CI, Remark
candidates_good_guy_table3 <- data.frame(matrix(ncol = 9, nrow = 10))
colnames(candidates_good_guy_table3) <- c("Gene ID", "Gene Description", "P-value", "Adjusted P-value", "HR*", "95% CI", "HR*", "95% CI", "Remark")
candidates_good_guy_table3[,c(1,3,4,5,7)] <- candidates_good_guy[,c(1,2,3,5,7)]
#colnames(candidates_good_guy)[1:4] <- c("Gene_id", "P_value", "Bonferroni P_value", "FDR P_value")
# # # Bonferroni Bonferroni_cutoff; [, c(1,2,4)]
attach(candidates_good_guy_table3) # removal of as.factor (P_value)
# in scientific notation: formatC of [, c(2)]
candidates_good_guy_table3$`P-value` <- formatC(as.numeric(as.character(`P-value`)), format = "e", digits = 2)
candidates_good_guy_table3$`Adjusted P-value` <- formatC(as.numeric(as.character(`Adjusted P-value`)), format = "f", digits = 3)
#candidates_good_guy$`FDR P_value` <- formatC(as.numeric(as.character(`FDR P_value`)), format = "f", digits = 4)
candidates_good_guy_table3 <- candidates_good_guy_table3[order(`Adjusted P-value`), ] #sorting by order(ascending)
detach(candidates_good_guy_table3)
# get HR from HNSCC_survivalAnalysis_marginS_ZNF557.Rda
# #$ find . -maxdepth 1 -name "HNSCC_survivalAnalysis_marginS_*.Rda" | xargs [command]
i_geneID <- candidates_good_guy_table3$`Gene ID`
untar_Rda <- data.frame(matrix(ncol = 1, nrow = length(i_geneID)))
for (iig in 1:length(i_geneID)) {
  untar_Rda[iig, 1] <- paste("HNSCC_survivalAnalysis_marginS_", i_geneID[iig], ".Rda", sep="")
  load(file=untar_Rda[iig, 1]) 
  # retrieving 95% CI from tableOS1
  ci95 <- tableOS1[18, c(3,4,7,8)] # `RNAseq(z-score)`
  uni_ci95 <- paste(ci95[1], "-", ci95[2], sep="")
  multi_ci95 <- paste(ci95[3], "-", ci95[4], sep="")
  candidates_good_guy_table3[iig, c(6,8)] <- c(uni_ci95, multi_ci95)
} # end of for loop
# get_Rda_pvalue <- function(geneName) {
#   load(file=paste("HNSCC_survivalAnalysis_marginS_", geneName, ".Rda", sep=""))
#   # load a list = c("tableChi1", "tableOS1", "tableRFS1", "OS_pvalue", "RFS_pvalue") in this .Rda
#   # a example: geneName <- "ZNF557"
#   #OS_pvalue$p_OS[which.min(OS_pvalue$p_OS)]
#   
#   if (nrow(OS_pvalue) > 0)  { # we hit this gene with P-value < 0.05 in KM plot
#     return(min(OS_pvalue$p_OS))} else {return(NA)}
# }


# export up to top 10 genes [1:30,], don't show FDR P-value
xlsx.addTable(wb, sheet, data = candidates_good_guy_table3[, -4], startCol=2,
              fontColor="darkblue", fontSize=12,
              rowFill=c("white", "white"), row.names = TRUE)

xlsx.addLineBreak(sheet, 1)
xlsx.addParagraph(wb, sheet, value=paste("Selection criteria:", "\n", "Bonferroni KM P-value < 0.05", "\n", "Cox's univariate & multivariate HR <=", good_FC), isItalic=FALSE, colSpan=4, 
                  rowSpan=4, fontColor="darkblue", fontSize=14)
xlsx.addLineBreak(sheet, 5)  # add two blank lines
## Write wb to disk ###
##..[2019/07/22] _marginS_ done
# Finalize the workbook to an Excel file and write the file to disk.
xlsx::saveWorkbook(wb, filenamex) # file name only, no path




# human curation [2020/03/08] ####
candidate_bad_unimulti_HR2p5 <- candidates_bad_guy # from Bonferroni correction
save(candidate_bad_unimulti_HR2p5, file=file.path(path_cohort, paste(sub("_", "", marginTag), "candidate_bad_unimulti_HR2p5", ".Rda", sep="")))
# as marginS_candidate_bad_unimulti_HR2p5.Rda
candidate_good_unimulti_HR0p5 <- candidates_good_guy # from Bonferroni correction
save(candidate_good_unimulti_HR0p5, file=file.path(path_cohort, paste(sub("_", "", marginTag), "candidate_good_unimulti_HR0p5", ".Rda", sep="")))
# as marginS_candidate_good_unimulti_HR0p5.Rda
# 
# retrieve .xlsx of each candidate
#{
# retrieve .xlsx from HNSCC.survival.marginS.20500.tar.gz
# survival.txt 20500 filename of .xlsx
#geneNameX <- as.vector(HNSCC_OS_marginS_multi_CoxHR2p5_Bonf$gene_id)
#}
# tar of all HNSCC_survivalAnalysis_marginPlus_*.xlsx of candidate list
# e.x. HNSCC_survivalAnalysis_marginS_ZZZ3.xlsx
# => marginS_candidate_xlsx.tar.gz
geneNameX <- c(as.vector(candidate_bad_unimulti_HR2p5[, c(1)]), as.vector(candidate_good_unimulti_HR0p5[, c(1)]))
# R4> geneNameX
# [1] "DKK1"      "CAMK2N1"   "STC2"      "PGK1"      "SURF4"    
# [6] "USP10"     "NDFIP1"    "FOXA2"     "STIP1"     "DKC1"     
# [11] "ZNF557"    "ZNF266"    "IL19"      "MYO1H"     "FCGBP"    
# [16] "LOC148709" "EVPLL"     "PNMA5"     "KIAA1683"  "NPB" ...bad + good guys

# >tar
#xlsx_list <- as.data.frame(paste(TCGA_cohort, "_survivalAnalysis", marginTag, geneNameX, ".xlsx", sep=""))
xlsx_list <- paste(TCGA_cohort, "_survivalAnalysis", marginTag, geneNameX, ".xlsx", sep="")
#xlsx_list <- paste("./", gsub("_", "", marginTag), "/", TCGA_cohort, "_survivalAnalysis", marginTag, geneNameX, ".xlsx", sep="")
#colnames(xlsx_list) <- "value" # rename V1 to value; for SparkR
#tmp <- tempfile() # name as "/tmp/RtmpRXBLFL/file5fb02526933a"
#writeLines(xlsx_list, con = tmp)
write.table(xlsx_list, file="xlsx_list.txt", sep = "",
          row.names = F, col.names = F, quote = FALSE)
#write.text(xlsx_list, file.path(path_cohort))
# or system(paste("ls", xlsx_list, "> xlsx_list.txt"))
system(paste("tar -czvf ", sub("_", "", marginTag), "candidate_xlsx.tar.gz -T xlsx_list.txt", sep=""))
# HNSCC.survival.marginS.20500.tar.gz
# "tar -czvf marginS_candidate_xlsx.tar.gz -T xlsx_list.txt"
system(paste("tar -xzvf HNSCC.survival.marginS.20500.tar.gz -T xlsx_list.txt"))
>apple$ tar -xzvf /Users/apple/Downloads/Rstudio_marginS_result/HNSCC.survival.marginS.20500.tar.gz  -T xlsx_list.txt
# how to delete them
#$ find . -maxdepth 1 -name "HNSCC_survivalAnalysis_marginS_*.Rda" | xargs rm



# human reviewing.... good and bad guys: its KM plot and Cox HR[2019/07/23]
# > Z-score > 0.8 should be applied at the end ####
# ranking tables by Z-score
zcut <- zcut # 0.8
# mean=0, Max.: 2.3034
load(file=file.path(path_cohort, paste(sub("_", "", marginTag), "candidate_bad_unimulti_HR2p5", ".Rda", sep="")))
load(file=file.path(path_cohort, paste(sub("_", "", marginTag), "candidate_good_unimulti_HR0p5", ".Rda", sep="")))
# bad guy
View(candidate_bad_unimulti_HR2p5)
# top 1 => Z-score 2.3
candidate_bad_unimulti_HR2p5[candidate_bad_unimulti_HR2p5$Gene_id=="DKK1", ]
#x       gene_id  p_value z_score uni_HR uni_P_value multi_HR  multi_P_value
#x 13084    PCTP 7.66e-05    2.30  4.286       0.041    4.915  0.027
# Gene_id  P_value Bonferroni P_value FDR P_value uni_HR
# 4952    DKK1 8.90e-08              0.001      0.0003  2.266
# uni_P_value multi_HR multi_P_value
# 4952       0.001    2.135             0
# Dickkopf WNT Signaling Pathway Inhibitor 1
# 
# good guy
View(candidate_good_unimulti_HR0p5)
# top 1 => Z-score 2.3
candidate_good_unimulti_HR0p5[candidate_good_unimulti_HR0p5$Gene_id=="ZNF557", ]
# Gene_id  P_value Bonferroni P_value FDR P_value uni_HR
# 10  ZNF557 8.64e-08              0.001      0.0003  0.465
# uni_P_value multi_HR multi_P_value
# 10       0.001    0.499             0
# Zinc Finger Protein 557

# > ROC and AUC for validation ####
# insinstall.packages("OptimalCutpoints")
library(OptimalCutpoints)
data(elas) # The elas data set was obtained from the Cardiology Department at the Galicia General Hospital (Santiago de Compostela, Spain). This study was conducted to assess the clinical usefulness of leukocyte elastase determination in the diagnosis of coronary artery disease (CAD).
#########################################################9##
# Youden Index Method ("Youden"): Covariate gender ###
# ###################################################9### 
optimal.cutpoint.Youden<-optimal.cutpoints(X = "elas", status = "status", 
                                           tag.healthy = 0, methods = "Youden", 
                                           data = elas, pop.prev = NULL, 
                                           categorical.cov = "gender", control = control.cutpoints(), 
                                           ci.fit = TRUE, conf.level = 0.95, trace = FALSE) 
summary(optimal.cutpoint.Youden)
# Change the method for computing the confidence interval
# of Sensitivity and Specificity measures, by control = ...
optimal.cutpoint.Youden<-optimal.cutpoints(X = "elas", status = "status", tag.healthy = 0, methods = "Youden", data = elas, pop.prev = NULL, categorical.cov = "gender",
                                           control = control.cutpoints(ci.SeSp = "AgrestiCoull"), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
summary(optimal.cutpoint.Youden)
# Compute the Generalized Youden Index
optimal.cutpoint.Youden<-optimal.cutpoints(X = "elas", status = "status", tag.healthy = 0, methods = "Youden", data = elas, pop.prev = NULL, categorical.cov = "gender",
                                           control = control.cutpoints(generalized.Youden = TRUE), ci.fit = TRUE, conf.level = 0.95, trace = FALSE)
summary(optimal.cutpoint.Youden)
plot.optimal.cutpoints(optimal.cutpoint.Youden)
summary.optimal.cutpoints(optimal.cutpoint.Youden)

# [DAVID] :-) convert to EntrezID ####
# the database for annotation, visualization and integrated discovery (DAVID), a web-based online bioinformatics resource (http://david.abcc.ncifcrf.gov)
# Gene Functional Classification (DAVID, ***PANTHER.db)
#Provide a rapid means to reduce large lists of genes into functionally related groups of genes to help unravel the biological content captured by high throughput technologies.
# RDAVIDWebService: 
# https://gist.github.com/svigneau/9699239; RDavidQuery.R
# how to query David from R, using the RDAVIDWebService package
library("RDAVIDWebService")
# To register, go to: http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html
# https://david.ncifcrf.gov/webservice/services/DAVIDWebService/authenticate?args0=d622101005@tmu.edu.tw
# => return: true
# 2015 updated: https://support.bioconductor.org/p/70091/
david <- DAVIDWebService$new("d622101005@tmu.edu.tw", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
david$is.connected()
getIdTypes(david) # "ENTREZ_GENE_ID"
# This implementation carries all the limitations of DWS as stated at DAVID’s Web site (http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html): 
# (i) gene or term cluster report will handle up to 3000 genes, 
# (ii) a user or computer can compute up to 200 jobs in a day and 
# (iii) DAVID Team reserves the right to suspend any improper uses of DWS without notice.
# Define foreground and background gene lists.
# foreground genes = a gene list (i.e., genes to be analyzed) 
# background genes (i.e., global gene background, the total genes in the genome, n=20500).
# The foreground list should be contained within the background list.

# >Gene ID Conversion Tool: AnnotationDbi
# convert HGNC gene symbol to Entrez gene ID; https://www.genenames.org/about/guidelines/#!/#tocAnchor-1-7
source("https://bioconductor.org/biocLite.R") # Annotation Database Interface, biomaRt
biocLite("AnnotationDbi") 
# https://www.gungorbudak.com/blog/2018/08/07/convert-gene-symbols-to-entrez-ids-in-r/
biocLite('org.Hs.eg.db')
# ls("package:AnnotationDbi")

library(org.Hs.eg.db)
columns(org.Hs.eg.db) # "SYMBOL" -> "ENTREZID"
# Human genomes include both protein-coding DNA genes (20500, 2%) and noncoding DNA (ncDNA, 1025000? 98%)
e2s_hs_genome_symbol <- toTable(org.Hs.egSYMBOL) # n=60118 whole(?) genome in gene symbol

## *** https://bioconductor.org/packages/release/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
biocLite("hgu95av2.db")
suppressPackageStartupMessages({library(hgu95av2.db)})
# AnnotationDbi::select() # returned 1:1 mapping between keys and columns
columns(hgu95av2.db) # "SYMBOL" -> "ENTREZID"
e2s_hs_genome_symbol <- keys(hgu95av2.db, keytype="SYMBOL") # n=60048

# ** converting and mapIds() # one column; select() # more columns
e2s_hs_genome_entrezIDs <- select(hgu95av2.db, keys=e2s_hs_genome_symbol, 
                                  columns=c('ENTREZID'), keytype="SYMBOL") # n=60048
TCGA_whole_genome_entrezIDs <- select(hgu95av2.db, keys=whole_genome, 
                              columns=c('ENTREZID'), keytype="SYMBOL") # n=20500 -> 20504, with NaNs
# table(is.na.data.frame(TCGA_whole_genome_entrezIDs$ENTREZID))
# FALSE  TRUE 
# 17800  2704
# TCGA_whole_genome_entrezIDs has 17800 genes as myBackgroundGenes
geneNameX_entrezIDs <- select(hgu95av2.db, keys=geneNameX, 
                              columns=c('ENTREZID'), keytype="SYMBOL") # n=90
#= convert2EntrezID(IDs=ensemblIDs, orgAnn="org.Hs.eg.db",
#                             ID_type="entrez_gene_id")
# mapping the missing gene symbol: https://www.genecards.org/
# Previous HGNC Symbol -> Current Gene Symbol -> ENTREZID
updated_geneNameX <- c(
	"ATHL1",	"PGGHG", 80162,
	"CCDC64",	"BICDL1", 92558,
	"CNIH",	"CNIH1", 10175,
	"DPH3B",	"DPH3P1", 100132911,
	"FAM10A4",	"ST13P4", 145165,
	"FAM55C",	"NXPE3", 91775,
	"KAZ",	"KAZN", 23254,
	"NDNL2",	"NSMCE3", 56160
)
# CNIH has CNIH1 CNIH2 CNIH3 CNIH4 family
updated_geneNameX <- data.frame(matrix(data=updated_geneNameX, byrow=T, ncol=3))
colnames(updated_geneNameX) <- c("old_SYMBOL", "SYMBOL", "ENTREZID")

geneNameX_entrezIDs_QC <- geneNameX_entrezIDs
u_index <- which(geneNameX_entrezIDs_QC$SYMBOL %in% updated_geneNameX$old_SYMBOL)
#geneNameX_entrezIDs_QC[u_index, 2] <- 
#  updated_geneNameX[which(geneNameX_entrezIDs_QC$SYMBOL[u_index]==as.character(updated_geneNameX$old_SYMBOL)), 3]
for (u_i in 1:length(u_index)) {
  geneNameX_entrezIDs_QC[which(geneNameX_entrezIDs_QC$SYMBOL==updated_geneNameX[u_i, 1]), 2] <- 
    as.numeric(as.character(updated_geneNameX[u_i, 3]))
}
# which(is.na(geneNameX_entrezIDs_QC$ENTREZID))



# submit to DAVID (prepare gene myForegroundGenes) ####
# myForegroundGenes, n=90, HNSCC risky gene list
myForegroundGenes <- geneNameX_entrezIDs_QC$ENTREZID
FG <- addList(david, myForegroundGenes, idType="ENTREZ_GENE_ID", listName="hazards", listType="Gene") #
#, species="Homo sapiens")
write.table(myForegroundGenes, file="DAVID_hnsccForegroundGenes.txt", sep = "",
            row.names = F, col.names = F, quote = FALSE)
# myBackgroundGenes <- whole_genome; 
# TCGA_whole_genome_entrezIDs has n=17800 genes as myBackgroundGenes
myBackgroundGenes <- TCGA_whole_genome_entrezIDs[complete.cases(TCGA_whole_genome_entrezIDs$ENTREZID), 2]
BG <- addList(david, myBackgroundGenes, idType="ENTREZ_GENE_ID", listName="all", listType="Background")
#, species="Homo sapiens")
write.table(myBackgroundGenes, file="DAVID_hnsccBackgroundGenes.txt", sep = "",
            row.names = F, col.names = F, quote = FALSE)
# error: java.net.SocketTimeoutException: Read timed out
# or Error: !listName %in% switch(listType[1], Gene = getGeneListNames(),  .... is not TRUE
# you cannot submit jobs with more than 3000 genes.
# read the R script: https://rdrr.io/bioc/RDAVIDWebService/src/R/DAVIDWebService-class.R
# 

# *** Gene functional classification (with cluster identification)
# to provide an initial glance of major biological functions associated with gene list
# {\cite{HuangNP2009} Nature Protocols: https://www.nature.com/articles/nprot.2008.211
# DAVID website,  
#   Click on the gene name that leads to individual gene reports for in-depth information about the gene.
#   Click on the red 'T' (functional related term reports) to list associated biology of the gene group.
#   Click on 'RG' (related genes) to list all genes functionally related to the particular gene group.
#   Click on the 'green icon' to invoke 2D (gene-to-term) view.
#   Create a new subgene list for further analysis on a subset of the genes.
# }
#***hand upload https://david.ncifcrf.gov/gene2gene.jsp
#Gene List Report, the related genes list (RG)
#The Kappa Statistic is a chance corrected measure of agreement between two sets of categorized data. Kappa result ranges from 0 to 1. The higher the value of Kappa, the stronger the agreement. If Kappa = 1, then there is perfect agreement. If Kappa = 0, then there is no agreement.

# https://david.ncifcrf.gov/content.jsp?file=FAQs.html#25 explanation in detail
# Specifiy annotation categories.
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# ***Get functional annotation chart as R object.
# https://david.ncifcrf.gov/summary.jsp
# ***Pathway: Proteoglycans in cancer => KEGG: https://david.ncifcrf.gov/kegg.jsp?path=hsa05205$Proteoglycans%20in%20cancer&termId=550028858&source=kegg
# ***kidney disease
FuncAnnotChart <- getFunctionalAnnotationChart(david)
# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "hnscc_DAVID_FuncAnnotChart.tsv")

# ***Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)
# Print functional annotation clustering to file (limited to 3000 genes).
getClusterReportFile(david, "hnscc_DAVID_FuncAnnotClust.tsv")


#GSEA ####
shinyPathview() 
# http://bioconductor.org/packages/devel/bioc/vignettes/GDCRNATools/inst/doc/GDCRNATools.html
# 
#Subramanian, Tamayo, et al. (2005, PNAS 102, 15545-15550) and Mootha, Lindgren, et al. (2003, Nat Genet 34, 267-273).
# gene set enrichment analysis
# GSEA/MSigDB
# The criteria for significant enrichment gene sets in 
#GSEA were: P<0.05 and false discovery rate (FDR) <0.25.
# msigdbr, GSEABase, RGSEA, GAGE
# https://bioconductor.org/help/course-materials/2009/SeattleJan09/GSEA/Category.R
# http://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf
# Pathview works with all types and species of KEGG pathways. We plan to support pathways from Reactome (Croft et al., 2011), NCI Pathway Interaction and other databases based on needs in the future install.packages("Pathview")
# http://bioconductor.org/packages/release/bioc/vignettes/pathview/inst/doc/pathview.pdf
biocLite("pathview")
library(pathview)
# x The BACA package

#> PANTHER.db ####
# *** http://pantherdb.org saved workspace, *** hand analysis with .tiff and pathway pie chart
#https://bioconductor.org/packages/release/data/annotation/vignettes/PANTHER.db/inst/doc/PANTHER.db.pdf
# nature protocol: https://www.nature.com/articles/s41596-019-0128-8
library(PANTHER.db)
browseVignettes("PANTHER.db")
pthOrganisms(PANTHER.db) # HUMAN
keytypes(PANTHER.db) 
# [1] "CLASS_ID"     "COMPONENT_ID" "ENTREZ"       "FAMILY_ID"   
#[5] "GOSLIM_ID"    "PATHWAY_ID"   "SPECIES"      "UNIPROT" 
## use select to extract some data
#keys_panther <- c("E1C9F4","O14618")
# View(keys(PANTHER.db)) # n=90742
keys_panther <- keys(myForegroundGenes, keytype = "ENTREZ")
cols_panther <- c("ENTREZ", "FAMILY_ID","GOSLIM_ID","FAMILY_TERM", "PATHWAY_ID", "CONFIDENCE_CODE")
ktype_panther <- "PATHWAY_ID"
select(PANTHER.db, keys_panther, cols_panther, ktype_panther)



# **** there is PANTHER PIE CHARTs [2019/07/25]
# *** new project: Smoking nicotine pathway in hnscc ####
# clinical feature: tobacco_exposure [high/low]
# => MYO10 in (nicotine) Nicotinic acetylcholine receptor signaling pathway => Myosin
# http://pantherdb.org/genes/gene.do?acc=HUMAN%7CHGNC%3D7593%7CUniProtKB%3DQ9HD67
# http://identifiers.org/panther.pathway/P00044 # from paxtoolsr
# https://apps.pathwaycommons.org/pathways?uri=http%3A%2F%2Fidentifiers.org%2Fpanther.pathway%2FP00044
# ***TFs pickup from panther as well

# https://www.pathwaycommons.org/pcviz/#pathsbetween/27102,23676,7009,55217,58488,6695,10857,81533,94081,3198,153768,219995,4723,23204,647033,56160,8815,23193,100132911,1717,55299,26509,51751,57222,145165,22800,7511,1852,10783,9922,83752,4494,9167,10175,60625,8841,55839,79366,909,55298,91775,55142,83540,25914,347902,3916,2074,5648,400935,728743,2145,728294,163050,353497,83416,65095,117581,23254,53358,116535,1101,92558,25897,84824,10497,643641,327657,7593,5600,30832,1240,116984,58491,7768,51043,5799,80317,339761,5329,79786,80162,962,924,29964,286133,7280,6909,90007,64976,345275
# CDEs dictionary of GDC: https://docs.gdc.cancer.gov/Data_Dictionary/gdcmvs/
#1) ***checking TCGA cohort: smoker, ex-smoker(quit for 15 years +/-) 重要, never-smoker; 
#   and smoking_history: exposure_intensity, exposure_duration, quit_for_year
#   CDE (Collection)=>
#   2955385 (caDSR): number_pack_years_smoked
#   2181650 (caDSR): 'Lifelong Non-Smoker', 'Current Smoker', 'Current Reformed Smoker for > 15 yrs', 
#     'Current Refomed Smoker for < or = 15 yrs', 'Current Reformed Smoker Duration Not Specified'
#   and smokeless_history: exposure_intensity, exposure_duration, quit_for_year, exposure_regularly
load(file=file.path("~/R", "HNSCC.clinical.Fire.Rda"))
colnames(HNSCC.clinical.Fire) # smoking features
#{# # alcohol c(5, 6, 39) 
  # [1] "alcohol_history_documented"           
  # [2] "amount_of_alcohol_consumption_per_day"
  # [3] "frequency_of_alcohol_consumption
  # => # tobacco/{smoking/smokless} c(66, 83:88, 95, 103)
  # [1] "number_pack_years_smoked"   # lifetime tobacco exposure defined as number of cigarettes smoked per day x number of years smoked divided by 20.
  # [2] "smokeless_tobacco_use_age_at_quit" 
  # [3] "smokeless_tobacco_use_age_at_start"
  # [4] "smokeless_tobacco_use_at_diag"     
  # [5] "smokeless_tobacco_use_per_day"     
  # [6] "smokeless_tobacco_use_regularly"   
  # [7] "stopped_smoking_year"              
  # [8] "tobacco_smoking_history"           #  1   2   3   4   5 
  # which:                                   122 178  73 140   2
  # 'Lifelong Non-Smoker', 1
  # 'Current Smoker', 2
  # 'Current Reformed Smoker for > 15 yrs', 3
  # 'Current Refomed Smoker for < or = 15 yrs', 4
  # 'Current Reformed Smoker Duration Not Specified', 5
  # ==> [tobacco_exposure] a new binomial variable (categorizing) as risk of [high/low]:
  #       level high by c(2, 4)
  #       level low by c(1, 3, 5)
  # 
  # [9] "year_of_tobacco_smoking_onset"
  # # HPV c(43:47); IHC or HPV
  # [1] "hpv_call_1"                "hpv_call_2"               
  # [3] "hpv_status"                "hpv_status_by_ish_testing"
  # [5] "hpv_status_by_p16_testing"
#}


# 2) nicotine gene list from panther (n=62)
# http://amp.pharm.mssm.edu/Harmonizome/gene_set/Nicotinic+acetylcholine+receptor+signaling+pathway/PANTHER+Pathways
# http://www.pantherdb.org/pathway/pathDetail.do?clsAccession=P00044
# 62 proteins participating in the Nicotinic acetylcholine receptor signaling pathway pathway from the PANTHER Pathways dataset.
# brought you by Harmonizome(\cite{})
nicotine_gene62 <- read.csv(stdin(), header=T, sep="\t")
View(nicotine_gene62)

# venn, geneNameX_entrezIDs_QC vs nicotine_gene62
# *** there is no intersection
venn_nicotine <- list(geneNameX_entrezIDs_QC$SYMBOL, nicotine_gene62$Symbol) # Gene Symbol
names_nicotine <- c(paste("the candidate of", TCGA_cohort, "biomarker ( n =", nrow(geneNameX_entrezIDs_QC), ")"), 
                    paste("the nicotine pathway ( n =", nrow(nicotine_gene62), ")"))
library(gplots)
tmp_ni <- venn(venn_nicotine, names=names_nicotine, show.plot=F) #library(gplots); the group count matrix alone
isect_nicotine <- attr(tmp_ni, "intersections")
detach(package:gplots)
library(venn)
tiff("Rplot11_venn_hnscc_nicotine.tiff", units="cm", width=5, height=5, res=300) 
venn(venn_nicotine, snames=names_nicotine,
     ilabels = T, counts = T, ellipse = FALSE, zcolor = "red, deeppink", opacity = 0.6, size = 15, cexil = 0.7, cexsn = 0.3, borders = TRUE)
# meta-language 1 0 or -
#title <- c(paste(TCGA_cohort, "survival analysis"), paste("(KM P-Value <=", signif(alpha_HNSCC, 3), ")")) #, collapse = "\n")
title <- paste(TCGA_cohort, "biomarker is not overlaped with the Nicotinic acetylcholine receptor signaling pathway")
text(500,900, labels = title[1], cex = 0.60) # (0,0) on bottom_left corner
text(500,855, labels = title[2], cex = 0.40) 
# n=0
dev.off()


# Pathway Components "Category ID": n=13, pathway ontology term, P as pathway:
nicotine_Components <- data.frame(matrix(c("P01098",
"P01097",
"P01096",
"P01095",
"P01094",
"P01093",
"P01092",
"P01091",
"P01090",
"P01089",
"P01088",
"P01087",
"P01086"), byrow = F, ncol=1)) # source: panther.db, Component Accession
#nicotine_gene_json <- fromJSON(file.path(path_cohort, "Nicotinic_acetylchol.json"), simplifyDataFrame=TRUE)
#View(nicotine_gene_json[[1]][nicotine_gene13_json[[1]][,2]=="bp:ProteinReference",1:6])
#


# 3) => categorizing the tobacco_exposure feature: smoking ratio (high/low) in cohort
# stored at hnscc_nicotine$tobacco_exposure
load(file="~/R/HNSCC.clinical.RNAseq.Fire.Rda") # n=521, as clean6_oscc
#hnscc_n521 <- clean6_oscc
# included $$.clinico_mRNA.Fire for survival analysis
# 
hnscc_nicotine <- clean6_oscc[ , c(1:9,11)]
# # inner join by merge
hnscc_nicotine <- merge(hnscc_nicotine, 
                        subset(HNSCC.clinical.Fire, select=c("tcga_participant_barcode", "tobacco_smoking_history")))
# in HNSCC.clinical.Fire: n=528
# [1] "tcga_participant_barcode" 
# [103] "tobacco_smoking_history"         #  1   2   3   4   5 
# which:                                   122 178  73 140   2
# 'Lifelong Non-Smoker', 1
# 'Current Smoker', 2
# 'Current Reformed Smoker for > 15 yrs', 3
# 'Current Refomed Smoker for < or = 15 yrs', 4
# 'Current Reformed Smoker Duration Not Specified', 5
# ==> [tobacco_exposure] a new binomial variable (categorizing) as risk of [high/low]:
#       level low by c(1, 3, 5)
#       level high by c(2, 4)
#
# in hnscc_nicotine, n=521
# [11] "tobacco_exposure"             1   2   3   4   5 
#   (as "tobacco_smoking_history)   117 177  73 139   2 
# categorizing
colnames(hnscc_nicotine)[11] <- "tobacco_exposure"
hnscc_nicotine$tobacco_exposure[hnscc_nicotine$tobacco_exposure %in% c(1, 3, 5)] <- 1 
# low as 1, recoding => n=192
hnscc_nicotine$tobacco_exposure[hnscc_nicotine$tobacco_exposure %in% c(2, 4)] <- 2 
# high as 2 => n=316
#smoking % (high/low) in cohort
hnscc_smoking_rate <- table(hnscc_nicotine$tobacco_exposure)[2]/nrow(hnscc_nicotine)
# 60.65% is in high risk smoking group (HNSCC cohort n=521)

# 4) => by its significant cutoff (high/low)
# kappa statistics? or Chi-square test


 
# 5) [hypothesis] to prove by bioinformatic way: 
# "Nicotinic acetylcholine receptor signaling pathway"
# nicotine consumption has oral cancer prognostic impact as well as tumor initiation
# ***article review...Embase or PubMed 
#{
#(no citation) Epithelial cell nicotinic acetylcholine receptor expression in head and neck squamous cell carcinoma pathogenesis
# Carracedo et al., Anticancer Research 2007
#classical nicotine- and alcohol-associated carcinoma.
res_df <- read.csv(file=file.choose(), header=T) # with abstract, PUI as embase
res_df <- rbind(res_df, read.csv(file=file.choose(), header=T)) # 35 articles
# get its PMID for PubCurator
# 
#}

# [R_mbp>] run n=521 cohort, survival analysis without RNAseq, and cutoff at smoking_feature high/low.
#  (update local Rstudio (for macOS) and R 3.6.1 for test analysis)
# $ curl -# --url  "https://download1.rstudio.org/desktop/macos/RStudio-1.2.1335.dmg" -o /Users/apple/RStudio-1.2.1335.dmg
# or for 2009 Macbook Pro: $ curl -# --url  "http://download1.rstudio.org/RStudio-1.1.463.dmg" -o /Users/apple/RStudio-1.1.463.dmg
# x $ brew upgrade r --with-java


# [R4>] Instance-4 at GCE re-run 6 days for cutoff finder with new tobacco_feature [nicotine pathway in hnscc]
# [tobacco_feature] since [2019/07/27] ####
# save(clean6_oscc, file="~/R/HNSCC.clinical.RNAseq.Fire.Rda") # n=521, 85Mb [2019/06/06]
load(file="~/R/HNSCC.clinical.RNAseq.Fire.Rda") # as clean6_oscc
# tobacco_feature should be place after margin
clean6_oscc_tobacco <- clean6_oscc[ , 1:8]
# inner join by merge
clean6_oscc_tobacco <- merge(clean6_oscc_tobacco, subset(hnscc_nicotine, select=c("tcga_participant_barcode", "tobacco_exposure")))
# , subset(HNSCC.clinical.Fire, select=c("tcga_participant_barcode", "tobacco_smoking_history"))
clean6_oscc_tobacco <- merge(clean6_oscc_tobacco, clean6_oscc[ , c(1, 9:ncol(clean6_oscc))])
# tobacco_exposure at column 9; then updated HNSCC.clinical.RNAseq.Fire.Rda
save(clean6_oscc_tobacco, file="~/R/HNSCC.clinical.RNAseq.Fire.Rda") # as clean6_oscc_tobacco
# $ scp  tex@35.201.169.0:~/R/HNSCC.clinical.RNAseq.Fire.Rda ~/hnscc_github
# HNSCC.clinical.RNAseq.Fire.Rda => download to MBP
# done: -rw-r--r-- 1 tex tex  85031069 Jul 26 16:19 HNSCC.clinical.RNAseq.Fire.Rda


# 
# Visualize genes on *** MSKCC: Pathway Commons databases, BioCarta and KEGG pathway maps ####
# http://www.pathwaycommons.org PaxtoolsR: Access Pathways from Multiple Databases through BioPAX and Pathway Commons
# the querying Pathway Commons (PC) molecular interaction database that are hosted by the Computational Biology Center at
# Memorial Sloan-Kettering Cancer Center (MSKCC).
# MSKCC: Pathway Commons databases
# include: BIND, BioGRID, CORUM, CTD, DIP, DrugBank, HPRD, HumanCyc, IntAct, KEGG, MirTarBase, Panther, PhosphoSitePlus, Reactome, RECON, TRANSFAC.
# receipe *** http://bioconductor.org/packages/release/bioc/vignettes/paxtoolsr/inst/doc/using_paxtoolsr.html
# "Using PaxtoolsR: A BioPAX and Pathway Commons Tutorial in R"
# 2 May 2019
# The Biological Pathway Exchange (BioPAX) format is a community-driven standard language to represent biological pathways
biocLite("paxtoolsr")
library(paxtoolsr)
help.search("paxtoolsr")
help(graphPc)
# the OWL (Web Ontology Language) file format; XML object.
searchResults <- searchPc(q = "nicotine", type = "pathway", verbose = TRUE)
# => http://www.pathwaycommons.org/pc2/search.xml?q=nicotine&page=0&type=pathway 
library(plyr)
searchResultsDf <- ldply(xmlToList(searchResults), data.frame)
simplifiedSearchResultsDf <- searchResultsDf[, c("name", "uri", "biopaxClass")]
(simplifiedSearchResultsDf$name[1:6])
# [1] Nicotine_degradation                              
# [2] Nicotine Metabolism                               
# [3] Nicotine Activity on Chromaffin Cells             
# [4] Nicotine pharmacodynamics pathway                 
# [5] Nicotine Activity on Dopaminergic Neurons         
# [6] Nicotinic acetylcholine receptor signaling pathway # http://identifiers.org/panther.pathway/P00044

# *** chapter 7: Common Data Visualization Pathways and Network Analysis
#7.1 Visualizing SIF Interactions from Pathway Commons using R Graph Libraries
library(igraph)
#7.2 Pathway Commons Graph Query
#7.3 Overlaying Experimental Data on Pathway Commons Networks
#7.4 Network Statistics => SIF Network Statistics




# IHC cross validation



# Tobacco HR plot under candidate gene signature ####
# [2019/08/11]
# load(file=file.path(path_ZSWIM2, paste(TCGA_cohort, "_OS", marginTag, "pvalueKM_candidate_cox.Rda", sep="")))
# as candidate_sample, candidate_cox, n_percent_Bonferroni,
# R4> colnames(candidate_sample)
# [1] "number"  "gene_id" "p_value"
# candidate_cox[[]] from aa to bb (full TableChi1, table 3 + KM P-value)
# at HNSCC_OS_marginS_pvalueKM_candidate_cox.Rda
# > Tobacco HR plot under candidate gene signature
# high in bad guy gene vs low in good guy gene


# > percentage of tobacco high on each gene Cox PH table 3
# HNSCC_OS_marginS_THREE_pvalue005_noCancerGene
# 

# > HR of tobacco high (和 margin or clinical T 比較)
# (HR when P-value <= 0.05)


## margin issue ####
## a new comparison table for impact genes ####
paste("HNSCC_survivalAnalysis_marginS_", geneName, ".Rda")
# a new comparison for margin issue
paste("HNSCC_survivalAnalysis_marginFree_", geneName, ".Rda")

##
oby <- read.table(header=TRUE, text='
                  Var1 Var2 Freq
                  0      1    36
                  1      1   91
                  0     2    41
                  1       2   88
                  ')
# x oby <- read.table(stdin(), header=TRUE) 
# 
# 
# 


# >Analysis finish; tar until here (refinement ok) ####
# R4> options(prompt="R4_plus>")
# tar
#{
#好用的工具 guidebook
# ***bash $ TODAY=`date +"%b %d"`;ls -l | grep "$TODAY"
# # list today's files
# https://www.howtogeek.com/248780/how-to-compress-and-extract-files-using-the-tar-command-on-linux/
#   https://www.gnu.org/software/tar/manual/tar.html

# =tar and scp from a list of files .xlsx 
# (download
# generated genes list .xlsx and ZSWIM2_free_archive.Rda; 
# analysis -> .tiff, HNSCC_OS_marginS_candidates_Venn.xlsx, HNSCC_OS_marginS_THREE_pvalue005_noCancerGene.Rda)
# $ ls HNSCC_survivalAnalysis_marginS_*.* > marginS_list.txt
# $ tar -czvf /home/tex/R/marginPlus_xlsx.tar.gz -T marginPlus_list.txt   # or  list as many directories
# error of tar: suggestion that Removing leading `/' from member names (files list)
# $ info tar # tar -t --list  —remove-file...
# 
# $ scp  tex@35.201.169.0:~/margin*_xlsx.tar.gz ./
# tex@instance-4:$ ~/R/LUAD_Peter_survival$ sudo mv LUAD_survivalAnalysis_marginS*.* ./survivalAnalysis_marginS/
# tar -xvf marginPlus_xlsx.tar HNSCC_survivalAnalysis_marginPlus_*.Rda
# $ tar -xzvf archive.tar.gz -C /tmp # to extract them to /tmp
# $ tar -xzf archive.tar.gz --overwrite
# $ ls -halt | grep -v ".Rda"
