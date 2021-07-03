# just for a venn of Figure 1
# TCGA 20500
HNSCC_OS_marginFree_pvalueBonferroni_sorted <- data.frame(matrix(seq(1,20500,1), nrow=20500))
HNSCC_OS_marginFree_pvalueBonferroni_sorted$gene_id <- data.frame(matrix(seq(1,20500, 1), nrow=20500))

# GSE65858 31000
HNSCC_OS_marginS_pvalueBonferroni_sorted <- data.frame(matrix(seq(1,20500,1), nrow=20500))
HNSCC_OS_marginS_pvalueBonferroni_sorted$gene_id <- data.frame(matrix(seq(20500-2,20500-2+20499, 1), nrow=20500))
  
venn_marginSFP <- list(HNSCC_OS_marginFree_pvalueBonferroni_sorted$gene_id, HNSCC_OS_marginS_pvalueBonferroni_sorted$gene_id) #, 
#                       HNSCC_OS_marginPlus_pvalueBonferroni_sorted$gene_id)
names_marginSFP <- c(paste("GSE65858"), paste("TCGA"))#, paste("margin[+]"))


library(venn)
library(ggplot2)
library(ggpolypath)
# https://cran.r-project.org/web/packages/venn/venn.pdf
#tiff("venn_marginSFP.tiff", units="cm", width=5, height=5, res=300) 
# # saving as .tiff (by tiff())
x <- list(GSE65858 = 1:31000, TCGA = seq(31000-2,31000-2+20499, 1))
venn(x, #snames=names_marginSFP,
     ilabels = F, counts = F, 
     plotsize = 15,
     ellipse = FALSE, opacity = 0.6, size = 20, 
     ilcs = 1, sncs = 1.5,
#     cexil = 2, #Character expansion for the intersection labels
#     cexsn = 4, #Character expansion for the set names
     borders = FALSE,
     box=FALSE,
     fontfamily = "sans",
     fontface = "bold",
     ggplot = FALSE,
#     , cex =2,
     zcolor = "deeppink", "lightgreen")#,
#     filename = 'rplot_venn_diagramm.png',
#     output=TRUE)
#      predefined colors if "style"
#http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
# meta-language 1 0 or -
text(centroid[1], centroid[2], labels = "3", cex = 3)
# 
title <- c(paste("HNSCC", "survival analyses: candidate genes"))
           #, paste("(KM P-Value <=", signif(Bonferroni_cutoff, 2), ")"))

text(500,950, labels = title, cex = 1.30) # (0,0) on bottom_left corner
#text(500,100, labels = title[2], cex = 1.20) 
#
library(ggforce)
#ggplot # 0.866
corr_x <- 1.3
df.venn <- data.frame(x = c(corr_x, -corr_x),
                      y = c(-0.5, -0.5),
                      labels = c("TCGA","GSE65858"))
ggplot(df.venn)+
  geom_circle(aes(x0 = x, y0 = y, r = 1.5, fill = labels, color=c("deeppink", "green")),
    alpha = .3, size = 1,
    colour = "grey50",
    show.legend=F) +
  scale_fill_manual(values = c('cornflowerblue', 'chartreuse'))+
  coord_fixed() +
  theme_void()+
annotate(geom="text", x=-corr_x, y=1.2, #1.9, 
         label="GSE65858",
         color="black", size = 5)+
#TCGA
  annotate(geom="text", x=corr_x, y=1.2, #1.9, 
           label="TCGA",
           color="black", size = 5)+
# intersection
annotate(geom="text", x=-corr_x, y=-0.5, #1.9, 
         label="30997",
         color="black", size = 5)+
annotate(geom="text", x=0, y=-0.5, #1.9, 
         label="3",
         color="red", size = 7)+
annotate(geom="text", x=corr_x, y=-0.5, #1.9, 
         label="20497",
         color="black", size = 5)
# exported as Rplot09_venn_GSE65858_TCGA.pdf
ggsave("Rplot09_venn_GSE65858_TCGA.pdf") 


# text(500,950, labels = paste(c("GSE65858", "TCGA")), cex = 1.30)
# c("deeppink", "lightgreen")