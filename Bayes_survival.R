# logistic regression
# binominal or multinomial (softmax function)
# https://www.datacamp.com/community/tutorials/logistic-regression-R?utm_source=adwords_ppc&utm_campaignid=1455363063&utm_adgroupid=65083631748&utm_device=c&utm_keyword=&utm_matchtype=b&utm_network=g&utm_adpostion=&utm_creative=278443377095&utm_targetid=aud-392016246653:dsa-429603003980&utm_loc_interest_ms=&utm_loc_physical_ms=9040379&gclid=CjwKCAjwq7aGBhADEiwA6uGZp1XeeIc8HsOi0QPz_wVNjS3ZRQcrDwOsTsTAC2nxts_27he0JNk_HxoC7WEQAvD_BwE
install.packages("ISLR")


# missing map
library(Amelia)
library(mlbench)
missmap(Smarket, col=c("blue", "red"), legend=FALSE)


# correlation plot
# A dot-representation was used where blue represents positive correlation and red negative. The larger the dot the larger the correlation. 
library(corrplot)
correlations <- cor(Smarket[,1:8])
corrplot(correlations, method="circle")

# There's a pairs() function which plots the variables in Smarket into a scatterplot matrix. In this case, Direction, your binary response, is the color indicator.
pairs(Smarket, col=Smarket$Direction)

# the density distribution of each variable broken down by Direction value. The density plot by Direction can help see the separation of Up and Down.  
library(caret)
x <- Smarket[,1:8]
y <- Smarket[,9]
scales <- list(x=list(relation="free"), y=list(relation="free"))
featurePlot(x=x, y=y, plot="density", scales=scales)
# You can see that the Direction values overlap for all of these variables, meaning that it's hard to predict Up or Down based on just one or two variables.
# => So, we need a logistic regression

# Logistics Regression
# family = binomial 
# it assumes a linear relationship between link function (***maximum likelihood) and independent variables in logit model
# p(X)(1−p(X)) the odds ratio
# logit = log(P(x)/(1-P(x))) = β0+β1X
# l(β0,β1)=p(X)(1−p(X))
glm.fit <- glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume, data = Smarket, family = binomial)

#returns the estimate, standard errors, z-score, and p-values on each of the coefficients
# It also gives you the null deviance (the deviance just for the mean) and the residual deviance (the deviance for the model with all the predictors). 
# There's a very small difference between the 2, along with 6 degrees of freedom.
summary() 

# training ####
# This will make predictions on the training data that you use to fit the model and give me a vector of fitted probabilities (0-1, up or down).
# Creating Training and Test Samples before your prediction
# 80% of the patients are selected at random (without replacement) as training set. The remaining 20% is used as validation set.
#attach(Smarket)
#train = Year<2005   # or by random 80%
#detach(Smarket)

# from HTLR; ####
# # or by random 80%
set.seed(1234)
library(HTLR)
# y (p=2000 genes); x (n=510 participants)
n <- 510
p <- 2000
dat <- gendata_MLR(n = n, p = p, NC=2) # y class as alive or dead

## x
means <- rbind( # 10 * 3
  c(0, 1, 0),
  c(0, 0, 0),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1),
  c(0, 0, 1)
) * 2 # matrix * 2

means <- rbind(means, matrix(0, p - 10, 3))

A <- diag(1, p) # diagonal

A[1:10, 1:3] <-
  rbind(
    c(1, 0, 0),
    c(2, 1, 0),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1),
    c(0, 0, 1)
  )
# muj as means
# SGM dim as p * p
# sd_g Numeric value indicating noise level
dat <- gendata_FAM(n, means, A, sd_g = 0.5, stdx = TRUE)
#

ggplot2::qplot(dat$y, bins = 6)
require(corrplot)
cor(dat$X[ , 1:11]) %>% corrplot::corrplot(tl.pos = "n", method="circle")


# 510 * 80%
# muj and means are removed (why?)
# training set for training
# test set for testing (never seen by model)
# => there is no Cross-validation (CV) by HTLR
dat <- split_data(dat$X, dat$y, n.train = n*0.8) # python style

# *** feature selection => includes gene signature that are related to a grouping structures with similar biological functions tend to have similar expressions. 
# For example, a group of co-regulated genes that are preferable to be taken into consideration in select features. 
# https://www.nature.com/articles/s41598-020-66466-z Scientific Reports
# 
# sparse regularization techniques such as the LASSO
# ## Logistic Regression feature selection by LASSO
# Hastie et al. (2017) concluded that LASSO outperformed two state-of-the-art FS algorithms, forward stepwise regression Weisberg (1980) and best-subset- selection Bertsimas et al. (2016).
# gene expression and survival (time-to-event data), LASSO using Cox proportional hazard regression model was suggested by Park and Hastie (2007)\cite{Tibshirani1997}\cite{Park2007a}. 
# the L_{1}$-regularized path algorithm for the Cox Proportional Hazards Models using the TCGA/GSE65858 survival data
 

## Tibshirani (1996)\cite{Tibshirani1996} suggested, in the statistical literature, Least Absolute Shrinkage and Selection Operation (LASSO); for Cox\cite{Tibshirani1997}
## (or gOMP https://www.biorxiv.org/content/10.1101/431734v1.full.pdf; R package MXM https://cran.r-project.org/web/packages/ MXM/index.html.)
##  the penalized minimization of an objective function (usually log-likelihood)
##  residual-based algorithm: The complexity, or number of operations, required by LASSO is O(np2) Rosset and Zhu (2007), where n and p denote the sample size and number of features, respectively. (e.x. n=270, p=30330);
##  This regularization (LASSO) has been successfully applied to several biomedical problems with p ≫ n and correlated features\cite{Branders2014}.
##  argmax  

## # HTLR: Bayesian Logistic Regression with Heavy-tailed Priors
# we don't find BayesHL (Bayesian Robit regression method with Hyper-LASSO priors), using HTLR instead
# https://cran.r-project.org/web/packages/HTLR/vignettes/simu.html
install.packages("HTLR")
library(HTLR)
library("bayesplot") # This is bayesplot version 1.8.1

##  Compute a Coxlogit model on the train\cite{Branders2014}.
##  fit HTLR ####
fit.t <- htlr(dat$x.tr, dat$y.tr)
fit.t2 <- htlr(X = dat$x.tr, y = dat$y.tr, 
               prior = htlr_prior("t", df = 1, logw = -20, sigmab0 = 1500), 
               iter = 4000, init = "bcbc", keep.warmup.hist = T)

summary(fit.t2, features = c(1:10, 100, 200, 1000, 2000), method = median)

# Plot interval estimates from posterior draws by bayesplot:
post.t <- as.matrix(fit.t2, k = 2)
## signal parameters
mcmc_intervals(post.t, pars = c("Intercept", "V1", "V2", "V3", "V1000"))
# Trace plot of MCMC draws:
as.matrix(fit.t2, k = 2, include.warmup = T) %>%
  mcmc_trace(c("V1", "V1000"), facet_args = list("nrow" = 2), n_warmup = 2000)



##  
glm.fit <- glm(Direction ~ Lag1 + Lag2 + Lag3 + Lag4 + Lag5 + Volume,
               data = Smarket,
               family = binomial,
               subset = train)


# predict with test set ####
# # a matrix two-by-two table
# A glance at the prediction accuracy:
y.class <- predict(fit.t2, dat$x.te, type = "class")
print(paste0("prediction accuracy of model t2  = ", 
             sum(y.class == dat$y.te) / length(y.class)))

predict(fit.t2, dat$x.te, type = "response") %>%
  evaluate_pred(y.true = dat$y.te)
# This function compares the prediction results returned by a classifier with ground truth, and finally gives a summary of the evaluation.
# $amlp
#> [1] 0.6433173 => 0.412

#
glm.probs <- predict(glm.fit,
                     newdata = Smarket[!train,],
                     type = "response")
#glm.probs <- predict(glm.fit,type = "response")
# as a classifier: I'll turn the probabilities into classifications by thresholding at 0.5. 
glm.pred <- ifelse(glm.probs > 0.5, "Up", "Down")


# Predict the class and survival on the test.
# prediction of test set (sensitivity and speicificity)
# a confusion matrix two-by-two table
Direction.2005 = Smarket$Direction[!train]
table(glm.pred, Direction.2005)
# The mean gives a proportion of 0.52 (a classification rate)
mean(glm.pred == Direction.2005)


## or MSE (training by alpha 10/10=1: LASSO)
## alpha 0: ridge regularization
## alpha 0.5: mixed above i.e., elastic net
fit10 <- cv.glmnet(x.train, y.train, type.measure="mse", 
          alpha=i/10,family="gaussian") # with cross-validation (cv.glmnet)
yhat10 <- predict(fit10, s=fit10$lambda.1se, newx=x.test)
mse10 <- mean((y.test - yhat10)^2)
# n <- 1000  # Number of observations
#p <- 5000  # Number of predictors included in model
#real_p <- 15  # Number of true predictors
# then LASSO is good at picking up a small signal through lots of noise. 
## *** regularized Cox model  with CV (x genes do not need a cutoff) 須要保密 *** ####
fit <- glmnet(x, y, family = "cox")
# C: the Harrell concordance index, only considers comparable pairs
set.seed(1)
# "coxnet" objects: cvfit
cv.fit <- cv.glmnet(x, y, family = "cox", nfolds=10, type.measure = "C") # C or deviance
#The survfit method (for plot) is available for cv.glmnet objects as well. By default, the s value chosen is the “lambda.1se” value stored in the CV object.
plot(survival::survfit(cv.fit, x = x, y = y2)) # default s = "lambda.1se"
# plot 看起來只是 KM plot



# 
# upgrade R to v4.1 updateR(admin_password = 'Admin user password')
library(devtools)
install_github('andreacirilloac/updateR')
library(updateR)
updateR()
# done

# 
# #######*** RLassoCox####
# how to select feature? by non-zero coefficients
# no R package LASSO; 
# RLassoCox, A reweighted Lasso-Cox by integrating gene interaction information
# \cite{Simon2011}\cite{Liu2021}
# https://bioconductor.org/packages/devel/bioc/vignettes/RLassoCox/inst/doc/RLassoCox.pdf
#BiocManager::install(version='devel')
#BiocManager::install(version = '3.12')
BiocManager::install("RLassoCox", type = "source", checkBuilt = TRUE)
install_github('weiliu123/RLassoCox')
#install.packages("/Users/texchi/Downloads/RLassoCox_1.1.0.tar.gz", repos = NULL, type="source")
#install.packages("glmnet") # basic lasso
library(RLassoCox)
citation("RLassoCox")
browseVignettes("RLassoCox")

# load TCGA HNSCC dataset
#load(file="HNSCC.clinical.RNAseq_tobacco.Fire.Rda") #clinical data only
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

# convert Gene.symbol (z.score_A1BG) to NCBI Entrez Gene: (hsa:)10846; PDE10A ####
BiocManager::install("org.Hs.eg.db", version = "3.13") # 2021
#library(org.Dr.eg.db) #Zebrafish
#DOI: 10.18129/B9.bioc.org.Hs.eg.db       
#Genome wide annotation for Human
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(clusterProfiler)
library(stringr) # extract last word
# removal of "z.score_"
gsymbol <- str_replace(colnames(mRNA_matrix), "z.score_", "")
# str_replace(), str_extract()
# by regular expression
#  "\\w+$" last word or [:alnum:]
# extract first word: "\\w+"

# p=20239 genes
#t <- c("lepa","lepr","lepb","leprot")
et <- bitr(gsymbol, fromType="SYMBOL", toType=(c("ENTREZID")), OrgDb="org.Hs.eg.db") # homo Sapiens
# ,"ALIAS", , "GENETYPE"(all are protein-coding)
# 15.94% of input gene IDs are fail to map
head(et) # p=17016
# manual convert: https://biit.cs.ut.ee/gprofiler/convert to Entrez_accesion number
# or package: library(gprofiler2)
View(gsymbol[!gsymbol %in% et$SYMBOL])
write.csv(gsymbol[!gsymbol %in% et$SYMBOL], file="~/Downloads/gsymbol.csv")
et2 <- read.csv(file="~/Downloads/gProfiler_hsapiens_6-21-2021_6-35-44\ PM.csv") # 514 converted

gsymbol <- data.frame("Gene.symbol" = gsymbol)
colnames(et)[1] <- colnames(gsymbol)[1]
colnames(et2)[c(1,2)] <- c(colnames(gsymbol)[1], "ENTREZID")
et2$ENTREZID <- ifelse(et2$ENTREZID=="NaN", NaN, et2$ENTREZID) # becoming true NaN (nil)
  
gsymbol0 <- gsymbol
# gsymbol <- gsymbol0
et3 <- rbind(et, et2[, c(1,2)])
#gsymbol <- merge(gsymbol, et2[, c(1,2)], by = "Gene.symbol", all.x=TRUE)
gsymbol <- merge(gsymbol, et3, by = "Gene.symbol", all.x=TRUE)

# p=20239
for (gi in 1:length(colnames(mRNA_matrix))) {
  # take Gene.symbol
  gi_symbol <- str_replace(colnames(mRNA_matrix)[gi], "z.score_", "")
  # convert it to ENTREZ_accesion
  colnames(mRNA_matrix)[gi] <- gsymbol$ENTREZID[gsymbol$Gene.symbol==gi_symbol] # hsa homo Sapiens
}
# p=17501
col_idx <- str_detect(colnames(mRNA_matrix), "\\d+")
mRNA_matrix0 <- mRNA_matrix
# # why only p=3765
# mRNA_matrix2 <- mRNA_matrix[, col_idx]
# ?? undefined columns selected
# 
# 1st column
mRNA_matrix1 <- subset(mRNA_matrix, select=c(colnames(mRNA_matrix)[1]))

col_idx[is.na(col_idx)] <- FALSE # eliminate the NA
# 2nd column ~
for (gi in 2:length(colnames(mRNA_matrix))) {
  if (col_idx[gi]==TRUE) {
    #print(gi)
    mRNA_matrix1 <-cbind(mRNA_matrix1, subset(mRNA_matrix, select=c(colnames(mRNA_matrix)[gi])))
    }
} # why only p=3765? (OK); error and stop at gi=4811 (NA :-)

mRNA_matrix <- mRNA_matrix1

# data cleaning####
#removal missing event status or time-to-event (2 NaN  cases at no 165 416)
table(complete.cases(mRNA_matrix))
which(is.na(survData$time))
# knock-out at 165 416: TCGA-CQ-A4CA, TCGA-H7-A6C4
mRNA_matrix[c(165, 416), 1] <- NaN
# and at 32: TCGA-BA-A6DF
survData[c(32), 2] <- NaN
survData <- survData[complete.cases(survData$time), ]
# missing at 32 x (imputation this missing data)
mRNA_matrix <- mRNA_matrix[complete.cases(mRNA_matrix$"1"), ]
save(survData, mRNA_matrix, file="TCGA_HNSCC_survival518_mRNA_matrix17501.Rda")

# n=518, p=17499
 # col by col cleaning of mRNA_matrix
 # # 1st column
mRNA_matrix3 <- subset(mRNA_matrix, select=c(colnames(mRNA_matrix)[1]))
# 2nd column
for (im in 2:length(mRNA_matrix)) {
  if (!anyNA(mRNA_matrix[, im])) {
    print(paste("complete values at column", im, "(name:",  colnames(mRNA_matrix)[im], ")", sep = ""))
    mRNA_matrix3 <-cbind(mRNA_matrix3, subset(mRNA_matrix, select=c(colnames(mRNA_matrix)[im])))
  }
}
mRNA_matrix <- mRNA_matrix3
# n=518, p=17499 -> 17345 (genes, features)
# removed all NaN missing value over 154 genes
save(survData, mRNA_matrix, file="TCGA_HNSCC_survival518_mRNA_matrix17345_noNaN.Rda")

# missing data checked by plot ####
# library(ggplot2)

x_nan <- 500
y_nan <- dim(mRNA_matrix)[1]
y_range <- c(1:y_nan)
for (i_nan in 1:(dim(mRNA_matrix)[2]/x_nan)) {
# +500*(i_nan-1)
x_range <- c((1+x_nan*(i_nan-1)):
               (x_nan+x_nan*(i_nan-1)))

# convert to vector by rows, as.vector():; unlist() by columns
color_nan <- c(t(mRNA_matrix[x_range, y_range]))
DAX <- rep(x_range, times=y_nan) # 1...10 *10
DAY <- rep(rev(y_range), each=x_nan) # reverse 10 ....1
mRNA_matrix4 <- data.frame(DAX=DAX, DAY=DAY, color_nan=color_nan)
#red_nan <- data.frame(DAX_nan=0, DAY_nan=0)
red_nan <- mRNA_matrix4[is.na(mRNA_matrix4$color_nan), 1:2]
ggplot() + 
  geom_point(data=mRNA_matrix4, aes(x=DAX, y=DAY, colour=color_nan)) + 
  theme(legend.position="none") + 
  scale_colour_gradient(low="white", high="black") +
# #A9C8F3, #0C2389
  geom_point(data=red_nan, aes(x=DAX, y=DAY), colour="red",size=5)
print(i_nan)
ggsave(file=paste("check_NaN_", i_nan, ".png", sep = ""), device = "png")

} # end of for loop i_nan
# 

# n=314, p=670 (genes, features)
#data(mRNA_matrix) # gene expression profiles
#data(survData) # survival information: 1 meaning the time is a failure time (dead), and zero a censoring time (alive).
data(dGMMirGraph) # The KEGG network constructed by the R package iSubpathwayMiner.
# list of cpd:C00399

# > training vs testing dataset ####
# cross-validation when training
set.seed(20210622) # 夏至; World ALS day
train.Idx <- sample(1:dim(mRNA_matrix)[1], floor(0.8*dim(mRNA_matrix)[1])) # sampling 80%: 416 vs 105
test.Idx <- setdiff(1:dim(mRNA_matrix)[1], train.Idx)
x.train <- mRNA_matrix[train.Idx ,] # clean6_oscc_tobacco
x.test <- mRNA_matrix[test.Idx ,] # clean6_oscc_tobacco
y.train <- survData[train.Idx,] # clean6_oscc_tobacco
y.test <- survData[test.Idx,] # clean6_oscc_tobacco
# go to cv.glmnet


# no cross validation
# glmnet (lasso) to train a mod (model)
# Calculating Cox p-value and performing random walk
#mod <- RLassoCox(x=x.train, y=y.train, globalGraph=dGMMirGraph)
#
#plot(mod$glmnetRes)

# 
#print(mod$glmnetRes) # a summary path at each step
# df represents the number of non-zero coefficients, # %Dev represents the percent (of null) deviance explained
# Lambda represents the value of λ

# The actual coefficients of genes
# at one or more λs (s=) (Why s and not lambda? :-)
#head(coef(mod$glmnetRes, s = 0.2))



# training by k-fold cross-validation
# cv.glmnet ####
# Solving or avoiding Overfitting ##
#  Cross-validation is perhaps the simplest and most widely used method for select one of models at specific 𝜆s to avoid overfitting.
cv.mod <- cvRLassoCox(x=x.train, y=y.train, globalGraph=dGMMirGraph, nfolds = 10) # or 5
#Calculating Cox p-value...Done [2021/06/22] 19:42
#Performing directed random walk...Done
#Performing cv.glmnet...Done
# # glmnet Package to fit ridge/lasso/elastic net models
# code should look like:
#fit.lasso <- glmnet(x.train, y.train, family="gaussian", alpha=1) # L1-regularization
#fit.ridge <- glmnet(x.train, y.train, family="gaussian", alpha=0) # L2-regularization
#fit.elnet <- glmnet(x.train, y.train, family="gaussian", alpha=.5) # mixed of above


# Deal with overfitting by 1) regularization (e.x. LASSO) and 2) cross validation:
# the optimal λ value and a cross validated error plot
# x-axis: log(lambda)
# y-axis: lowest CV error or highest C-index
# Model developers tune the overall impact of the regularization term by multiplying its value by a scalar known as lambda (also called the regularization rate)
# 
# # Evaluation of Overfitting ##
# For time-to-event outcome, the concordance index (C-index) Harrell et al. (1996) is the standard performance metric for model assessment in survival analysis.
# C-index [0.5-1] measures the percentage of pairs of subjects correctly ordered (risk) by the model in terms of their expected survival time.
# The predictive performances are computed in terms of accuracy and C-index, respectively for the classification and the survival prediction. 
# When there are no censored values, the C-index is equivalent to the Area Under the Curve (AUC) (in the case-control outcome scenario); or accuracy.

# selected gene signature:
# To gain some biological sense of the selected features (genes) we downloaded all gene-disease associations from DisGenet, which integrates data from expert curated repositories, KEGG, animal models and the scientific literature.


# Y-axis for the cross-validation error curve (red dotted line) (partial likelihood deviance)
#plot(cv.mod$glmnetRes, xlab = "log(lambda)")

# how to?
# find my candidate genes (features)####
# res <- RLassoCox(x=trainSmpl, y=survData[trainSmpl.Idx ,], globalGraph=dGMMirGraph)
# to yield a  glmnetRes (response) and a topological weight vector PT of genes
# It shows the path of the coefficient of each gene and L1-norm when λ varies. The axis above indicates the number of nonzero coefficients at the current λ
# lasso 就是要找出 coefficient of regression non-zero 的基因
plot(cv.mod$glmnetRes) # xlab = "log(λ)" in Greece

# ***
print(cv.mod$glmnetRes) # a summary path at each step
# df represents the number of non-zero coefficients, 
# %Dev represents the percent (of null) deviance explained
# Lambda represents the value of λ



# The first hit of  λ can be obtained:
cv.mod$glmnetRes$lambda.min
## [1] 0.4701213 -> 0.2312771
## 
## *** the most regularized model with CV-error within 1 standard deviation of the minimum.
cv.mod$glmnetRes$lambda.1se # my candidates: a gene signature with 14 genes
## [1] 0.7842097 -> 0.5863717
## 


# The actual coefficients of genes
# at one or more λs (s=) (Why s and not lambda? :-)
#head(coef(cv.mod$glmnetRes, s = "lambda.min"))

# if lambda.min
# 62 genes selected by it's coefficient
coef.min <- coef(cv.mod$glmnetRes, s = "lambda.min") 
View(coef.min)
# thus, the selected features (62 genes) and their coefficients:
nonZeroIdx <- which(coef.min[,1] != 0) 
features <- rownames(coef.min)[nonZeroIdx] 
#features == Entrez accession
features.coef <- coef.min[nonZeroIdx] 
names(features.coef) <- features # named vector
features.coef <- as.data.frame(features.coef)

# converting to Gene.symbol
features.coef2 <- bitr(rownames(features.coef), 
                          fromType="ENTREZID",
                          toType=(c("SYMBOL")),
                          OrgDb="org.Hs.eg.db") 
# candidates under lambda.min
cv.mod.candidate <- cbind(features.coef, features.coef2)
colnames(cv.mod.candidate) <- c("features.coef", "ENTREZID", "Gene.symbol")


# if lambda.1se: the most regularized model with CV-error within 1 standard deviation of the minimum.
# 14 genes
coef.1se <- coef(cv.mod$glmnetRes, s = "lambda.1se") 
View(coef.1se)
# thus, the selected features (62 genes) and their coefficients:
nonZeroIdx <- which(coef.1se[,1] != 0) 
features <- rownames(coef.1se)[nonZeroIdx] 
#features == Entrez accession
features.coef <- coef.1se[nonZeroIdx] 
names(features.coef) <- features # named vector
features.coef <- as.data.frame(features.coef)

# converting to Gene.symbol
features.coef2 <- bitr(rownames(features.coef), 
                       fromType="ENTREZID",
                       toType=(c("SYMBOL")),
                       OrgDb="org.Hs.eg.db") 
# candidates under lambda.1se
cv.mod.candidate.1se <- cbind(features.coef, features.coef2)
colnames(cv.mod.candidate.1se) <- c("features.coef", "ENTREZID", "Gene.symbol")
cv.mod.candidate.1se[order(cv.mod.candidate.1se$features.coef, decreasing = TRUE), 3]
#[1] "POLR2C" "EFNB2"  "PLAU"   "HK1"    "WNT7A"  "HPRT1"  "APP"   
#[8] "SC5D"   "EZR"    "RRAGA"  "S1PR4"  "DUSP16" "CRLF2"  "SH3BP2"
save(cv.mod.candidate.1se, cv.mod.candidate, file="RLassoCox_TCGA_candidate14_signatures.Rda")
# final 4 gene signature: HK1, EZR, RRAGA, and DUSP16


# TCGA
# [2021/06/22] code copied from Validation_GSE65858_survival.R
#> TCGA Top3 CAMK2N1/IL19/FCGBP and cv.mod.candidate (64) in TCGA #### 
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
gene_labeled_golden <- c("CAMK2N1", "IL19", "FCGBP")
gene_labeled_silver <- c(
  "DKK1", #CAMK2N1, 
  "STC2", "PGK1", "SURF4", "USP10", "NDFIP1", "FOXA2", "STIP1", "DKC1", 
  "ZNF557", "ZNF266", #IL19, 
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
# 62 genes
load(file="RLassoCox_TCGA_candidate62.Rda")
gene_labeled <- cv.mod.candidate$Gene.symbol
#
# 4 gene signature: #### 
# "EZR" has P value >0.05
gene_labeled_4 <- cv.mod.candidate$Gene.symbol[c(17, 23, 28, 52)]
# TCGA HNSCC
# # going for volcano plot by ggplot() ####
library(ggplot2)
library(ggrepel) # function "geom_label_repel"
cut_pvalue <- 0.05
ggplot(hazards, aes(y = multi_HR, x = plog10P)) +   # x = plog10P  #Creation of ggplot graph with pre-determined aesthetics
  geom_point(shape = 21, alpha = 0.8, aes(size = log10P, fill = multi_HR)) + 
  scale_fill_distiller(palette = "RdYlGn", trans="log", limits = c(min(hazards$multi_HR), max(hazards$multi_HR))) + # min(hazards$uni_HR)
  #Diverging  BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral
  geom_vline(xintercept = log10(cut_pvalue), color="red", lty=2) + 
  geom_hline(yintercept = 0.5, lty = 2) +
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
    expand = expansion(mult = 0.2)) +
  
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
  ) +  
  # ### ggplot() for gene_labeled_silver
  # # hazards
  # geom_label_repel(data=subset(hazards, Gene.Symbol %in% gene_labeled_silver),
  #   aes(label=Gene.Symbol),
  #   color = "black",
  #   size = 3,
  #   nudge_x       = -3.2 - subset(hazards, 
  #                                 Gene.Symbol %in% gene_labeled_silver)$plog10P,
  #   nudge_y       = 1.2 - subset(hazards, 
  #                                Gene.Symbol %in% gene_labeled_silver)$multi_HR,
  #   segment.size  = 0.2,
  #   segment.color = "blue",
  #   direction     = "y",
  #   hjust         = 0,
  #   max.overlaps = Inf
  # ) +
  
  # gene_labeled (other than CAMK2N1) => 62 Bayes survival (Cox-Lasso) model
  # *** 4 gene signature
  # hazards
  geom_text_repel(data=subset(hazards, Gene.Symbol %in% gene_labeled_4),
                  aes(label=Gene.Symbol),
                  color = "grey50", #factor(golden)), #  == TRUE
                  #                   data          = subset(dat, wt < 3),
                  nudge_x       = 1 - subset(hazards, Gene.Symbol %in% gene_labeled_4)$plog10P,
                  #                   segment.color = "grey50",
                  direction     = "y",
                  hjust         = 1,
                  
                  # x  color = ifelse(hazards$Gene.Symbol %in% gene_labeled_golden, "red", "grey50"),
                  #                   scale_color_manual(values = cols),
                  # breaks = c("TRUE", "FALSE"),
                  # guide = "none"), # better than color=
                  size=5, box.padding = unit(0.4, "lines"), 
                  segment.alpha=0.8, segment.size = 0.5,
                  segment.color = "black",
                  max.overlaps = Inf)  + 
  theme_bw() + 
  theme(text=element_text(family="sans"),
        axis.title=element_text(size=18))


##
## EFNB2, PLAU, LASP1
##
ggsave("Rplot_TCGA_HNSCC_CoxHR_CAMK2N1_signature4_FDRKM.pdf") #, width = 12, height = 8, dpi = 84)


# GSE65858
#> ggplot(Cox uni_HR of GSE65858) ####
# 62 genes 
load(file="RLassoCox_TCGA_candidate62.Rda")
gene_labeled <- cv.mod.candidate$Gene.symbol

gene_labeled <- gene_labeled[c(-7,-36,-62,-35,-60,-9)] # "APP", "YWHAG", "LIF" "LASP1" "TPM3" "ATF6" removal
load(file="HNSCC_OS_GSE65858_pvalueCox_HR.Rda") # uni_HR => osHR
hazard0 <- subset(osHR, select=c(Gene.Symbol, uni_HR)) # Cox.P-value 
load(file="fdr7614_KMsurvival_GSE65858.Rda") # as fdr100_KMsurvival
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

gene_labeled_golden <- c("CAMK2N1", "IL19", "FCGBP")

# from pvalueTex without gene_labeled_golden
# # **** skip below ***
# gene_labeled <- c(
#   #  "DKK1", 
#   #"CAMK2N1", 
#   # "STC2", "PGK1", "SURF4", #"USP10", "NDFIP1",
#   #"FOXA2", "STIP1", "DKC1", 
#   #                  "ZNF557", "ZNF266", 
#   #"IL19", "FCGBP",
#   # "MYO1H", "LOC148709", "EVPLL", "PNMA5", "KIAA1683", 
#   #"NPB",
#   
#   #    "ATP13A4",    
#   "PLAU",      # it is in Bayes Cox-Lasso
#   #"ZNF662",     
#   #"POMP",       
#   #"DOT1L",     
#   #"EPHX3" ,     #"IL34",       
#   "FAM3D",     
#   "MSMB", "LYPD2", "RBM11", # *** top KM P value in GSE65858 
#   #"SERPINE1",   #"CALML5",    
#   # "RASIP1",    # "FUT6",       
#   #"BCAR3",     
#   # "ST6GALNAC1", #"AQP1",   
#   #"AIG1",        
#   # "FAM151B"
#   # "ABCB1",      
#   #"SMPX"
#   #"MASP1" ,   
#   # "GRIA3"
#   # CoxHR > 1.81:
#   "ADAMTSL2",
#   "BMP6" , "DUSP6" ,
#   "FXYD5" ,   "GPX8",
#   #"HERPUD2",
#   "IFIT2" ,  
#   #"KANK4",   
#   "MYO1B" ,
#   "NFKBID",   "PDIA5" ,   "PMEPA1",  "RAB43"  , 
#   # "RET" ,  
#   "SPP1" , 
#   "TCP11L1",  "TIMM10"  ,
#   "TNFAIP6" , # "TNPO1"  ,
#   #"TPD52L2", 
#   "TRAPPC10",
#   "TTYH3" ,   "XYLT1" 
# )

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



# * FDR correction of KM P value of GSE65858 ##
# unadjusted P value < 0.05 (534 probes) -> fdr
#p.adjust(HNSCC_OS_marginS_pvalue_sorted$p_value, method="bonferroni")
hazards$FDR.KMpvalue[hazards$KM.Pvalue < 0.05] <- 
  p.adjust(hazards$KM.Pvalue[hazards$KM.Pvalue < 0.05], method = "fdr")
#R4> hazards[hazards$Gene.Symbol %in% gene_labeled_golden, c(1,3,10,6,2)]

#save(hazards, file="HNSCC_GSE65858_FDRKMpvalue_7614probes.Rda")

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
  scale_y_continuous(
    expand = expansion(mult = 0.2)) +
  
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
                  size=2, box.padding = unit(0.4, "lines"), 
                  segment.alpha=0.8, segment.size = 0.1,
                  max.overlaps = Inf)  + 
  theme_bw() + 
  theme(text=element_text(family="sans"),
        axis.title=element_text(size=18))
# thanks https://ggrepel.slowkow.com/articles/examples.html
ggsave(file="Rplot_GSE65858_HNSCC_CoxHR_CAMK2N1_Bayes62_FDRKM.pdf")
# "APP" removed for clear show


## TCGA KM plot of signature ####
## R code from SurvExpress
## version <- "runweb.r v1.13 Sep/111/2015"
## by http://bioinformatica.mty.itesm.mx:8080/Biomatec/SurvivaX.jsp
## Concordance= 0.633 (C-index)
# Symbol Id	Beta	p Value
# HK1	1.413	0.002069226
# EZR	0.973	0.016782975
# RRAGA	0.925	0.021660706
# DUSP16	-0.745	0.035048687

# 

# plot.kaplan prepare ####
# # call here
# if (xfunc == "cox-survival-groups") 有二段
# plot.kaplan(data[,tr,drop=FALSE], time=time[tr], status=status[tr], model=1:nrow(data), risk.groups=groups, 
#             col=xcolor, col2=survcolor,
#             main=paste(paste(predictor,ifelse(train.is.dif,"(train only)",""), ifelse(useWeights,"(not using weights)","")),paste(symb,collapse="|"),sep="\n"), 
#             symbols=symb,
#             draw.text=fitting.info,
#             draw.main=fitting.info)
data0 <- data
tr0 <- tr
time0 <- time
status0 <- status

# symb as "HK1"    "EZR"    "RRAGA"  "DUSP16" signature
#save(data0, time0, status0, tr0, symb, file="SurvExp_signature4.Rda")
load(file="SurvExp_signature4.Rda")

data <- data0[,tr0,drop=FALSE]
time <- time0[tr0]
status <- status0[tr0]
model <- 1:nrow(data0)
tr <- 1:ncol(data0)
tr.betas <- tr


risk.groups=2
col=1:risk.groups+1
col2=1:risk.groups+1
# predictor
main=paste(paste("KM plot",ifelse(train.is.dif,"(train only)",""), ifelse(useWeights,"(not using weights)","")), "for ", paste(symb,collapse="|"),sep="")

roundcoef<-3
symbols<-symb
betas<-NULL
betas.p<-rep(1, length(betas))
draw.text<-TRUE
draw.main<-TRUE
draw<-TRUE
needMaximizeGroups<-FALSE
univariateBetas<-FALSE

# removal of function rapping:
# plot.kaplan <- function(data, time, status, model, tr=1:ncol(data),
#                         risk.groups=2, col=1:risk.groups+1,
#                         col2=1:risk.groups+1,
#                         main="", 
#                         tr.betas=tr, roundcoef=3, symbols=NULL, 
#                         betas=NULL, betas.p=rep(1, length(betas)),
#                         draw.text=TRUE, draw.main=TRUE, draw=TRUE,
#                         needMaximizeGroups=FALSE, univariateBetas=FALSE) {
  #col[1:risk.groups] <- col[risk.groups:1]
  library(survival)

# glmnet handles ties in survival time with the Breslow approximation. 
# This is different from survival package’s coxph function, whose default tie-handling method is the Efron approximation
  ci <- NA
  betas.supplied <- ! is.null(betas)
  coxp <- betas.p
  estimated <- FALSE
  cox <- NULL
  try (
    cox <- coxph(Surv(time[tr.betas], status[tr.betas]) ~ .,data.frame(t(data[model,tr.betas,drop=FALSE])), method="breslow")
  ) # or ties = "breslow"
  if (is.null(betas)) {
    estimated <- is.null(cox)
    if (estimated  ||  univariateBetas) {
      betas <- c()
      coxp <- c()
      for (i in 1:length(model)) {
        cox <- coxph(Surv(time[tr.betas], status[tr.betas]) ~ x,data.frame(x=data[model[i],tr.betas]), method="breslow")
        betas <- c(betas, cox$coefficients)
        coxp <- c(coxp, summary(cox)$coefficients[,5])
      }
    } else {
      betas <- cox$coefficients
      coxp <- summary(cox)$coefficients[,5]
    }
  }
  coxo <- cox
  scoxo <- summary(coxo)
  if (any(is.na(betas))) 
    betas[is.na(betas)] <- 0
  ci <- concordance.index(concordance.cox(betas, time[tr], status[tr], data[model, tr, drop=FALSE]))
  p.i <- betas %*% data[model, , drop=FALSE]
  
  maxRisk <- toupper(getParam("max_risk", "YES")) %in% c("1","YES","TRUE")
  if (maxRisk || needMaximizeGroups) {
    mrg <- maximizeRiskGroups(time, status, p.i, risk.groups, tr=tr.betas)
    pi.cluster <- if (maxRisk) mrg$cluster.max else mrg$cluster.quantile
  } else {
    qcuts <- quantile(p.i, seq(0,1,length.out=risk.groups+1, na.rm=TRUE))
    pi.cluster <- cut(p.i, breaks=qcuts, labels=FALSE, include.lowest=TRUE)
    mrg <- NULL
  }
  
  ### hazard-ratio estimation of the PI as predictor, THIS WILL ALWAYS BE 1
  ### see Assessment of performance of survival prediction models for cancer prognosis - Chen - Chen - BMC Medical Research Methodology 2012.pdf
  ### Instead, it is better to estimate the hazard ratio by the two-groups
  ### That is, fitting a cox model to the cluster as predictor
  #hr <- coxph(Surv(time[tr.betas], status[tr.betas]) ~ PI, data.frame(PI=p.i[tr.betas]), method="breslow") 
  hr <- coxph(Surv(time, status) ~ cluster, data=data.frame(time=time[tr], status=status[tr], cluster=pi.cluster[tr]), method="breslow") 
  shr <- summary(hr)
  
  
  surdif <- survdiff(Surv(time, status) ~ cluster, data=data.frame(time=time[tr], status=status[tr], cluster=pi.cluster[tr]))
  p <- 1-pchisq(surdif$chisq,df=risk.groups-1)
  
  vars <- if (is.null(rownames(data))) model else rownames(data)[model]
  xm <- data.frame("Symbol:Id"=vars)
  xm$Beta <- round(betas,roundcoef)
  xm$pValue <- coxp
  xcol <- matrix(1, ncol=risk.groups*2+3, nrow=length(model)+1)
  xcol[-1,2:3] <- ifelse(coxp < 0.05, 16, ifelse(coxp < 0.1, 8, 1))
  xcol[ 1,2:3] <- 16
  res <- list()
  cox <- list()
  cox.ci <- list()
  cox.cirefit <- list()
  censxrisk <- rep(0,risk.groups)
  survcox <- list()
  xnames <- c()
  for (k in 1:risk.groups) {
    xcol[1,2:3+k*2] <- col[k]
    w <- which(pi.cluster[tr] == k)
    censxrisk[k] <- sum(status[tr[w]] == 0)
    coxk <- NULL
    xnames <- c(xnames, paste(c("Beta","pValue"),k,sep=""))
    #if (draw.text) {
    try(coxk <- coxph(Surv(time[tr[w]], status[tr[w]]) ~ ., data.frame(t(data[model, tr[w], drop=FALSE])), method="breslow"))
    #}
    cox.ci[[k]] <- cox.cirefit[[k]] <- NA
    if(!is.null(coxk)) {
      cox[[k]] <- coxk 
      survcox[[k]] <- NULL
      try(survcox[[k]] <- survfit(cox[[k]]))
      res[[k]] <- data.frame(coef=round(cox[[k]]$coefficients,roundcoef), 
                             p=summary(cox[[k]])$coefficients[,5])
      xm <- cbind(xm, res[[k]])
      xcol[-1,2:3+k*2] <- ifelse(res[[k]]$p < 0.05, col[k], ifelse(res[[k]]$p < 0.1, rgb2rgb(col2rgb(col[k])/2), 1))
      ## con modelo del refit
      try (cox.cirefit[[k]] <- concordance.index(concordance.cox(coxk, time[tr[w]], status[tr[w]], data[model, tr[w], drop=FALSE])) )
    } else {
      cox[[k]] <- NA
      xm <- cbind(xm, matrix(NA, ncol=2, nrow=length(model)))
    }
    ## con modelo "original" (de todo el grupo)
    try (cox.ci[[k]] <- concordance.index(concordance.cox(coxo, time[tr[w]], status[tr[w]], data[model, tr[w], drop=FALSE])) )
  }
  colnames(xm) <- c(colnames(xm)[1:3], xnames)
  
  #r2 <- scoxo$rsq["rsq"] 
  #r2max <- scoxo$rsq["maxrsq"]
  ## if betas.supplied then tr is supposed to be test in this case, so r2 is about the correlation in test
  ## these 2 functions made the estimation in any case
  r2 <- cox.r2(time[tr], status[tr], betas, data[model, tr, drop=FALSE])
  r2max <- cox.r2.max(time[tr], status[tr], betas, data[model, tr, drop=FALSE])
  
  pp = par(mar=c(5.1+max(0,risk.groups-2), 4.1, 5.1, 2.1))
  on.exit(par(pp))
  
# plot hereafter ####
# code from Validation_GSE65858
  pdf_fn <- paste("rplot_KMplot_HK1.pdf", sep="")
  pdf(pdf_fn) 
  # save PDF, size 746 x 431 pixel(?)
surv_signature4 <-  survfit(Surv(time, status) ~ cluster, 
          data=data.frame(time=time[tr], status=status[tr], cluster=pi.cluster[tr]))

##

##

plot(surv_signature4, lty=1, xscale=365.25, xmax=6000, col=c("blue","red"), 
       sub=paste("Kaplan-Meier P Value =", format(p, digits=3)), 
       main=paste("OS in TCGA for ", paste(symb,collapse="|"),sep=""), 
       ylab="Percent Survival", 
       xlab="Years")
# legend
legend("topright", legend=c(paste("low(",surv_signature4$n[1], ")"), paste("high(",surv_signature4$n[2], ")")),
         lty=1:1, col=c("blue","red"))
  
dev.off() 
#
# go line 500 for TCGA plot


##
# spared code from SurvExpress (a function)
#  if (draw && length(unique(pi.cluster[tr])) > 1)
plot(survfit(Surv(time, status) ~ cluster, 
             data=data.frame(time=time[tr], status=status[tr], cluster=pi.cluster[tr])), 

###
getParam <- function(pname, default=NULL) {
xV <- params[pname]
if (is.null(xV)  ||  is.na(xV)  || (is.character(xV) && length(xV) == 1 && nchar(xV) == 0))
  xV <- params[make.names(pname)]
if (is.null(xV)  ||  is.na(xV) ||  (is.character(xV) && length(xV) == 1 && nchar(xV) == 0))
  xV <- default
xV
}  #end of function
####################i###








########################3e###
#x skip predict ####
#lp <- predict(object = res, newx = testSmpl)
# object is trained Cox lasso model
# k-fold validation by cvRLassoCox
# cv.mod <- cvRLassoCox(x=x.train, y=y.train, globalGraph=dGMMirGraph, nfolds = 10) # or 5

#lp <- predict(object = mod, newx = x.test, s = c(0.1, 0.2)) 
lp <- predict.cvRLassoCox(object = cv.mod, newx = x.test, s = cv.mod$glmnetRes$lambda.min)
head(lp)

