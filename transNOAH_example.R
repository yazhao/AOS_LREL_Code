###############################################
# This script takes the transNOAH cancer data set 
# and applies the test by calculating the 
# corresponding T statistic from
# the LREL project, and then proceeds
# to draw conclusions. Additionally, it applies
# comparison tests from Zhong and Chen 2011 and
# Cui et al (2018)
###############################################

  # and was wondering if it would be possible for you to help me better understand the numerical example especially with regards to the data set used:
  # Your paper mentions that "We apply the method developed in Section 2 to this dataset with the goal of testing particular associations between BRCA1 and a few other genes conditional on all remaining genes present in the study" with results reported in Table 3. Are the Genes listed in Table 3 comprehensive (ie you only tested associations between BRCA1 and those genes)?To my understanding, the transNOAH data has several sequences that are associated with each gene (there are 2 sequences associated with BRCA1, for example). Were the expression of these sequences combined in some way, or were specific ones picked for each gene, or was some other procedure used?

# libraries
library(GEOquery)
library(data.table)
source("simulation_code/simulations_functions.R")

# data loading for transNOAH to be done once and saving
Sys.setenv("VROOM_CONNECTION_SIZE" = 500000)
readr::local_edition(1)
transNOAH = GEOquery::getGEO("GSE50948")
save(transNOAH, file = "transNOAH.RData")

# loading semi-processed transNOAH data
load("transNOAH.rData")

# processing data
tn = transNOAH$GSE50948_series_matrix.txt.gz
data = assayData(tn)$exprs
View(Biobase::fData(tn))

GOstep1 = Biobase::fData(tn)$`Gene Ontology Molecular Function`
uniqueGO1 = unique(unlist(strsplit(unlist(GOstep1), "[^0-9]+")))
uniqueGO = sort(uniqueGO1[nchar(uniqueGO1) == 7])

# get only the subset of info where gene title contains cancer
#xsubset = Biobase::fData(tn)[which(Biobase::fData(tn)$`Gene Title` %like% "*cancer*"|Biobase::fData(tn)$`Target Description` %like% "*cancer*"),]$ID
GOterms = rep(0,length(uniqueGO))
for(i in 1:length(uniqueGO)){
  xsubset = Biobase::fData(tn)[which(Biobase::fData(tn)$`Gene Ontology Molecular Function` %like% uniqueGO[i]),]$ID
  GOterms[i] = length(xsubset)
}
HDGO = data.frame("GO_term" = uniqueGO[which(GOterms > 158)], "terms" = GOterms[which(GOterms > 158)])
HDGO = HDGO[order(HDGO$terms),]

# loading the data
yrows = which(rownames(data) %in% c("204531_s_at", "211851_x_at"))

y1 = data[yrows[1],]
y2 = data[yrows[2],]
ybar1 = mean(y1)
ybar2 = mean(y2)

x1 = data[-(which(rownames(data) %in% c("204531_s_at", "211851_x_at"))),]

results = data.frame("GO_ID" = HDGO$GO_term, "LREL_y1" = NA, "LREL_y2" = NA, "ZC_y1" = NA, "ZC_y2" = NA,  "LREL_y1_p" = NA, "LREL_y2_p" = NA, "ZC_y1_p" = NA, "ZC_y2_p" = NA)

for(j in 1:125){#nrow(HDGO)){
  
  xpick = Biobase::fData(tn)[which(Biobase::fData(tn)$`Gene Ontology Molecular Function` %like% HDGO$GO_term[j]),]$ID

  x = x1[which(rownames(x1) %in% xpick),]
  n = ncol(x)
  c = nrow(x)/ncol(x)
  p = nrow(x)
  xbar = rowSums(x)/n
  
  z.mat1 = matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  
  # generating the z matrix from the real data
  for(i in 1:n){
    z.mat1[,i] = (x[,i] - xbar) * (y1[i] - ybar1)
  }
  
  z.mat2 = matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
  
  # generating the z matrix from the real data
  for(i in 1:n){
    z.mat2[,i] = (x[,i] - xbar) * (y2[i] - ybar2)
  }
  
  # Calculating the l, k, and alpha needed for the T and T_0
  l.n = n^(5/4) * log(n)
  k.n = ((p)/log((p)))^(1/2)
  alpha = rep(1, (p))/sqrt((p))
  
  # calculating the T statistics for y1 and y2
  Tstats = calculate_statistic_new(z = z.mat1, x = x, y = y1, k = k.n, alpha = alpha, T0 = FALSE)
  results$LREL_y1[j] = Tstats[1]
  results$LREL_y1_p[j] = 2*(1-pnorm(abs(Tstats[1])))
  # 1.671916, p-value 0.09454082 for cancer stuff only
  
  Tstats2 = calculate_statistic_new(z = z.mat2, x = x, y = y2, k = k.n, alpha = alpha, T0 = FALSE)
  results$LREL_y2[j] = Tstats2[1]
  results$LREL_y2_p[j] = 2*(1-pnorm(abs(Tstats2[1])))
  # 20.53035, p-value basically 0 for cancer stuff only
  
  # Comparisons
  # Zhong and Chen 2011
  TstatsZC1 = zhongchen2011(X = x, y = y1, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE)
  results$ZC_y1[j] = TstatsZC1[2]
  results$ZC_y1_p[j] = 2*(1-pnorm(abs(TstatsZC1[2])))
  # 1.835148, p-value 0.06648373 for cancer stuff only
  
  TstatsZC2 = zhongchen2011(X = x, y = y2, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE)
  results$ZC_y2[j] = TstatsZC2[2]
  results$ZC_y2_p[j] = 2*(1-pnorm(abs(TstatsZC2[2])))
  # 25.1242, p-value basically 0 for cancer stuff only
}

results$lsig1 = ifelse(results$LREL_y1_p < 0.05, "Yes","No")
results$lsig2 = ifelse(results$LREL_y2_p < 0.05, "Yes","No")
results$zsig1 = ifelse(results$ZC_y1_p < 0.05, "Yes","No")
results$zsig2 = ifelse(results$ZC_y2_p < 0.05, "Yes","No")
results$diff1 = ifelse(results$lsig1 == results$zsig1, "No", "Yes")

# Cui 2018
# RUN MANUALLY BECAUSE OF RCV
#Tstats3 = cui2018(X = x, y = y, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE)
#Tstats3[2]
#2*(1-pnorm(abs(Tstats3[2])))
#normalized.T
#182.4385
#2*(1-pnorm(abs(182.4385)))