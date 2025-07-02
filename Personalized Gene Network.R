rm(list = ls())

library(Dozer)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(dplyr)
library(knitr)
library(foreach)
library(doParallel)
library(cluster)
library(Rtsne)
library(R.matlab)
library(Rcpp)
library(igraph)
library(mstknnclust)

############################################
## Prepare Network Dataset For Testing
############################################
load(system.file("extdata", "Jerber_demo.rda", package = "Dozer"))
theme_set(theme_pubr(base_size = 12))

path = paste0(tempdir(), '/dozer_tmp')
if (! file.exists(path)){
  dir.create(path)
}
cl <- makeCluster(detectCores()) 
registerDoParallel(cl)
## Load noise ratio into a "gene by donor" matrix.
noise_ratio_gene_by_donor = foreach (i = 1:nrow(donor_info), .packages = c('Matrix'), .combine = 'cbind') %dopar%{
  donor = donor_info$donor_id[i]
  data = matrix(counts[ , metadata$donor_id == donor], nrow = nrow(counts))
  meta = metadata[metadata$donor_id == donor, ]
  ribo = unlist(lapply(rownames(data), FUN=function(x){substr(x,1,3)}))%in%c('RPL','RPS')
  ribosomal_perc = colSums(data[ribo,])/colSums(data)
  ## If there are several sample_index presented in one dataset, regress it out.
  res = Dozer::compute_gene_correlation(data, covs = data.frame(batch_label = meta$sample_id, nFeature = colSums(data>0), ribosomal_perc))
  ## Saver co-expression matrices to file
  save(res, file = paste0(path, donor, '-coexpr.rda'))
  res$ratio[,1]
}
stopCluster(cl)

keep_gene = rowMeans(noise_ratio_gene_by_donor) < .9

gene_name = rownames(counts)[keep_gene]
## The number of genes passed filtering with noise ratio.
sum(keep_gene)
#> [1] 1130

PersonalGeneNet = list()
for(i in 1:nrow(donor_info)){
  donor = donor_info$donor_id[i]
  load(paste0(path, donor, '-coexpr.rda'))
  # hard-thresholding
  network_i = abs(res$network[keep_gene, keep_gene])
  q = quantile(network_i[upper.tri(network_i)], .95)
  network_i[network_i<q] = 0
  network_i[network_i>0] = 1
  PersonalGeneNet[[i]] = network_i
}


############################################
## Load NEPTUNE Functions
############################################
code_dir  = "/Users/hanfeng/Dropbox/Miao/Aim1 Paper/V5 Resubmit/Code/"                   ### Directory for NEPTUNE_Main
Sys.setenv("CPLUS_INCLUDE_PATH" = "/Library/Developer/CommandLineTools/SDKs/MacOSX15.2.sdk/usr/include/c++/v1")
source(paste0(code_dir, "NEPTUNE_Updated_Main.R"))

############################################
## Below Should be Input
############################################
## Compute average networks in each donor group
group = NA
idx1 = which(donor_info$phenotype=='Success')
idx0 = which(donor_info$phenotype=='Failure')
group[idx1] = 1 
group[idx0] = 0

table(donor_info$phenotype)
prop.table(table(donor_info$phenotype))
table(group)
#####
N_Sample1 = length(idx0)    ## Sample size 1
N_Sample2 = length(idx1)    ## Sample size 2
N_Total = N_Sample1 + N_Sample2       ## Total Sample Size

## Read Data
WorkSample = list()
for (i_case in 1:N_Sample1){
  print(idx0[i_case])
  WorkSample[[i_case]] = PersonalGeneNet[[idx0[i_case]]]
}
for (i_case in 1:N_Sample2){
  print(idx1[i_case])
  WorkSample[[(N_Sample1 + i_case)]] = PersonalGeneNet[[idx1[i_case]]]
}

############################################
## Main Function
############################################
set.seed(11)
n_perm = 1000  ## Empirical p Resolution
delta  = 0.5  ## Correction term for 0

### Get Parmters from Data  
Vsize = nrow(WorkSample[[1]])                       ## Vertex Size from Data
Whalf = Half_Width(Vsize)                           ## Calculate Half Width Based on Vertex Size
Eigen_Mat = matrix(0, nrow = Vsize, ncol = N_Total) ## Matrix Save Eigen Values From Subjects by Column to Avoid Repeated Calculation

for (idx in 1:N_Total){
  Eigen_Mat[, idx] <- Lap_Eigen(WorkSample[[idx]])
}

## QIM Distance
Sample_QIM <- QIMMatrix_Calc(WorkSample, N_Total, N_Sample1, Whalf, Eigen_Mat, l=1)
temp       <- PCalc(Sample_QIM, N_Sample1, n_perm)
p_QIM  <- (sum(temp > temp[1]) + delta)/(n_perm + delta)
p_QIM

### Highdim Distance
Sample_MR  <- MR_MAT(Sample_QIM)
temp       <- PCalc(Sample_MR, N_Sample1, n_perm)
p_MR   <- (sum(temp > temp[1]) + delta)/(n_perm + delta)
p_MR

############################################
## 2D Illustrations
############################################
set.seed(68)
result_mstknn <- mst.knn(Sample_QIM)
plot(result_mstknn$network, vertex.size=18,
     vertex.color=(workindex+3),
     layout=igraph::layout.fruchterman.reingold(result_mstknn$network, niter=10000),
     main=paste("MST-kNN Clustering Based on QIM \n Colored by Phenotype",sep="" ),
     vertex.label = NA)

set.seed(68)
result_mstknn <- mst.knn(Sample_MR)
plot(result_mstknn$network, vertex.size=18,
     vertex.color=(workindex+3),
     layout=igraph::layout.fruchterman.reingold(result_mstknn$network, niter=10000),
     main=paste("MST-kNN Clustering Based on MR \n Colored by Phenotype",sep="" ),
     vertex.label = NA)

############################################
## Main Function (Or simply call)
############################################
a = NEPTUNE_Test(WorkSample, N_Sample1, N_Sample2, n_perm = 1000, delta = 0.5, seed=11)