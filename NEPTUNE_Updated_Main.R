############################################
## Functions Section
############################################
library(Rcpp)

## Debugging
# N is the total sample size

# Integral for the Spectrum Distance Calculation
cppFunction('double integral(NumericVector b1, NumericVector b2, double whalf){
              int n = b1.size();
              double total = 0;
              for(int i=0; i<n; ++i){
                for(int j=0; j<n; ++j){
                  if((std::abs(b1[i]- b2[j]) > 0)){
                    total += ((M_PI - atan(-b1[i]/whalf) - atan(-b2[j]/whalf)) / whalf
                    + log((pow(whalf, 2) + pow(b1[i], 2)) / (pow(whalf, 2) + pow(b2[j], 2))) / (b1[i] - b2[j])) 
                    / (pow(b1[i] - b2[j], 2) + 4 * pow(whalf, 2));
                  } else {
                    total += ((M_PI / 2 - atan(-b1[i] / whalf)) 
                    + (whalf * b1[i] / (pow(whalf, 2) + pow(b1[i], 2)))) 
                    / (2 * pow(whalf, 3));
                  }      
                }
              }
              return total;
            }')

# Calculate Hamming Distance for Two Matrices
## Input: AI1, AI2 (N by N adjacency matrix of networks) 
## Output: Hamming Distance (Number) 
HamD_Calc <- function(AI1, AI2){
  nrow <- nrow(AI1)
  adjcons <- max(max(AI1), max(AI2)) - min(min(AI1), min(AI2))
  # Calculate Elementwise Average Absolute Difference 
  return (sum(abs(AI1 - AI2)) / nrow / (nrow - 1)) / adjcons
}

# Calculate Laplacian Matrix for Given Matrix and Obtain Ordered Eigen Values
## Input: AI1 (N by N adjacency matrix of networks)
## Output: List containing: evalue (eigen values), Lap (Laplacian matrix)  
Lap_Eigen <- function(AI1){
  ## Calculate Laplacian Matrix 
  LI1 <- diag(rowSums(AI1)) - AI1
  ##eigen decomposition on laplacian matrix
  v <- eigen(LI1, symmetric=TRUE, only.values = TRUE)$values
  return (v)
}

## Calculating the distance of empty and full network with given half-width (called by Half_Width)
## Input: gamma (number, half-width), N (number, node size of networks)
## Output: specdis_max (number, need to make it equal to 1 to calculate half-width)
Max_Dis_Calc <- function(gamma, N){
  ## This formula is found in the appendix of 
  ## <The HIM glocal metric and kernel for network comparison and classification>
  intgrl <- pi / 2 + atan(sqrt(N) / gamma)
  specdis_max <- (pi * gamma)^(-1) + (2 * gamma)^(-1) * (intgrl)^(-2) * 
    (gamma * sqrt(N) / (gamma^2 + N) + intgrl) - 
    (4 * gamma) / (pi * intgrl) / (4 * gamma^2 + N) * 
    (intgrl + pi / 2 - gamma / sqrt(N) * 
       log(gamma^2 / (gamma^2 + N))) - 1
  return(specdis_max)
}


## Calculating half width of matrix A
## Input: A (N by N adjacency matrix of networks), wmax (number, start value of half-width) 
## Output: wtemp (number, half-width that makes distance of empty and full network be 1)
Half_Width <- function(N, wmax = 1){
  # binary search
  wtemp <- 0.5
  # case that true half-width greater than 1
  if (Max_Dis_Calc(wmax, N) > 0){
    print("need change maxwidth from 1 to bigger value!")
  } else {
    wmin <- 0
    delta <- Max_Dis_Calc(wtemp, N)
    while(abs(delta) > 10^-15){
      if (delta > 0){
        wmin  <- wtemp
        wtemp <- (wmax + wmin)/2
        delta <- Max_Dis_Calc(wtemp, N)
      } else {
        wmax  <- wtemp
        wtemp <- (wmax + wmin)/2
        delta <- Max_Dis_Calc(wtemp, N)
      } 
    }
    return (wtemp)
  }
}

## Calculating Ipsen-Mikhailov distance
## Input: whalf (number, half-width calculated by Half_Width), 
##        AI1, AI2 (N by N adjacency matrix of networks) 
## Output: Ipsen-Mikhailov distance (Number)
## This formula is found in the appendix of 
## <The HIM glocal metric and kernel for network comparison and classification>
IpsD_Calc <- function(whalf, lam_Mat, idx1, idx2){
  lam1 <- lam_Mat[, idx1]
  lam2 <- lam_Mat[, idx2]
  lam1[lam1 < 10^-15] <- 0
  lam2[lam2 < 10^-15] <- 0
  # Calculate Eigen values and discard the smallest one (which is 0)
  w1 <- sqrt(lam1[-length(lam1)])
  w2 <- sqrt(lam2[-length(lam2)])
  ## Due to the form of [(w1, w2)T (w1, w2)], divided to 3 main parts
  Inter_sum    <- 0
  Quadra_sum_1 <- 0
  Quadra_sum_2 <- 0
  # Calculate constant coeffcients to make each integral to 1
  K1 <- 1 / sum((atan(w1 / whalf) + pi / 2))
  K2 <- 1 / sum((atan(w2 / whalf) + pi / 2))
  ## Part1: (w1)T (w2)
  Inter_sum <- - 2 * integral(w1, w2, whalf) * K1 * K2 * whalf^2
  ## Part2: (w1)T (w1)
  Quadra_sum_1 <- integral(w1, w1, whalf) * K1^2 * whalf^2
  ## Part3: (w2)T (w2)
  Quadra_sum_2 <- integral(w2, w2, whalf) * K2^2 * whalf^2
  return(sqrt(max(Quadra_sum_1 + Quadra_sum_2 + Inter_sum, 0)))
}


## QIM Matrix calculation 
## Input: SampleGen (list of N by N adjacency matrix of networks, 
##                   the first few belong to one group while rest the other group), 
##        nsample1 (number, number of matrix belong to the first group), 
##        l (number, weight to be put on global distance)
## Output: QIMMatrix (N by N distance matrix)
QIMMatrix_Calc <- function(tempdata, N, n1, whalf, lam_Mat, l=1){
  n2  <- N - n1
  QIMMatrix <- matrix(0, N, N) 
  for (i in 1:(N-1)){
    #print(i)
    tempmat1 <- tempdata[[i]]
    for (j in (i+1):N){
      #print(j)
      tempmat2 <- tempdata[[j]]
      QIMMatrix[i, j] <- HamD_Calc(tempmat1, tempmat2) + IpsD_Calc(whalf, lam_Mat, i, j)
      ## symmetric matrix
      QIMMatrix[j, i] <- QIMMatrix[i, j]
    }
  }
  tempmat1 <- tempmat2 <- c()
  return(QIMMatrix)
}

## HM Matrix calculation (No Global Information)
## Output: HMMatrix (N by N distance matrix)
HMMatrix_Calc <- function(tempdata, N, n1){
  n2  <- N - n1
  HMMatrix <- matrix(0, N, N) 
  for (i in 1:(N-1)){
    #print(i)
    tempmat1 <- tempdata[[i]]
    for (j in (i+1):N){
      #print(j)
      tempmat2 <- tempdata[[j]]
      HMMatrix[i, j] <- HamD_Calc(tempmat1, tempmat2)
      ## symmetric matrix
      HMMatrix[j, i] <- HMMatrix[i, j]
    }
  }
  tempmat1 <- tempmat2 <- c()
  return(HMMatrix)
}

## MR Matrix calculation (Only need call this to calculate final proximity(probability) matrix)
## Input: DIMMatrix (N by N matrix with element represent QIM distance between each pair)
## Output: MRMatrix (N by N matrix with element represent remoteness between each pair) 
cppFunction('NumericMatrix MR_MAT(NumericMatrix QIM){
              int N = QIM.nrow();
              NumericMatrix MR(N, N);

              for(int i = 0; i < N; i++){
                for(int j = i + 1; j < N; j++){
                  double counts = 0;
                  for(int k = 0; k < N; k++){
                    if(((QIM(k, i) > QIM(i, j)) && (QIM(k, j) > QIM(i, j)) && k!=i && k!=j)){
                      counts += 1;
                    }
                  }
                  MR(i, j) = 1 - (counts / (N - 2));
                  MR(j, i) = 1 - (counts / (N - 2));
                }
              }
              return (MR);
            }')

## F_MR test statistics calculation (Called by pvalue_Calc)
## Input: MRMatrix (N by N probability matrix), 
##        v1 (Vector of number, indicators of subjects that belong to the first group), 
##        v2 (Vector of number, indicators of subjects that belong to the other group)
## Output: F_MR (number, test statistics of our test)
cppFunction('NumericVector PCalc(NumericMatrix DisMat, int n1, int nperm){
              srand(1);
              int N = DisMat.nrow();
              int n2 = N - n1;
              NumericVector pperm(nperm + 1);
              NumericVector idlist(N);
              std::iota (std::begin(idlist), std::end(idlist), 0);

              for(int k = 0; k < (nperm+1); k++){
                double sum1 = 0;
                double sum2 = 0;
                double sum12 = 0;
                for(int i = 0; i < n1; i++){
                  for(int j = i + 1; j < n1; j++){
                    sum1 += pow(DisMat(idlist[i], idlist[j]), 2); 
                  }
                }
                for(int i = n1; i < N; i++){
                  for(int j = i + 1; j < N; j++){
                    sum2 += pow(DisMat(idlist[i], idlist[j]), 2); 
                  }
                }
                for(int i = 0; i < n1; i++){
                  for(int j = n1; j < N; j++){
                    sum12 += pow(DisMat(idlist[i], idlist[j]), 2); 
                  }
                }
                pperm[k] = ((sum12 + sum1 + sum2)/N - sum1 / n1 + sum2 / n2) / (sum1 / n1 + sum2 / n2) * (N - 2);
                for (int idx = 0; idx < N - 1; idx++){
                  int idy = idx + rand() % (N - idx);
                  std::swap(idlist[idx], idlist[idy]);
                }
              }
            return pperm;
            }')


############################################
## Main Function
############################################
## Input: WorkSample List of N by N network matrix for both groups, 
##        the first N_Sample1 networks from the first group 
##        the later N_Sample2 networks from the second group
##        N_Sample1  number of subjects belong to group 1, 
##        N_Sample2  number of subjects belong to group 2
##        n_perm number of permutation for test; delta continuous correction; seed random seed
## Output: List of Test p-values and pairwise distance matrix 
##        p_QIM , Sample_QIM, Original NEPTUNE test p-value and pairwise distance based on QIM
##        p_MR  , Sample_MR,  High-dimensional NEPTUNE test p-value and pairwise distance based on MR
##        Diagnosis plot for hubness will be plot. If too many points overlapped then hubness issue exists. 
NEPTUNE_Test <- function(WorkSample, N_Sample1, N_Sample2, n_perm = 1000, delta = 0.5, seed=11) {
  set.seed(seed)
  workindex = c(rep(1, N_Sample1), rep(2, N_Sample2))
  
  N_Total <- length(WorkSample)
  Vsize <- nrow(WorkSample[[1]])                       # Vertex size
  Whalf <- Half_Width(Vsize)                           # Compute half-width
  Eigen_Mat <- matrix(0, nrow = Vsize, ncol = N_Total) # Eigenvalue matrix
  
  # Calculate Laplacian eigenvectors
  for (idx in 1:N_Total) {
    Eigen_Mat[, idx] <- Lap_Eigen(WorkSample[[idx]])
  }
  
  # QIM Distance
  Sample_QIM <- QIMMatrix_Calc(WorkSample, N_Total, N_Sample1, Whalf, Eigen_Mat, l = 1)
  temp_qim <- PCalc(Sample_QIM, N_Sample1, n_perm)
  p_QIM <- (sum(temp_qim > temp_qim[1]) + delta) / (n_perm + delta)
  
  # High-dimensional (MR) Distance
  Sample_MR <- MR_MAT(Sample_QIM)
  temp_mr <- PCalc(Sample_MR, N_Sample1, n_perm)
  p_MR <- (sum(temp_mr > temp_mr[1]) + delta) / (n_perm + delta)
  
  # 2D Illustrations
  result_mstknn_qim <- mst.knn(Sample_QIM)
  plot(result_mstknn_qim$network, vertex.size = 18,
       vertex.color = (workindex + 3),
       layout = igraph::layout.fruchterman.reingold(result_mstknn_qim$network, niter = 10000),
       main = "MST-kNN Clustering Based on QIM \n Colored by Phenotype",
       vertex.label = NA)
  
  # Return p-values
  return(list(p_QIM = p_QIM, p_MR = p_MR,
              Sample_QIM = Sample_QIM, Sample_MR = Sample_MR))
}