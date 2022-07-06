#' Multi-class correlation (MCC)
#'
#' This function is used to calculate the multi-class correlation (MCC) between two studies of gene expression, where the gene expressions are assumed to follow Gaussian distribution.
#'
#' @param study1.data a data frame or matrix with rows for the genes and columns for the samples. (study 1)
#' @param study1.label a vector of the class labels for study 1.
#' @param study2.data a data frame or matrix with rows for the genes and columns for the samples. (study 2)
#' @param study2.label a vector of the class labels for study 1.
#' @param w.est a vector of pre-defined class weights that add up to 1. If not specified, it will be calculated from the data.
#'
#' @return multi-class correlation (MCC) between study 1 and study 2.
#' @export
#'
#' @references
#' Lu, Shuya, Jia Li, Chi Song, Kui Shen, and George C. Tseng. 2010. “Biomarker Detection in the Integration of Multiple Multi-Class Genomic Studies.” Bioinformatics  26 (3): 333–40. https://doi.org/10.1093/bioinformatics/btp669.
#'
#' @examples
MCC <- function(study1.data, study1.label, study2.data, study2.label, w.est = NULL){
  if(ncol(study1.data) != ncol(study2.data)){stop("The number of features should be equal between two studies.")}
  n.feature <- ncol(study1.data)

  if(is.null(w.est)){w.est <- unname(c(table(c(study1.label, study2.label)))/length(c(study1.label, study2.label)))}
  study1.mu.grp  <- aggregate(study1.data, list(study1.label), FUN=mean)[,-1]
  study1.var.grp <- aggregate(study1.data, list(study1.label), FUN=var)[,-1]
  study2.mu.grp  <- aggregate(study2.data, list(study2.label), FUN=mean)[,-1]
  study2.var.grp <- aggregate(study2.data, list(study2.label), FUN=var)[,-1]

  study1.mu  <- apply(study1.mu.grp, 2, function(x) sum(x * w.est))
  study1.var <- sapply(1:n.feature, function(i) sum(w.est*(study1.mu.grp[,i]^2 + study1.var.grp[,i])) - (study1.mu[i]^2))
  study2.mu  <- apply(study2.mu.grp, 2, function(x) sum(x * w.est))
  study2.var <- sapply(1:n.feature, function(i) sum(w.est*(study2.mu.grp[,i]^2 + study2.var.grp[,i])) - (study2.mu[i]^2))

  MCC_stat <- sapply(1:n.feature, function(i) (sum(w.est * study1.mu.grp[,i] * study2.mu.grp[,i]) - (study1.mu[i] * study2.mu[i]))/(sqrt(study1.var[i] * study2.var[i])))

  return(MCC_stat = MCC_stat)
}


#' Minimum multi-class correlation (min-MCC)
#'
#' This function is used to calculate the minimum multi-class correlation (min-MCC) among multiple studies of gene expression, where the gene expressions are assumed to follow Gaussian distribution.
#'
#' @param study.data.list a list of data frames or matrices with rows for the genes and columns for the samples.
#' @param study.label.list a list of vectors of the class labels.
#' @param w.est a vector of pre-defined class weights that add up to 1. If not specified, it will be calculated from the data.
#'
#' @return minimum multi-class correlation (min-MCC) for multiple studies.
#' @export
#' @references
#' Lu, Shuya, Jia Li, Chi Song, Kui Shen, and George C. Tseng. 2010. “Biomarker Detection in the Integration of Multiple Multi-Class Genomic Studies.” Bioinformatics  26 (3): 333–40. https://doi.org/10.1093/bioinformatics/btp669.
#'
#' @examples
minMCC <- function(study.data.list, study.label.list, w.est = NULL){
  if(!all(unlist(lapply(study.data.list, ncol)) == ncol(study.data.list[[1]]))){stop("The number of features should be equal between two studies.")}
  if(is.null(w.est)){w.est <- unname(c(table(unlist(study.label.list)))/length(unlist(study.label.list)))}

  n.feature <- ncol(study.data.list[[1]])
  n.study <- length(study.data.list)
  cmb <- expand.grid(1:n.study, 1:n.study)
  cmb <- cmb[cmb$Var1 < cmb$Var2,]
  MCC_tbl <- data.frame(matrix(nrow = nrow(cmb), ncol = n.feature))

  for (c in 1:nrow(cmb)){
    MCC_tbl[c, ] <- MCC(study1.data = study.data.list[[cmb[c,1]]], study1.label = study.label.list[[cmb[c,1]]],
                        study2.data = study.data.list[[cmb[c,2]]], study2.label = study.label.list[[cmb[c,2]]],
                        w.est = w.est)
  }

  minMCC_stat <- apply(MCC_tbl, 2, min)
  return(minMCC_stat = minMCC_stat)
}


#' Multi-class total correlation (MC-TC)
#'
#' This function is used to calculate the multi-class total correlation (MC-TC) among multiple studies of gene expression, where the gene expressions are assumed to follow Gaussian distribution.
#'
#' @param study.data.list a list of data frames or matrices with rows for the genes and columns for the samples.
#' @param study.label.list a list of vectors of the class labels.
#' @param w.est a vector of pre-defined class weights that add up to 1. If not specified, it will be calculated from the data.
#'
#' @return Multi-class total correlation (MC-TC) for multiple studies.
#' @export
#' @references
#' Lu, Shuya, Jia Li, Chi Song, Kui Shen, and George C. Tseng. 2010. “Biomarker Detection in the Integration of Multiple Multi-Class Genomic Studies.” Bioinformatics  26 (3): 333–40. https://doi.org/10.1093/bioinformatics/btp669.
#'
#' @examples
MCTC <- function(study.data.list, study.label.list, w.est = NULL){
  if(!all(unlist(lapply(study.data.list, ncol)) == ncol(study.data.list[[1]]))){stop("The number of features should be equal between two studies.")}
  if(is.null(w.est)){w.est <- unname(c(table(unlist(study.label.list)))/length(unlist(study.label.list)))}

  ## Obtain the between study correlation
  n.feature <- ncol(study.data.list[[1]])
  n.study <- length(study.data.list)
  cmb <- expand.grid(1:n.study, 1:n.study)
  cmb <- cmb[cmb$Var1 < cmb$Var2,]
  MCC_tbl <- data.frame(matrix(nrow = nrow(cmb), ncol = n.feature))

  for (c in 1:nrow(cmb)){
    MCC_tbl[c, ] <- MCC(study1.data = study.data.list[[cmb[c,1]]], study1.label = study.label.list[[cmb[c,1]]],
                        study2.data = study.data.list[[cmb[c,2]]], study2.label = study.label.list[[cmb[c,2]]],
                        w.est = w.est)
  }
  MCC_tbl[MCC_tbl < 0] = 0 ## A trunctated term

  ## Obtain the variance for each study
  study.var.grp <- lapply(1:n.study, function(s) aggregate(study.data.list[[s]], list(study.label.list[[s]]), FUN=var)[,-1])
  study.mu.grp  <- lapply(1:n.study, function(s) aggregate(study.data.list[[s]], list(study.label.list[[s]]), FUN=mean)[,-1])
  study.var.tbl <- data.frame(matrix(nrow = n.study, ncol = n.feature))

  for(s in 1:n.study){
    study.mu <-  apply(study.mu.grp[[s]], 2, function(x) sum(x * w.est))
    study.var <- sapply(1:n.feature, function(i) sum(w.est*(study.mu.grp[[s]][,i]^2 + study.var.grp[[s]][,i])) - (study.mu[i]^2))
    study.var.tbl[s,] <- study.var
  }

  ## Obtain the total correlation (TC) and covariance matrix for each feature
  TC_features = c()
  for(f in 1:n.feature){
    Sigma <- matrix(nrow = n.study, ncol = n.study)
    diag(Sigma) <- study.var.tbl[,f]

    for (c in 1:nrow(cmb)){
      off_dig <- sqrt(study.var.tbl[cmb[c,1],f] * study.var.tbl[cmb[c,2],f]) * MCC_tbl[c,f]
      Sigma[cmb[c,1], cmb[c,2]] = off_dig
      Sigma[cmb[c,2], cmb[c,1]] = off_dig
    }

    ### total correlation
    TC <- -log(det(Sigma)/det(diag(diag(Sigma))))/2
    TC_features <- c(TC_features, TC)
  }

  return(MCTC_stat = TC_features)
}


#' A function to combine the MCTC and minMCC
#'
#' @keywords internal
#' @param study.data.list a list of data frames or matrices with rows for the genes and columns for the samples.
#' @param study.label.list a list of vectors of the class labels.
#' @param w.est a vector of pre-defined class weights that add up to 1. If not specified, it will be calculated from the data.
#'
#' @return
#' @export
#'
#' @examples
mc.stat <- function(study.data.list, study.label.list, w.est = NULL){
  if(is.null(w.est)){w.est <- unname(c(table(unlist(study.label.list)))/length(unlist(study.label.list)))}

  MCTC_stat   = MCTC(study.data.list, study.label.list, w.est = w.est)
  minMCC_stat = minMCC(study.data.list, study.label.list, w.est = w.est)

  return(list(MCTC_stat = MCTC_stat, minMCC_stat = minMCC_stat))
}


#' Permutation test for min-MCC, MCTC, and post-hoc MCMI analysis
#'
#' @param study.data.list a list of data frames or matrices with rows for the genes and columns for the samples.
#' @param study.label.list a list of vectors of the class labels.
#' @param w.est a vector of pre-defined class weights that add up to 1. If not specified, it will be calculated from the data.
#' @param n.perm permutation times
#' @param n.parallel number of cores used for parallel computation
#'
#' @return
#' @export
#'
#' @examples
mc.permute.test <- function(study.data.list, study.label.list, w.est = NULL, n.perm = 1000, n.parallel = 50){
  if(!all(unlist(lapply(study.data.list, ncol)) == ncol(study.data.list[[1]]))){stop("The number of features should be equal between two studies.")}
  if(is.null(w.est)){w.est <- unname(c(table(unlist(study.label.list)))/length(unlist(study.label.list)))}
  feature.names <- colnames(study.data.list[[1]])


  ## Obtain the permuted labels
  n.feature <- ncol(study.data.list[[1]])
  n.study   <- length(study.data.list)
  size.study <- sapply(1:n.study, function(s) length(study.label.list[[s]]))
  perm.label <- list()
  for(p in 1:n.perm){
    perm <- lapply(1:n.study, function(s) sample(study.label.list[[s]]))
    perm.label <- c(perm.label, list(perm))
  }

  ### an overall analysis
  ## permutation stat
  perm.res <- mclapply(1:n.perm, function(p) mc.stat(study.data.list, perm.label[[p]], w.est = w.est), mc.cores = n.parallel)
  MCTC_perm <- t(sapply(1:n.perm, function(p) perm.res[[p]]$MCTC_stat))
  minMCC_perm <- t(sapply(1:n.perm, function(p) perm.res[[p]]$minMCC_stat))
  stat <- mc.stat(study.data.list, study.label.list, w.est = w.est)


  ## TC table
  MCTC_stat <- stat$MCTC_stat
  MCTC_pval <- sapply(1:n.feature, function(s) mean(MCTC_perm[,s] >= stat$MCTC_stat[s]))
  MCTC_qval <- p.adjust(MCTC_pval, method = "fdr")
  MCTC_tbl <- data.frame(stat = MCTC_stat,
                         pval = MCTC_pval,
                         qval = MCTC_qval)
  rownames(MCTC_tbl) <- feature.names

  ## minMCC table
  minMCC_stat <- stat$minMCC_stat
  minMCC_pval <- sapply(1:n.feature, function(s) mean(minMCC_perm[,s] >= stat$minMCC_stat[s]))
  minMCC_qval <- p.adjust(minMCC_pval, method = "fdr")
  minMCC_tbl <- data.frame(stat = minMCC_stat,
                           pval = minMCC_pval,
                           qval = minMCC_qval)
  rownames(minMCC_tbl) <- feature.names


  ### post-hoc pair-wise analysis
  pairs_MI_tbl_list <- list()
  allpairs = combn(1:length(study.data.list), 2)
  for(i in 1:ncol(allpairs)){
    pairs_perm <- mclapply(1:n.perm, function(p) mc.stat(study.data.list[allpairs[,i]], perm.label[[p]][allpairs[,i]], w.est = w.est), mc.cores = n.parallel)
    pairs_perm <- t(sapply(1:n.perm, function(p) pairs_perm[[p]]$MCTC_stat))

    pairs_stat <- mc.stat(study.data.list[allpairs[,i]], study.label.list[allpairs[,i]], w.est = w.est)
    pairs_MI_stat <- pairs_stat$MCTC_stat
    pairs_MI_pval <- sapply(1:n.feature, function(s) mean(pairs_perm[,s] >= pairs_stat$MCTC_stat[s]))
    pairs_MI_qval <- p.adjust(pairs_MI_pval, method = "fdr")

    pairs_MI_tbl <- data.frame(stat = pairs_MI_stat,
                               pval = pairs_MI_pval,
                               qval = pairs_MI_qval)
    rownames(pairs_MI_tbl) <- feature.names

    pairs_MI_tbl_list <- c(pairs_MI_tbl_list, list(pairs_MI_tbl))
  }
  names(pairs_MI_tbl_list) <- paste(combn(1:length(study.data.list), 2)[1,], combn(1:length(study.data.list), 2)[2,], sep = ":")

  return(list(MCTC_tbl = MCTC_tbl,
              minMCC_tbl = minMCC_tbl,
              pairs_MI_tbl_list = pairs_MI_tbl_list))
}
