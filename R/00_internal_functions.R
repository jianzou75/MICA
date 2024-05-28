library(parallel)
library(igraph)

# Multi-class Correlation (MCC)
## mcc function works for one single biomarker (gene)
## study1.data is a vector for one gene expression
mcc <- function(study1.data, study1.label, study2.data, study2.label){
  if(length(unique(study1.label)) != length(unique(study2.label))){
    stop("The number of classes should be equal between two studies.")
  }

  w.est <- 1/(length(unique(study1.label)))
  study1.mu.grp  <- aggregate(study1.data, list(study1.label), FUN=mean)[,-1]
  study1.var.grp <- aggregate(study1.data, list(study1.label), FUN=var)[,-1]
  study2.mu.grp  <- aggregate(study2.data, list(study2.label), FUN=mean)[,-1]
  study2.var.grp <- aggregate(study2.data, list(study2.label), FUN=var)[,-1]

  study1.mu  <- mean(study1.mu.grp)
  study1.var <- w.est * sum(study1.mu.grp^2 + study1.var.grp) - (study1.mu^2)
  study2.mu  <- mean(study2.mu.grp)
  study2.var <- w.est * sum(study2.mu.grp^2 + study2.var.grp) - (study2.mu^2)

  stat <- ifelse(any(c(study1.var, study2.var) <= 0), 0,
                 ((w.est * sum(study1.mu.grp * study2.mu.grp)) - (study1.mu * study2.mu))/(sqrt(study1.var * study2.var)))
  return(mcc = stat)
}

## mcc.matr function works for gene expression matrix (gene X sample)
## study1.data.matrix is a matrix for gene expression
mcc.matr <- function(study1.data.matrix, study1.label, study2.data.matrix, study2.label){
  n.feature <- nrow(study1.data.matrix)
  stat <- sapply(1:n.feature, function(i) mcc(study1.data.matrix[i, ], study1.label,
                                              study2.data.matrix[i, ], study2.label))
  return(mcc = stat)
}


# Minimum multi-class Correlation (min-MCC)
## min.mcc function works for one single biomarker (gene)
min.mcc <- function(study.data.list, study.label.list){
  n.study <- length(study.data.list)
  cmb <- expand.grid(1:n.study, 1:n.study)
  cmb <- cmb[cmb$Var1 < cmb$Var2,]
  mcc_tbl <- c()

  for (c in 1:nrow(cmb)){
    mcc_tbl <- c(mcc_tbl, mcc(study1.data = study.data.list[[cmb[c,1]]], study1.label = study.label.list[[cmb[c,1]]],
                              study2.data = study.data.list[[cmb[c,2]]], study2.label = study.label.list[[cmb[c,2]]]))
  }
  stat <- min(mcc_tbl)
  return(min.mcc = stat)
}

## min.mcc.matr function works for gene expression matrix (gene X sample)
## study1.data.matrix.list is a list of matrices for gene expression
min.mcc.matr <- function(study.data.matrix.list, study.label.list){
  n.study <- length(study.data.matrix.list)
  n.feature <- nrow(study.data.matrix.list[[1]])
  cmb <- expand.grid(1:n.study, 1:n.study)
  cmb <- cmb[cmb$Var1 < cmb$Var2,]
  mcc_tbl <- data.frame(matrix(nrow = nrow(cmb), ncol = n.feature))

  for (c in 1:nrow(cmb)){
    mcc_tbl[c, ] <- mcc.matr(study1.data = study.data.matrix.list[[cmb[c,1]]], study1.label = study.label.list[[cmb[c,1]]],
                             study2.data = study.data.matrix.list[[cmb[c,2]]], study2.label = study.label.list[[cmb[c,2]]])
  }
  stat <- apply(mcc_tbl, 2, min)
  return(min.mcc = stat)
}


# mutual information (mi+)
## mi function works for one single biomarker (gene)
## study1.data is a vector for one gene expression
mi <- function(study1.data, study1.label, study2.data, study2.label){
  mcc.value <- mcc(study1.data, study1.label, study2.data, study2.label)
  stat <- ifelse(mcc.value < 0, 0, -log(1 - mcc.value^2)/2)
  return(mi = stat)
}

## mi.matr function works for gene expression matrix (gene X sample)
## study1.data.matrix is a matrix for gene expression
mi.matr <- function(study1.data.matrix, study1.label, study2.data.matrix, study2.label){
  n.feature <- nrow(study1.data.matrix)
  stat <- sapply(1:n.feature, function(i) mi(study1.data.matrix[i, ], study1.label,
                                             study2.data.matrix[i, ], study2.label))
  return(mi = stat)
}

# Generalized mutual information (gmi+)
## gmi function works for one single biomarker (gene)
gmi <- function(study.data.list, study.label.list){
  n.study <- length(study.data.list)
  cmb <- expand.grid(1:n.study, 1:n.study)
  cmb <- cmb[cmb$Var1 < cmb$Var2,]
  correlation <- c()

  ### Calculate the off-diagonal elements
  for (c in 1:nrow(cmb)){
    correlation <- c(correlation, mcc(study1.data = study.data.list[[cmb[c,1]]], study1.label = study.label.list[[cmb[c,1]]],
                                      study2.data = study.data.list[[cmb[c,2]]], study2.label = study.label.list[[cmb[c,2]]]))
  }
  correlation[correlation < 0] = 0 ### A trunctated term

  ### Calculate the diagonal elements
  study.var.grp <- lapply(1:n.study, function(s) aggregate(study.data.list[[s]], list(study.label.list[[s]]), FUN=var)[,-1])
  study.mu.grp  <- lapply(1:n.study, function(s) aggregate(study.data.list[[s]], list(study.label.list[[s]]), FUN=mean)[,-1])
  study.var.tbl <- c()
  w.est <- 1/(length(unique(study.label.list[[1]])))
  for(s in 1:n.study){
    study.mu <-  mean(study.mu.grp[[s]])
    study.var <- w.est * sum(study.mu.grp[[s]]^2 + study.var.grp[[s]]) - (study.mu^2)
    study.var[study.var <= 0] <- 0
    study.var.tbl <- c(study.var.tbl, study.var)
  }

  ### Obtain the covariance matrix
  Sigma <- matrix(nrow = n.study, ncol = n.study)
  diag(Sigma) <- study.var.tbl
  for (c in 1:nrow(cmb)){
    off_dig <- sqrt(study.var.tbl[cmb[c,1]] * study.var.tbl[cmb[c,2]]) * correlation[c]
    Sigma[cmb[c,1], cmb[c,2]] = off_dig
    Sigma[cmb[c,2], cmb[c,1]] = off_dig
  }

  ### Calculate generalized MI
  study_rm <- apply(Sigma, 1, function(x) all(x == 0)) & apply(Sigma, 2, function(x) all(x == 0))
  if(length(study_rm) == 2 & any(study_rm)){
    stat <- 0
  } else{
    if(any(study_rm)){
      Sigma <- Sigma[-which(study_rm), -which(study_rm)]
    }
    stat <- -log(det(Sigma)/det(diag(diag(Sigma))))/2
  }
  return(gmi = stat)
}

## gmi.matr function works for gene expression matrix (gene X sample)
## study.data.matrix.list is a list of matrices for gene expression
gmi.matr <- function(study.data.matrix.list, study.label.list){
  n.feature <- nrow(study.data.matrix.list[[1]])
  stat <- sapply(1:n.feature, function(s) gmi(lapply(study.data.matrix.list, function(x) x[s,]),
                                              study.label.list))
  return(gmi = stat)
}


# mutual information concordance analysis (MICA)
## mica works function works for one single biomarker (gene)
## the min-mcc statistics is also included
mica <- function(study.data.list, study.label.list, n.perm = 500, p.threshold = 0.05, n.parallel = 40, post.hoc = TRUE){

  ### generate the permuted labels
  n.study   <- length(study.data.list)
  perm.label <- list()
  for(p in 1:n.perm){
    perm <- lapply(1:n.study, function(s) sample(study.label.list[[s]]))
    perm.label <- c(perm.label, list(perm))
  }

  ### gmi and min-mcc (global-wide pre-selection)
  perm.gmi <- mclapply(1:n.perm, function(p) gmi(study.data.list, perm.label[[p]]), mc.cores = n.parallel)
  perm.minmcc <- mclapply(1:n.perm, function(p) min.mcc(study.data.list, perm.label[[p]]), mc.cores = n.parallel)
  stat.gmi <- gmi(study.data.list, study.label.list)
  stat.minmcc <- min.mcc(study.data.list, study.label.list)
  pval.gmi <- mean(c(perm.gmi) >= stat.gmi)
  pval.minmcc <- mean(c(perm.minmcc) >= stat.minmcc)

  ### post-hoc investigation
  cmb <- NULL
  if(post.hoc == TRUE & pval.gmi < p.threshold){
    cmb <- expand.grid(1:n.study, 1:n.study)
    cmb <- cmb[cmb$Var1 < cmb$Var2,]
    stat.mi.list <- c()
    pval.mi.list <- c()
    for (c in 1:nrow(cmb)){
      stat <- mi(study1.data = study.data.list[[cmb[c,1]]], study1.label = study.label.list[[cmb[c,1]]],
                 study2.data = study.data.list[[cmb[c,2]]], study2.label = study.label.list[[cmb[c,2]]])
      perm.stat <- mclapply(1:n.perm,
                            function(p) mi(study1.data = study.data.list[[cmb[c,1]]], study1.label = perm.label[[p]][[cmb[c,1]]],
                                           study2.data = study.data.list[[cmb[c,2]]], study2.label = perm.label[[p]][[cmb[c,2]]]),
                            mc.cores = n.parallel)
      pval <- mean(c(perm.stat) >= stat)
      stat.mi.list <- c(stat.mi.list, stat)
      pval.mi.list <- c(pval.mi.list, pval)
    }
    cmb$pval <- pval.mi.list
    cmb$stat <- stat.mi.list
  }

  return(list(gmi.stat = stat.gmi, gmi.pval = pval.gmi,
              minmcc.stat = stat.minmcc, minmcc.pval = pval.minmcc,
              pairwise.tbl = cmb))
}

