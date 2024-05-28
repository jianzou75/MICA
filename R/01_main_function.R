#' Full function for Mutual Information Concordance Analysis (MICA)
#'
#' `mica.full()` is the main function for MICA, which runs the MICA and min-MCC simultaneously.
#' MICA first detects biomarkers with concordant multi-class patterns across partial or all of the omics studies using a global test by mutual information.
#' A post hoc analysis is then performed for each detected biomarkers and identify studies with concordant pattern.
#'
#'
#' @param study.data.matrix.list list of matrices for gene expression with genes in rows and samples in columns
#' @param study.label.list list of vectors for sample labels
#' @param n.perm number of permutations
#' @param p.threshold p-value threshold for post-hoc analysis
#' @param n.parallel number of parallel cores
#' @param post.hoc logical for post-hoc analysis
#'
#' @import parallel
#'
#' @return
#' A list with the following components:
#'  \item{gmi.stat}{statistics for generalized one-sided corrected mutual information (gMI+)}
#'  \item{minmcc.stat}{statistics for min-MCC}
#'  \item{pval.gmi}{p-values for gMI+}
#'  \item{pval.minmcc}{p-values for min-mcc}
#'  \item{pairwise.tbl}{table for post-hoc pairwise comparsion}
#'
#' @export
#'
#' @examples
mica.full <- function(study.data.matrix.list , study.label.list,
                      n.perm = 500, p.threshold = 0.05, n.parallel = 1, post.hoc = TRUE){

  ### generate the permuted labels
  n.study   <- length(study.data.matrix.list)
  n.feature <- nrow(study.data.matrix.list[[1]])
  perm.label <- list()
  for(p in 1:n.perm){
    perm <- lapply(1:n.study, function(s) sample(study.label.list[[s]]))
    perm.label <- c(perm.label, list(perm))
  }

  ### gmi and min-mcc (global-wide pre-selection)
  perm.gmi <- mclapply(1:n.perm, function(p) gmi.matr(study.data.matrix.list, perm.label[[p]]), mc.cores = n.parallel)
  perm.minmcc <- mclapply(1:n.perm, function(p) min.mcc.matr(study.data.matrix.list, perm.label[[p]]), mc.cores = n.parallel)
  stat.gmi <- gmi.matr(study.data.matrix.list, study.label.list)
  stat.minmcc <- min.mcc.matr(study.data.matrix.list, study.label.list)
  pval.gmi <- sapply(1:n.feature, function(s) mean(unlist(perm.gmi) >= stat.gmi[s]))
  pval.minmcc <- sapply(1:n.feature, function(s) mean(unlist(perm.minmcc) >= stat.minmcc[s]))
  qval.gmi <- p.adjust(pval.gmi, method = "BH")
  qval.minmcc <- p.adjust(pval.minmcc, method = "BH")

  ### post-hoc investigation
  cmb <- NULL; pairwise.tbl <- NULL
  preselect.genes <- which(pval.gmi < p.threshold)
  if(post.hoc == TRUE){
    study.data.matrix.list.filter <- lapply(study.data.matrix.list, function(l) l[pval.gmi < p.threshold, ])
    cmb <- expand.grid(1:n.study, 1:n.study)
    cmb <- cmb[cmb$Var1 < cmb$Var2,]
    stat.mi.list <- matrix(nrow = nrow(cmb), ncol = length(preselect.genes))
    pval.mi.list <- matrix(nrow = nrow(cmb), ncol = length(preselect.genes))
    for (c in 1:nrow(cmb)){
      stat <- mi.matr(study1.data = study.data.matrix.list.filter[[cmb[c,1]]], study1.label = study.label.list[[cmb[c,1]]],
                      study2.data = study.data.matrix.list.filter[[cmb[c,2]]], study2.label = study.label.list[[cmb[c,2]]])
      perm.stat <- mclapply(1:n.perm,
                            function(p) mi.matr(study1.data = study.data.matrix.list.filter[[cmb[c,1]]], study1.label = perm.label[[p]][[cmb[c,1]]],
                                                study2.data = study.data.matrix.list.filter[[cmb[c,2]]], study2.label = perm.label[[p]][[cmb[c,2]]]),
                            mc.cores = n.parallel)
      pval <- sapply(1:length(preselect.genes),
                     function(g) mean(sapply(perm.stat, function(x) x[g]) >= stat[g]))
      stat.mi.list[c, ] <- stat
      pval.mi.list[c, ] <- pval
    }
    pairwise.tbl <- list()
    for (g in 1:length(preselect.genes)){
      cmb$pval <- pval.mi.list[,g]
      cmb$stat <- stat.mi.list[,g]
      pairwise.tbl <- c(pairwise.tbl, list(cmb))
    }
  }

  return(list(gmi.stat = stat.gmi, gmi.pval = pval.gmi,
              minmcc.stat = stat.minmcc, minmcc.pval = pval.minmcc,
              pairwise.tbl = pairwise.tbl))
}
