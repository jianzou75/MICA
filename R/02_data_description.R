#' Mouse metabolism data list
#'
#' Bulk expression profiles are measured in mice with three genotypes
#' (wild-type, LCAD knock-out, and VLCAD knock-out).
#' LCAD deficiency is associated with impaired fatty acid oxidation, and VLCAD deficiency is associated with energy metabolism disorders in children. Microarray experiments were conducted on tissues from 12 mice (four mice per genotype) including brown fat, liver, heart, and skeletal
#'
#' @format A list of data frames with rows for the genes (Unigene ID) and columns for the samples.
#'
#' @references
#' Lu, S., Li, J., Song, C., Shen, K., & Tseng, G. C. (2010). Biomarker detection in the integration of multiple multi-class genomic studies. Bioinformatics, 26(3), 333-340.
#'
"metabolism_data"


#' Mouse metabolism label list
#'
#' A list of vectors for the labels of \link{metabolism_data}.
#'
#' @format A list of vectors (wt, VLCAD, LCAD) for the genotypes of the samples.
#'
#' @references
#' Lu, S., Li, J., Song, C., Shen, K., & Tseng, G. C. (2010). Biomarker detection in the integration of multiple multi-class genomic studies. Bioinformatics, 26(3), 333-340.
#'
"metabolism_label"


#' Estrogene data list
#'
#' The dataset comes from the EstroGene project, which integrates estrogen-related datasets, including RNA-Seq and microarray data.
#' This dataset focuses on gene expression studies involving estrogen receptor-positive (ER+) samples treated with estradiol (E2) at doses greater than 1nM. The samples were categorized by cell line, sequencing technology, and treatment duration (short, medium, long).
#' Data normalization was performed using TMM and ComBat methods. The resulting pooled studies include MCF7 microarray, MCF7 RNA-Seq, and T47D RNA-Seq
#'
#' @format A list of data frames with rows for the genes and columns for the samples.
#'
#' @references
#' Li, Z., Li, T., Yates, M. E., Wu, Y., Ferber, A., Chen, L., ... & Lee, A. V. (2023). The EstroGene database reveals diverse temporal, context-dependent, and bidirectional estrogen receptor regulomes in breast cancer. Cancer Research, 83(16), 2656-2674.
"estrogene_data"


#' Estrogene label list
#'
#' A list of vectors for the labels of \link{estrogene_data}.
#'
#' @format A list of vectors (short, medium, long) for the treatment duration of the samples.
#'
#' @references
#' Li, Z., Li, T., Yates, M. E., Wu, Y., Ferber, A., Chen, L., ... & Lee, A. V. (2023). The EstroGene database reveals diverse temporal, context-dependent, and bidirectional estrogen receptor regulomes in breast cancer. Cancer Research, 83(16), 2656-2674.
#'
"estrogene_label"


