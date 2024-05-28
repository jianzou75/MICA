#' Mouse metabolism data list
#'
#' Bulk expression profiles are measured in mice with three genotypes
#' (wild-type, LCAD knock-out, and VLCAD knock-out).
#' LCAD deficiency is associated with impaired fatty acid oxidation, and VLCAD deficiency is associated with energy metabolism disorders in children. Microarray experiments were conducted on tissues from 12 mice (four mice per genotype) including brown fat, liver, heart, and skeletal
#'
#' @format A list of data frames with rows for the samples and columns for the genes (Unigene ID).
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
