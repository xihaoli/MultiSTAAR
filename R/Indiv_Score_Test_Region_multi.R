#' Multi-trait score test for individual variants in a given variant-set
#'
#' The \code{Indiv_Score_Test_Region_multi} function takes in genotype and the object from fitting the null
#' model to analyze the associations between multiple (quantitative) phenotypes and
#' all individual variants in a given variant-set by using score test.
#' @param genotype an n*p genotype matrix (dosage matrix) of the target sequence,
#' where n is the sample size and p is the number of genetic variants.
#' @param obj_nullmodel an object from fitting the null model, which is the
#' output from either \code{\link{fit_null_glmmkin_multi}} function for unrelated or
#' related samples. Note that \code{\link{fit_null_glmmkin_multi}}
#' is a wrapper of the \code{\link{glmmkin}} function from the \code{\link{GMMAT}} package.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @return a data frame with p rows corresponding to the p genetic variants in the given variant-set
#' and (k+1) columns: \code{Score1},...,\code{Scorek} (the score test statistics for each phenotype)
#' and \code{pvalue} (the multivariate score test p-value).
#' If a variant in the given variant-set has minor allele frequency = 0 or
#' greater than \code{rare_maf_cutoff}, the corresponding row will be \code{NA}. If a variant in
#' the given variant-set has degenerate covariance matrix across multiple phenotypes, the p-value will be set as 1.
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @export

Indiv_Score_Test_Region_multi <- function(genotype,obj_nullmodel,
                                          rare_maf_cutoff=0.01,rv_num_cutoff=2){

  if(class(genotype)[1] != "matrix" && !(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix")){
    stop("genotype is not a matrix!")
  }

  if(dim(genotype)[2] == 1){
    stop(paste0("Number of rare variant in the set is less than 2!"))
  }

  results <- data.frame(matrix(NA, nrow = dim(genotype)[2], ncol = obj_nullmodel$n.pheno + 1))
  colnames(results) <- c(paste0("Score",seq_len(obj_nullmodel$n.pheno)),"pvalue")

  if(!is.null(attr(class(genotype), "package")) && attr(class(genotype), "package") == "Matrix"){
    genotype <- as.matrix(genotype)
  }
  genotype <- matrix_flip(genotype)
  MAF <- genotype$MAF
  RV_label <- as.vector((MAF<rare_maf_cutoff)&(MAF>0))
  Geno_rare <- genotype$Geno[,RV_label]

  rm(genotype,MAF)
  gc()

  if(sum(RV_label) >= rv_num_cutoff){
    G <- as(Geno_rare,"dgCMatrix")
    G <- Diagonal(n = obj_nullmodel$n.pheno) %x% G
    rm(Geno_rare)
    gc()

    if(!obj_nullmodel$sparse_kins){
      P <- obj_nullmodel$P
      P_scalar <- sqrt(dim(P)[1])
      P <- P*P_scalar

      residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
      residuals.phenotype <- residuals.phenotype*sqrt(P_scalar)

      results[RV_label,] <- do.call(c,Indiv_Score_Test_SMMAT_multi(G,P,residuals.phenotype,obj_nullmodel$n.pheno))
    }else{
      Sigma_i <- obj_nullmodel$Sigma_i
      Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
      cov <- obj_nullmodel$cov

      residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)

      results[RV_label,] <- do.call(c,Indiv_Score_Test_SMMAT_multi_sparse(G,Sigma_i,Sigma_iX,cov,residuals.phenotype,obj_nullmodel$n.pheno))
    }

    return(results)
  }else{
    stop(paste0("Number of rare variant in the set is less than ",rv_num_cutoff,"!"))
  }

}

