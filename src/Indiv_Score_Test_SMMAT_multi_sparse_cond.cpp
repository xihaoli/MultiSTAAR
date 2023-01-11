// [[Rcpp::depends(RcppArmadillo)]]

// #define ARMA_64BIT_WORD 1
#include <STAAR.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
  
using namespace STAAR;
using namespace Rcpp;

// [[Rcpp::export]]
List Indiv_Score_Test_SMMAT_multi_sparse_cond(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::mat X_adj, arma::vec residuals, int n_pheno=1)
{
  int i,k;
  int p = G.n_cols;
  
  // number of markers
  int pp = p/n_pheno;
  
  // Uscore
  arma::rowvec Uscore = trans(residuals)*G;
  
  arma::vec pvalue;
  pvalue.zeros(pp);
  
  arma::mat Uscore_cov;
  Uscore_cov.zeros(n_pheno,n_pheno);
  
  arma::uvec id_single;
  id_single.zeros(n_pheno);
  
  double test_stat = 0;
  
  int q = Sigma_iX.n_cols;
  
  // arma::mat tSigma_iX_G;
  // tSigma_iX_G.zeros(q,p);
  
  arma::mat Cov;
  Cov.zeros(p,p);
  
  int PX_adj_col = X_adj.n_cols;
  int PX_adj_row = Sigma_i.n_rows;
  
  arma::mat PX_adj;
  PX_adj.zeros(PX_adj_row,PX_adj_col);
  PX_adj = Sigma_i*X_adj - Sigma_iX*cov*(trans(Sigma_iX)*X_adj);
  Cov = trans(Sigma_i*G)*G - trans(trans(Sigma_iX)*G)*cov*trans(Sigma_iX)*G - trans(trans(X_adj)*G)*inv(trans(X_adj)*X_adj)*(trans(PX_adj)*G) - trans(trans(PX_adj)*G)*inv(trans(X_adj)*X_adj)*(trans(X_adj)*G) + trans(trans(X_adj)*G)*inv(trans(X_adj)*X_adj)*trans(PX_adj)*X_adj*inv(trans(X_adj)*X_adj)*(trans(X_adj)*G);
  
  arma::mat quad;
  quad.zeros(1,1);
  
  for(i = 0; i < pp; i++)
  {
    for(k = 0; k < n_pheno; k++)
    {
      id_single(k) = k*pp+i;
    }
    
    Uscore_cov = Cov(id_single,id_single);
    
    if (arma::det(Uscore_cov) == 0)
    {
      pvalue(i) = 1;
    }
    else
    {
      quad = trans(Uscore(id_single))*inv(Uscore_cov)*Uscore(id_single);
      test_stat = quad(0,0);
      pvalue(i) = R::pchisq(test_stat,n_pheno,false,false);
    }
    
  }
  
  return List::create(Named("Uscore") = trans(Uscore), Named("pvalue") = pvalue);
}

