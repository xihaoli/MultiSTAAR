// [[Rcpp::depends(RcppArmadillo)]]

#include <STAAR.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace STAAR;
using namespace Rcpp;

// [[Rcpp::export]]
List Indiv_Score_Test_SMMAT_multi(arma::sp_mat G, arma::mat P, arma::vec residuals, int n_pheno=1)
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

	arma::mat Cov;
	Cov.zeros(p,p);

	Cov = trans(P*G)*G;

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

