// [[Rcpp::depends(RcppArmadillo)]]

#include <STAAR.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace STAAR;
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec MultiSTAAR_O_SMMAT(arma::sp_mat G, arma::mat P, arma::vec residuals, arma::mat weights_B, arma::mat weights_S, arma::mat weights_A, arma::vec mac, int mac_thres=10, int n_pheno=1)
{
  int i,k,kk,kk2;
  int j;
  double sum0 = 0.0;
  double sumw = 0.0;
  double sumx = 0.0;
  
  // ACAT
  arma::rowvec sum0acat;
  sum0acat.zeros(n_pheno);
  
  // Burden	
  arma::rowvec sum0burden;	
  sum0burden.zeros(n_pheno);
  
  bool lower = false;
  bool logp = false;
  
  
  // int vn = G.n_rows;
  int un = G.n_cols;
  int wn = weights_B.n_cols;
  
  arma::vec res;
  res.zeros(3*wn);
  
  arma::mat Cov;
  Cov.zeros(un,un);
  
  arma::mat Covw;
  Covw.zeros(un,un);
  
  // cov: ACAT
  arma::mat Covwacat;
  Covwacat.zeros(n_pheno,n_pheno);
  // cov: Burden	
  arma::mat Covwburden;	
  Covwburden.zeros(n_pheno,n_pheno);
  
  arma::rowvec x = trans(residuals)*G;
  
  arma::vec eigenvals;
  eigenvals.zeros(un);
  
  
  arma::mat Wleft;
  Wleft.zeros(un,un);
  
  arma::mat Wright;
  Wright.zeros(un,un);
  
  double c1 = 0.0;
  double c2 = 0.0;
  double c4 = 0.0;
  double l = 0.0;
  
  int ii;
  
  Cov = trans(P*G)*G;
  
  
  // ACAT
  arma::uvec id_veryrare = arma::find(mac <= mac_thres);
  arma::uvec id_common = arma::find(mac > mac_thres);
  
  // Burden	
  arma::uvec id_all = arma::find(mac > 0);
  
  int n0 = id_veryrare.size()/n_pheno;
  int n1 = id_common.size()/n_pheno;
  int uun = un/n_pheno;
  
  arma::vec pseq;
  pseq.zeros(uun);
  
  arma::vec wseq;
  wseq.zeros(uun);
  
  arma::uvec id_common_single;
  id_common_single.zeros(n_pheno);
  
  arma::mat quad;
  quad.zeros(1,1);
  
  // double SSR = 0.0;
  // double SST = 0.0;
  
  for(k = 0; k < n1; k++)
  {
    for(kk = 0; kk < n_pheno; kk++)
    {
      id_common_single(kk) = id_common(k) + kk*uun;
    }
    quad = trans(x(id_common_single))*inv(Cov(id_common_single,id_common_single))*x(id_common_single);
    pseq(k) = quad(0,0);
    pseq(k) = R::pchisq(pseq(k),n_pheno,lower,logp);
  }
  
  
  int n = x.size();
  int n0_plus_n1 = n/n_pheno;
  
  for(i = 0; i < wn; i++)
  {
    // SKAT
    Wright.each_row() = trans(weights_S.col(i));
    Wleft.each_col() = weights_S.col(i);
    
    Covw = Wleft%Cov%Wright;
    
    sum0 = 0.0;
    for (k = 0; k < n; k++)
    {
      sum0 = sum0 + pow(x(k), 2) *pow(weights_S(k,i), 2);
    }
    
    eigenvals = arma::eig_sym(Covw);
    for(j = 0; j < un; j++)
    {
      if(eigenvals(j) < 1e-8)
      {
        eigenvals(j) = 0.0;
      }
    }
    
    
    res(i) = Saddle(sum0,eigenvals);
    
    if(res(i)== 2)
    {
      c1 = 0.0;
      c2 = 0.0;
      c4 = 0.0;
      
      for(ii = 0; ii < un; ii++)
      {
        c1 = c1 + Covw(ii, ii);
      }
      
      Covw = Covw*Covw;
      
      for(ii = 0; ii < un; ii++)
      {
        c2 = c2 + Covw(ii, ii);
      }
      
      Covw = Covw*Covw;
      
      for(ii = 0; ii < un; ii++)
      {
        c4 = c4 + Covw(ii, ii);
      }
      
      sum0 = (sum0 - c1)/sqrt(2*c2);
      l = pow(c2,2)/c4;
      res(i) = R::pchisq(sum0*sqrt(2*l)+l,l,lower,logp);
    }
    
    // Burden
    Wright.each_row() = trans(weights_B.col(i));
    Wleft.each_col() = weights_B.col(i);
    
    Covw = Wleft%Cov%Wright;
    
    sum0burden.zeros(n_pheno);	
    for(kk = 0; kk < n_pheno; kk++)	
    {	
      for (k = 0; k < n0_plus_n1; k++)	
      {	
        sum0burden(kk) = sum0burden(kk) + x(k+kk*uun)*weights_B(k,i);	
      }	
    }	
    
    for(kk = 0; kk < n_pheno; kk++)	
    {	
      for(kk2 = 0; kk2 < n_pheno; kk2++)	
      {	
        Covwburden(kk,kk2) = arma::accu(Covw.submat(id_all.subvec(kk*n0_plus_n1,(kk+1)*n0_plus_n1-1),	
                                                    id_all.subvec(kk2*n0_plus_n1,(kk2+1)*n0_plus_n1-1)));	
      }	
    }	
    
    quad = sum0burden*inv(Covwburden)*trans(sum0burden);	
    res(wn + i) = quad(0,0);	
    res(wn + i) = R::pchisq(res(wn + i),n_pheno,lower,logp);
    
    // ACAT
    for(k = 0; k < n1; k++)
    {
      wseq(k) = weights_A(id_common(k),i);
    }
    
    if(n0 == 0)
    {
      res(2*wn + i) = CCT_pval(pseq(arma::span(0,n1-1)),wseq(arma::span(0,n1-1)));
    }else
    {
      sum0acat.zeros(n_pheno);
      sumw = 0.0;
      for(kk = 0; kk < n_pheno; kk++)
      {
        for (k = 0; k < n0; k++)
        {
          sum0acat(kk) = sum0acat(kk) + x(id_veryrare(k)+kk*uun)*weights_B(id_veryrare(k),i);
        }
      }
      for(kk = 0; kk < n_pheno; kk++)
      {
        for(kk2 = 0; kk2 < n_pheno; kk2++)
        {
          Covwacat(kk,kk2) = arma::accu(Covw.submat(id_veryrare.subvec(kk*n0,(kk+1)*n0-1),
                                                    id_veryrare.subvec(kk2*n0,(kk2+1)*n0-1)));
        }
      }
      for (k = 0; k < n0; k++)
      {
        sumw = sumw + weights_A(id_veryrare(k),i);
      }
      
      quad = sum0acat*inv(Covwacat)*trans(sum0acat);
      pseq(n1) = quad(0,0);
      pseq(n1) = R::pchisq(pseq(n1),n_pheno,lower,logp);
      wseq(n1) = sumw/n0;
      res(2*wn + i) = CCT_pval(pseq(arma::span(0,n1)),wseq(arma::span(0,n1)));
    }
    
    
  }
  
  
  return res;
}

