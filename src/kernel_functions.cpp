#include <RcppArmadillo.h>
using namespace Rcpp;

//[[Rcpp::export]]
arma::vec kernel_gaussian(arma::vec cnvData, double x0, double delta){
  
  cnvData -= x0;
  
  for(arma::vec::iterator res_it = cnvData.begin(); res_it != cnvData.end(); ++res_it){
    *res_it = R::dnorm4(*res_it, 0, pow(delta,0.5), FALSE);
  }
    
  return cnvData;
}

//[[Rcpp::export]]
arma::vec kernel_epanechnikov(arma::vec cnvData, double x0, double delta){
  
  cnvData -= x0;
  
  for(arma::vec::iterator res_it = cnvData.begin(); res_it != cnvData.end(); ++res_it){
    *res_it = std::pow(delta,-0.5) * 0.75 * (1 - std::pow((*res_it),2.0)/delta) * (std::abs(*res_it) <= std::pow(delta,0.5));
  }
  
  return cnvData;
}


//[[Rcpp::export]]
arma::vec kernel_indicator(arma::vec cnvData, double x0, double delta){
  
  cnvData -= x0;
  
  for(arma::vec::iterator res_it = cnvData.begin(); res_it != cnvData.end(); ++res_it){
    *res_it = pow(delta,-0.5) * 0.5 * (std::abs(*res_it) <= std::pow(delta,0.5));
  }
  
  return cnvData;
}
