#ifndef _KERNEL_FUNCTIONS_H_
#define _KERNEL_FUNCTIONS_H_
#include <RcppArmadillo.h>
using namespace Rcpp;

arma::vec kernel_gaussian(arma::vec cnvData, double x0, double delta);

arma::vec kernel_epanechnikov(arma::vec cnvData, double x0, double delta);

arma::vec kernel_indicator(arma::vec cnvData, double x0, double delta);

#endif