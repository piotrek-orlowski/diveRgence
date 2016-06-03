// Generalized CLT for Ito semimartingales -- Monte Carlo inference
//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
//'@useDynLib diveRgence
//'@export
// [[Rcpp::export]]
arma::vec mcCltInference(arma::vec& rdivDerivX, arma::vec& rdivDerivZ, double rdivCont, arma::vec& spotVolPlus, arma::vec& spotVolMinus, int nSampl){
  
  // Matrix dimensions
  int nObs = rdivDerivX.n_elem;
  
  // Create limiting sum vector
  arma::vec limitVarSum(nSampl);
  
  // Set RNG with R
  RNGScope rngScope;
  
  // Start filling -- first the square root of cont variance
  // arma::vec uContVol(1);
  // uContVol.fill(pow(rdivCont,0.5));
  rdivCont = pow(rdivCont,0.5);
    
  // Take square roots of spot volatilities
  spotVolPlus = arma::sqrt(spotVolPlus);
  spotVolMinus = arma::sqrt(spotVolMinus);  
    
  for(int nIter = 0; nIter < nSampl; nIter++){
    // Create matrices of random numbers;
    arma::vec uPlus(nObs, arma::fill::randn);
    arma::vec uMinus(nObs, arma::fill::randn);
    arma::vec kpp(nObs, arma::fill::randu);
    arma::vec oneLessKpp(nObs);
    arma::vec uCont(1L, arma::fill::randn);
    
    oneLessKpp.fill(1.0);
    oneLessKpp -= kpp;  
    
    // Create limiting variable matrix
    arma::vec limitVar(nObs);
    
    // Calculate the limiting variable
    limitVar = uPlus % arma::sqrt(kpp) % spotVolMinus % rdivDerivX;
    limitVar += uMinus % arma::sqrt(oneLessKpp) % spotVolPlus % rdivDerivX;
    limitVar -= uMinus % arma::sqrt(kpp) % spotVolMinus % rdivDerivZ;
    
    limitVarSum(nIter) = arma::sum(limitVar);
    // add variation coming from the cont part;
    limitVarSum(nIter) += rdivCont * uCont(0);    
  }

  // // Create matrices filled with inputs
//   arma::mat mat_rdivDerivX(nObs, nSampl);
//   arma::mat mat_rdivDerivZ(nObs, nSampl);
//   arma::mat mat_spotVolPlus(nObs, nSampl);
//   arma::mat mat_spotVolMinus(nObs, nSampl);

// // Write data vectors to large matrices
//   for(int nn=0; nn < nSampl; nn++){
//     mat_rdivDerivX.col(nn) = rdivDerivX.col(0);
//     mat_rdivDerivZ.col(nn) = rdivDerivZ.col(0);
//     mat_spotVolPlus.col(nn) = spotVolPlus.col(0);
//     mat_spotVolMinus.col(nn) = spotVolMinus.col(0);
//   }
// 
// // Calculate the limiting variable
//   mat_limitVar = uPlus % arma::sqrt(kpp) % mat_spotVolMinus % mat_rdivDerivX;
//   mat_limitVar += uMinus % arma::sqrt(oneLessKpp) % mat_spotVolPlus % mat_rdivDerivX;
//   mat_limitVar -= uMinus % arma::sqrt(kpp) % mat_spotVolMinus % mat_rdivDerivZ;
//   
// // Sum the matrix for each MC sample;
//   vec_limitVarSum = sum(mat_limitVar, 0).t();

// // add variation coming from the cont part;
//   vec_limitVarSum += uContVol % uCont;
  
// RETURN
  return(limitVarSum);
}