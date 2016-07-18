#include <RcppArmadillo.h>
#include "kernel_functions.h"
#include "../inst/include/attribute_manipulators.h"
using namespace Rcpp;

//[[Rcpp::export]]
arma::vec spotVolBaseJump_cpp(double spotPoint, const arma::vec& rdataSq, const arma::vec& rdataAbs, const arma::vec& rdataInd, double tRange, const arma::vec& timeStampYears, double avgVol, double referenceTime, bool sepLR, double timeDelta, double yearLength, std::string kernelType = "gaussian"){
  
  arma::vec timeKernel;
  if(kernelType == "gaussian"){
    timeKernel = kernel_gaussian(timeStampYears, (spotPoint-referenceTime)/(86400.0*yearLength), pow(timeDelta,1.5));
    // timeKernel = kernel_gaussian(timeStampYears, (spotPoint-referenceTime)/(86400.0*yearLength), 100*pow(timeDelta,2.0)); 
  } else if(kernelType == "epanechnikov"){
    timeKernel = kernel_epanechnikov(timeStampYears, (spotPoint-referenceTime)/(86400.0*yearLength), pow(timeDelta,1.5));
    // timeKernel = kernel_epanechnikov(timeStampYears, (spotPoint-referenceTime)/(86400.0*yearLength), 100*pow(timeDelta,2.0)); 
  } else if(kernelType == "indicator"){
    timeKernel = kernel_indicator(timeStampYears, (spotPoint-referenceTime)/(86400.0*yearLength), pow(timeDelta,1.5));
    // timeKernel = kernel_indicator(timeStampYears, (spotPoint-referenceTime)/(86400.0*yearLength), 100*pow(timeDelta,2.0)); 
  }
  
  arma::vec timeKernelMinus, timeKernelPlus;

  double resultMinus, resultPlus;
  
  if(sepLR){
    arma::vec timeKernelMinus = timeKernel % (rdataInd < spotPoint);
    arma::vec timeKernelPlus = timeKernel % (rdataInd >= spotPoint);
    
    resultMinus = arma::accu(timeKernelMinus % rdataSq % (rdataAbs < 2.5 * pow(avgVol,0.5) * pow(timeDelta,0.4999)));
    resultMinus /= arma::accu(timeKernelMinus.subvec(0,timeKernelMinus.n_elem-2) % arma::diff(timeStampYears) );
    
    resultPlus = arma::accu(timeKernelPlus % rdataSq % (rdataAbs < 2.5 * pow(avgVol,0.5) * pow(timeDelta,0.4999)));
    resultPlus /= arma::accu(timeKernelPlus.subvec(1,timeKernelPlus.n_elem-1) % arma::diff(timeStampYears) );
  } else {
    resultMinus = arma::accu(timeKernel % rdataSq % (rdataAbs < 2.5 * pow(avgVol,0.5) * pow(timeDelta,0.4999)));
    resultMinus /= arma::accu(timeKernel.subvec(0,timeKernel.n_elem-2) % arma::diff(timeStampYears) );
    
    resultPlus = resultMinus;
  }

  arma::vec res(2);
  res(0) = resultMinus;
  res(1) = resultPlus;
  return res;
}