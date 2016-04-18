#include <Rcpp.h>
#include <vector>
#include <string>
#include <iterator>
#include <cmath>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/conversion.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include "../inst/include/attribute_manipulators.h"
using namespace Rcpp;

// [[Rcpp::depends(BH)]]

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector rMPVcpp(const Rcpp::NumericVector& rdata, double mNum, double pPow, double yearDays){
  // effective power for individual increments
  double pEff = pPow/mNum;
  
  // normalisation for absolute value in normal distribution
  // pi^(-0.5) * 2^(pEff/2) * gamma((pEff+1)/2)
  double muPow = std::pow(M_PI,-0.5) * std::pow(2.0, pEff/2.0) * std::tgamma((pEff + 1.0)/2.0);
  
  // declare variable for index
  std::vector<boost::posix_time::ptime> rdataIndex;
  
  // strip dates, store posix index in rdataIndex at the same time
  std::vector<boost::gregorian::date> stripDate = stripDatesFromXtsAttribute(rdata, rdataIndex);
  
  // std vector copy of input
  std::vector<double> rdata_std(rdata.begin(), rdata.end());
  
  // declare output
  Rcpp::NumericVector outRes;
  
  int mNumInt = mNum;
  
  for(std::vector<boost::gregorian::date>::reverse_iterator it = stripDate.rbegin(); it != stripDate.rend(); ++it){
    // get local date as posix
    boost::posix_time::ptime locDate_posix(*it,boost::posix_time::seconds(0.0));
    
    // declare local value vector
    std::vector<double> rdataSub_reversed(0);
    std::vector<boost::posix_time::ptime> locTimeStamp_reversed;
    std::vector<boost::posix_time::ptime>::reverse_iterator bck = rdataIndex.rbegin();
    
    // collect local date values, delete them from the input as you go; 
    while(*bck >= locDate_posix & !rdataIndex.empty()){
      rdataSub_reversed.push_back(rdata_std.back());
      rdata_std.pop_back();
      locTimeStamp_reversed.push_back(rdataIndex.back());
      rdataIndex.pop_back();
      ++bck;
    }
    
    // calculate mean time increment
    boost::posix_time::time_duration meanTimeIncrement_td = locTimeStamp_reversed.front() - locTimeStamp_reversed.back();
    double meanTimeIncrement = meanTimeIncrement_td.total_seconds() + meanTimeIncrement_td.fractional_seconds();
    meanTimeIncrement /= locTimeStamp_reversed.size();
    meanTimeIncrement /=(86400.0 * yearDays);
    
    // take powers of elements of rdataSub
    for(std::vector<double>::iterator rdataIt = rdataSub_reversed.begin(); rdataIt != rdataSub_reversed.end(); ++rdataIt){
      *rdataIt = std::pow(std::abs(*rdataIt),pEff);
    }
    std::vector<double> rdataAccu(1);
    rdataAccu[0] = 0.0;
    
    for(std::vector<double>::reverse_iterator rdataIt = rdataSub_reversed.rbegin(); rdataIt != std::next(rdataSub_reversed.rend(),-mNumInt); ++rdataIt){
      rdataAccu[0] += std::accumulate(rdataIt, std::next(rdataIt,mNumInt), 1.0, std::multiplies<double>());
    }
    // normalise
    rdataAccu[0] *= pow(meanTimeIncrement,1.0 - pPow/2.0);
    rdataAccu[0] *= pow(muPow,-mNum);
    double nObs = rdataSub_reversed.size();
    rdataAccu[0] *= nObs/(nObs - mNum + 1);
    
    std::vector<double> locStampVec(1);
    locStampVec[0] = boost::posix_time::to_time_t(locTimeStamp_reversed.front());
    Rcpp::NumericVector locOut = createXts(rdataAccu, locStampVec);
    outRes = resizeNV_withAttributes(outRes, locOut);
  }

  // Reverse the output vector
  std::reverse(outRes.begin(), outRes.end());
  // also reverse its index
  Rcpp::NumericVector finalIndex = outRes.attr("index");
  std::reverse(finalIndex.begin(), finalIndex.end());
  outRes.attr("index") = finalIndex;
  
  return outRes;
}
