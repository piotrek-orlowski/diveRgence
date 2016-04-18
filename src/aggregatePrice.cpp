// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <boost\date_time\posix_time\posix_time.hpp>
#include <boost\date_time\posix_time\conversion.hpp>
#include <boost\date_time\gregorian\gregorian.hpp>
#include "..\inst\include\attribute_manipulators.h"
using namespace Rcpp;

// This function aggregates an xts R time series to a given frequency by the 
// `lasttick' method. It handles multiple days and aggregation is done on a
// within-day basis.

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector aggregatePrice_Xts(Rcpp::NumericVector rdata, std::string period_, int numPeriods_, Rcpp::NumericVector dayStart_, Rcpp::NumericVector dayEnd_){
  
  // Check for rdata being an xts object with POSIX time
  
  // declare result
  Rcpp::NumericVector resultVec(0);
  
  // Strip datetimes as vector of integers
  Rcpp::NumericVector rdataIndex_rcpp(rdata.attr("index"));
  
  // Convert to Boost POSIX
  // First declare std::vector of boost posix times
  std::vector<boost::posix_time::ptime> rdataIndex(rdataIndex_rcpp.length());
  
  for(int kk = 0; kk < rdataIndex.size(); kk++){
    rdataIndex[kk] = boost::posix_time::from_time_t(rdataIndex_rcpp[kk]);
  }
  
  // get unique dates: first define a set, insert all into set, this will automatically get rid of redundant elements
  std::set<boost::gregorian::date> stripDates_set;
  for(std::vector<boost::posix_time::ptime>::iterator it = rdataIndex.begin(); it != rdataIndex.end(); ++it){
    stripDates_set.insert(it->date());
  }
  
  std::vector<boost::gregorian::date> stripDates;
  stripDates.assign(stripDates_set.begin(), stripDates_set.end());
 
  // loop over stripDates; for each date subset rdata for the requisite day,
  // find the starting and ending time values in the day, create limiting times
  // from dayStart_ and dayEnd_ and check if they're not too far from the day's
  // official values -- if yes, move them to the factual values.
  for(std::vector<boost::gregorian::date>::iterator it = stripDates.begin(); it != stripDates.end(); ++it){
    // get local date as posix
    boost::posix_time::ptime locDate_posix(*it,boost::posix_time::seconds(0.0));
    time_t locDate_seconds = boost::posix_time::to_time_t(locDate_posix);
    
    // subset rdata
    Rcpp::LogicalVector rdataSub_selector = rdataIndex_rcpp >= locDate_seconds & rdataIndex_rcpp < (locDate_seconds + 86400);
    Rcpp::NumericVector rdataSubIndex_rcpp = rdataIndex_rcpp[rdataSub_selector];
    
    // convert general day start and end times to posixs
    boost::posix_time::ptime locDayStart_posix(*it,boost::posix_time::time_duration(dayStart_[0],dayStart_[1],dayStart_[2]));
    boost::posix_time::ptime locDayEnd_posix(*it,boost::posix_time::time_duration(dayEnd_[0],dayEnd_[1],dayEnd_[2]));
    
    time_t locDayStart_seconds = boost::posix_time::to_time_t(locDayStart_posix);
    time_t locDayEnd_seconds = boost::posix_time::to_time_t(locDayEnd_posix);
    
    // get first and last time stamp of day of data
    time_t trueDayStart_seconds = *rdataSubIndex_rcpp.begin();
    time_t trueDayEnd_seconds = rdataSubIndex_rcpp[rdataSubIndex_rcpp.size()-1];
    
    // compare with required timestamps and choose what's reasonable
    trueDayStart_seconds = std::max(trueDayStart_seconds, locDayStart_seconds);
    trueDayEnd_seconds = std::min(trueDayEnd_seconds, locDayEnd_seconds);
    
    // further subset rdata
    rdataSub_selector = rdataIndex_rcpp >= trueDayStart_seconds & rdataIndex_rcpp <= trueDayEnd_seconds;
    Rcpp::NumericVector rdataSub_rcpp = rdata[rdataSub_selector];
    rdataSubIndex_rcpp = rdataIndex_rcpp[rdataSub_selector];
    std::vector<double> rdataSubIndex(rdataSubIndex_rcpp.begin(), rdataSubIndex_rcpp.end());
    
    std::vector<double> rdataSub(rdataSub_rcpp.begin(), rdataSub_rcpp.end());
    
    // make aggregation grid
    int gridStep = 0;
    if(!period_.compare(std::string("seconds"))){
      gridStep = numPeriods_;
    } else if(!period_.compare(std::string("minutes"))){
      gridStep = numPeriods_ * 60L;
    } else if(!period_.compare(std::string("hours"))){
      gridStep = numPeriods_ * 3600L;
    }
    
    int numSteps = (trueDayEnd_seconds - trueDayStart_seconds) / gridStep;
    
    std::vector<double> timeStamps(2+numSteps);
    timeStamps[0] = trueDayStart_seconds;
    timeStamps[numSteps+1] = trueDayEnd_seconds;
    for(int kk = 1; kk < numSteps+1; kk++){
      timeStamps[kk] = trueDayStart_seconds + kk * gridStep;
    }
    // Remove last stamp if equal to penultimate stamp
    if(timeStamps.back() == timeStamps[timeStamps.size()-2]){
      timeStamps.pop_back();
    }
    
    std::vector<double> rdataOutLoc_reversed;
    std::vector<double> rdataOutStamps_reversed;
    
    // loop over the grid from the end and always find the latest value in date vector, pick corresponding price value
    for(std::vector<double>::reverse_iterator bck = timeStamps.rbegin(); bck != timeStamps.rend(); ++bck){
      double locTimeStamp = *bck;
      
      while(rdataSubIndex.back() > locTimeStamp){
        rdataSubIndex.pop_back();
        rdataSub.pop_back();
      }
      rdataOutLoc_reversed.push_back(rdataSub.back());
      rdataOutStamps_reversed.push_back(locTimeStamp);
    }
    
    // The vectors that you got are reversed. Put them in the right order
    std::vector<double> rdataOutLoc(rdataOutLoc_reversed.rbegin(), rdataOutLoc_reversed.rend());
    std::vector<double> rdataOutStamps(rdataOutStamps_reversed.rbegin(), rdataOutStamps_reversed.rend());
    
    // Create xts-type object
    Rcpp::NumericVector resultLoc = createXts(rdataOutLoc, rdataOutStamps);
    
    // Resize result, bind
    resultVec = resizeNV_withAttributes(resultVec, resultLoc);
  }
  
  return resultVec;
}
