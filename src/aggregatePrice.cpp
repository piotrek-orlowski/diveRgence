// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
#include <vector>
#include <string>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/conversion.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include "../inst/include/attribute_manipulators.h"
using namespace Rcpp;

// This function aggregates an xts R time series to a given frequency by the 
// `lasttick' method. It handles multiple days and aggregation is done on a
// within-day basis.

////' @export
// [[Rcpp::export]]
Rcpp::NumericVector aggregatePrice_Xts(Rcpp::NumericVector& rdata, std::string period_, int numPeriods_, Rcpp::NumericVector dayStart_, Rcpp::NumericVector dayEnd_, Rcpp::IntegerVector aggr_vec, bool pad = true, double pad_arg = 0.0){
  
  // Strip datetimes as vector of integers
  Rcpp::IntegerVector rdataIndex_rcpp(rdata.attr("index"));
  std::vector<int> rdataIndex_time_t(rdataIndex_rcpp.begin(), rdataIndex_rcpp.end());

  // check prescribed aggr dates
  int agVecLen = aggr_vec.size();

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

  // Put rdata in std vector
  std::vector<double> rdata_std(rdata.begin(), rdata.end());

  // Make very big grid
  int gridStep = 0;
  if(!period_.compare(std::string("seconds"))){
    gridStep = numPeriods_;
  } else if(!period_.compare(std::string("minutes"))){
    gridStep = numPeriods_ * 60L;
  } else if(!period_.compare(std::string("hours"))){
    gridStep = numPeriods_ * 3600L;
  }

  boost::posix_time::ptime dayStart_posix(*stripDates.rbegin(),boost::posix_time::time_duration(dayStart_[0],dayStart_[1],dayStart_[2]));
  boost::posix_time::ptime dayEnd_posix(*stripDates.rbegin(),boost::posix_time::time_duration(dayEnd_[0],dayEnd_[1],dayEnd_[2]));  
  time_t trueDayStart_seconds = boost::posix_time::to_time_t(dayStart_posix);
  time_t trueDayEnd_seconds = boost::posix_time::to_time_t(dayEnd_posix);

  int numSteps = (trueDayEnd_seconds - trueDayStart_seconds) / gridStep;
  
  int gridSize = (numSteps + 1) * stripDates.size();
  
  int where_on_grid = 0L;
  
  std::vector<double> rdataOut(gridSize);
  std::vector<int> timeGrid(gridSize);
  
  // Loop over stripDates and create a time grid
  for(std::vector<boost::gregorian::date>::iterator it = stripDates.begin(); it != stripDates.end(); ++it){
    
    boost::posix_time::ptime locDate_posix(*it,boost::posix_time::seconds(0.0));
    
    // convert general day start and end times to posixs
    boost::posix_time::ptime locDayStart_posix(*it,boost::posix_time::time_duration(dayStart_[0],dayStart_[1],dayStart_[2]));
    boost::posix_time::ptime locDayEnd_posix(*it,boost::posix_time::time_duration(dayEnd_[0],dayEnd_[1],dayEnd_[2]));
    
    trueDayStart_seconds = boost::posix_time::to_time_t(locDayStart_posix);
    trueDayEnd_seconds = boost::posix_time::to_time_t(locDayEnd_posix);

    std::vector<int> timeStamps(2+numSteps);
      
    timeStamps[0] = trueDayStart_seconds;
    timeStamps[numSteps+1] = trueDayEnd_seconds;
    for(int kk = 1; kk < numSteps+1; kk++){
      timeStamps[kk] = trueDayStart_seconds + kk * gridStep;
    }

    // Remove last stamp if equal to penultimate stamp
    if(timeStamps.back() == timeStamps[timeStamps.size()-2]){
      timeStamps.pop_back();
    }
    for(int kk = where_on_grid; kk < where_on_grid + timeStamps.size(); kk++){
      timeGrid[kk] = timeStamps[kk-where_on_grid];
    }
    where_on_grid = where_on_grid + timeStamps.size();
  }
  
  if(agVecLen > 0L){
    timeGrid.resize(aggr_vec.size());
    rdataOut.resize(aggr_vec.size());
    timeGrid = Rcpp::as<std::vector<int>>(aggr_vec);
  }
  
  // Loop backwards over the time grid and pop unnecessary values from rdata_std, write the necessary ones into rdataOut
  where_on_grid = 0L;
  
  // if pad_na = true, put pad_args in times after the last observation
  if(pad){
    for(std::vector<int>::reverse_iterator it = timeGrid.rbegin(); it != timeGrid.rend(); ++it){
      
      int locTimeStamp = *it;
      while(rdataIndex_time_t.back() > locTimeStamp){
        rdataIndex_time_t.pop_back();
        rdata_std.pop_back();
      }
      
      // if(locTimeStamp > rdataIndex_time_t.back() + gridStep){
      if(rdataIndex_time_t.back() <= *std::next(it,1L)){
        rdataOut[where_on_grid] = pad_arg;
      } else {
        rdataOut[where_on_grid] = rdata_std.back();  
      }  
      
      ++where_on_grid;
    }
  } else {
    for(std::vector<int>::reverse_iterator it = timeGrid.rbegin(); it != timeGrid.rend(); ++it){
      
      int locTimeStamp = *it;
      while(rdataIndex_time_t.back() > locTimeStamp){
        rdataIndex_time_t.pop_back();
        rdata_std.pop_back();
      }
      
      rdataOut[where_on_grid] = rdata_std.back();  
      
      ++where_on_grid;
    }
  }
  
  std::vector<double> rdataOut_ordered(rdataOut.rbegin(), rdataOut.rend());
  
  NumericVector resultVec = createXts(rdataOut_ordered, timeGrid);
  
  return resultVec;
}
