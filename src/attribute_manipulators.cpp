#include <Rcpp.h>
#include <vector>
#include <string>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/posix_time/conversion.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
using namespace Rcpp;

// [[Rcpp::depends(BH)]]

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector createXts(std::vector<double> values_, std::vector<int> stamps_) {
  
  Rcpp::NumericVector values(values_.begin(), values_.end());
  Rcpp::NumericVector stamps(stamps_.begin(), stamps_.end());
  
  stamps.attr("tzone")    = "UTC";         // the index has attributes
  CharacterVector tclass = Rcpp::CharacterVector::create("POSIXct","POSIXt");
  stamps.attr("tclass")   = tclass;
  
  values.attr("dim") = Rcpp::IntegerVector::create(values.size(),1);
  values.attr("index") = stamps;
  CharacterVector klass  = Rcpp::CharacterVector::create("xts", "zoo");
  values.attr("class")       = klass;
  values.attr(".indexCLASS") = tclass;
  values.attr("tclass")      = tclass;
  values.attr(".indexTZ")    = "UTC";
  values.attr("tzone")       = "UTC";
  
  return values;
}

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector resizeNV_withAttributes(Rcpp::NumericVector vec1, Rcpp::NumericVector vec2){
  
  std::vector<double> outVec_std;
  for(Rcpp::NumericVector::iterator it = vec1.begin(); it != vec1.end(); ++it){
    outVec_std.push_back(*it);
  }
  for(Rcpp::NumericVector::iterator it = vec2.begin(); it != vec2.end(); ++it){
    outVec_std.push_back(*it);
  }
  // Rcpp::NumericVector outVec(outVec_std.begin(), outVec_std.end());
  
  std::vector<int> index_std;
  if(vec1.size() > 0){
    Rcpp::NumericVector index_vec1 = vec1.attr("index");  
    for(Rcpp::NumericVector::iterator it =index_vec1.begin(); it != index_vec1.end(); ++it){
      index_std.push_back(*it);
    }
  }
  
  if(vec2.size() > 0){
    Rcpp::NumericVector index_vec2 = vec2.attr("index");
    for(Rcpp::NumericVector::iterator it = index_vec2.begin(); it != index_vec2.end(); ++it){
      index_std.push_back(*it);
    }  
  }
  
  Rcpp::NumericVector outVec = createXts(outVec_std, index_std);
  
  return outVec;
}

std::vector<boost::gregorian::date> stripDatesFromXtsAttribute(const Rcpp::NumericVector &rdata, std::vector<boost::posix_time::ptime> &index){
  // Strip datetimes as vector of integers
  Rcpp::NumericVector rdataIndex_rcpp(rdata.attr("index"));
  
  // Convert to Boost POSIX
  // First declare std::vector of boost posix times
  std::vector<boost::posix_time::ptime> rdataIndex(rdataIndex_rcpp.length());
  
  for(int kk = 0; kk < rdataIndex.size(); kk++){
    rdataIndex[kk] = boost::posix_time::from_time_t(rdataIndex_rcpp[kk]);
  }
  index = rdataIndex;
  
  // get unique dates: first define a set, insert all into set, this will automatically get rid of redundant elements
  std::set<boost::gregorian::date> stripDates_set;
  for(std::vector<boost::posix_time::ptime>::iterator it = rdataIndex.begin(); it != rdataIndex.end(); ++it){
    stripDates_set.insert(it->date());
  }
  
  std::vector<boost::gregorian::date> stripDates;
  stripDates.assign(stripDates_set.begin(), stripDates_set.end());
  
  return stripDates;
}
