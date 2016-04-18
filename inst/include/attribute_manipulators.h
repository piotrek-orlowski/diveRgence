#ifndef ATTRIBUTE_MANIPULATORS_H
#define ATTRIBUTE_MANIPULATORS_H

#include <Rcpp.h>
#include <vector>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

Rcpp::NumericVector createXts(std::vector<double> values_, std::vector<double> stamps_);

Rcpp::NumericVector resizeNV_withAttributes(Rcpp::NumericVector vec1, Rcpp::NumericVector vec2);

std::vector<boost::gregorian::date> stripDatesFromXtsAttribute(const Rcpp::NumericVector& rdata, std::vector<boost::posix_time::ptime>& index);

#endif
