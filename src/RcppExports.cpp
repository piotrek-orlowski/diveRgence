// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// aggregatePrice_Xts
Rcpp::NumericVector aggregatePrice_Xts(Rcpp::NumericVector& rdata, std::string period_, int numPeriods_, Rcpp::NumericVector dayStart_, Rcpp::NumericVector dayEnd_, Rcpp::IntegerVector aggr_vec, bool pad, double pad_arg);
RcppExport SEXP diveRgence_aggregatePrice_Xts(SEXP rdataSEXP, SEXP period_SEXP, SEXP numPeriods_SEXP, SEXP dayStart_SEXP, SEXP dayEnd_SEXP, SEXP aggr_vecSEXP, SEXP padSEXP, SEXP pad_argSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector& >::type rdata(rdataSEXP);
    Rcpp::traits::input_parameter< std::string >::type period_(period_SEXP);
    Rcpp::traits::input_parameter< int >::type numPeriods_(numPeriods_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type dayStart_(dayStart_SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type dayEnd_(dayEnd_SEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type aggr_vec(aggr_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type pad(padSEXP);
    Rcpp::traits::input_parameter< double >::type pad_arg(pad_argSEXP);
    __result = Rcpp::wrap(aggregatePrice_Xts(rdata, period_, numPeriods_, dayStart_, dayEnd_, aggr_vec, pad, pad_arg));
    return __result;
END_RCPP
}
// createXts
Rcpp::NumericVector createXts(std::vector<double> values_, std::vector<int> stamps_);
RcppExport SEXP diveRgence_createXts(SEXP values_SEXP, SEXP stamps_SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< std::vector<double> >::type values_(values_SEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type stamps_(stamps_SEXP);
    __result = Rcpp::wrap(createXts(values_, stamps_));
    return __result;
END_RCPP
}
// resizeNV_withAttributes
Rcpp::NumericVector resizeNV_withAttributes(Rcpp::NumericVector vec1, Rcpp::NumericVector vec2);
RcppExport SEXP diveRgence_resizeNV_withAttributes(SEXP vec1SEXP, SEXP vec2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type vec2(vec2SEXP);
    __result = Rcpp::wrap(resizeNV_withAttributes(vec1, vec2));
    return __result;
END_RCPP
}
// mcCltInference
arma::vec mcCltInference(arma::vec& rdivDerivX, arma::vec& rdivDerivZ, double rdivCont, arma::vec& spotVolPlus, arma::vec& spotVolMinus, int nSampl);
RcppExport SEXP diveRgence_mcCltInference(SEXP rdivDerivXSEXP, SEXP rdivDerivZSEXP, SEXP rdivContSEXP, SEXP spotVolPlusSEXP, SEXP spotVolMinusSEXP, SEXP nSamplSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec& >::type rdivDerivX(rdivDerivXSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type rdivDerivZ(rdivDerivZSEXP);
    Rcpp::traits::input_parameter< double >::type rdivCont(rdivContSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type spotVolPlus(spotVolPlusSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type spotVolMinus(spotVolMinusSEXP);
    Rcpp::traits::input_parameter< int >::type nSampl(nSamplSEXP);
    __result = Rcpp::wrap(mcCltInference(rdivDerivX, rdivDerivZ, rdivCont, spotVolPlus, spotVolMinus, nSampl));
    return __result;
END_RCPP
}
// kernel_gaussian
arma::vec kernel_gaussian(arma::vec cnvData, double x0, double delta);
RcppExport SEXP diveRgence_kernel_gaussian(SEXP cnvDataSEXP, SEXP x0SEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type cnvData(cnvDataSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    __result = Rcpp::wrap(kernel_gaussian(cnvData, x0, delta));
    return __result;
END_RCPP
}
// kernel_epanechnikov
arma::vec kernel_epanechnikov(arma::vec cnvData, double x0, double delta);
RcppExport SEXP diveRgence_kernel_epanechnikov(SEXP cnvDataSEXP, SEXP x0SEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type cnvData(cnvDataSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    __result = Rcpp::wrap(kernel_epanechnikov(cnvData, x0, delta));
    return __result;
END_RCPP
}
// kernel_indicator
arma::vec kernel_indicator(arma::vec cnvData, double x0, double delta);
RcppExport SEXP diveRgence_kernel_indicator(SEXP cnvDataSEXP, SEXP x0SEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::vec >::type cnvData(cnvDataSEXP);
    Rcpp::traits::input_parameter< double >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    __result = Rcpp::wrap(kernel_indicator(cnvData, x0, delta));
    return __result;
END_RCPP
}
// rMPVcpp
Rcpp::NumericVector rMPVcpp(const Rcpp::NumericVector& rdata, double mNum, double pPow, double yearDays);
RcppExport SEXP diveRgence_rMPVcpp(SEXP rdataSEXP, SEXP mNumSEXP, SEXP pPowSEXP, SEXP yearDaysSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type rdata(rdataSEXP);
    Rcpp::traits::input_parameter< double >::type mNum(mNumSEXP);
    Rcpp::traits::input_parameter< double >::type pPow(pPowSEXP);
    Rcpp::traits::input_parameter< double >::type yearDays(yearDaysSEXP);
    __result = Rcpp::wrap(rMPVcpp(rdata, mNum, pPow, yearDays));
    return __result;
END_RCPP
}
// spotVolBaseJump_cpp
arma::vec spotVolBaseJump_cpp(double spotPoint, const arma::vec& rdataSq, const arma::vec& rdataAbs, const arma::vec& rdataInd, double tRange, const arma::vec& timeStampYears, double avgVol, double referenceTime, bool sepLR, double timeDelta, double yearLength, std::string kernelType);
RcppExport SEXP diveRgence_spotVolBaseJump_cpp(SEXP spotPointSEXP, SEXP rdataSqSEXP, SEXP rdataAbsSEXP, SEXP rdataIndSEXP, SEXP tRangeSEXP, SEXP timeStampYearsSEXP, SEXP avgVolSEXP, SEXP referenceTimeSEXP, SEXP sepLRSEXP, SEXP timeDeltaSEXP, SEXP yearLengthSEXP, SEXP kernelTypeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< double >::type spotPoint(spotPointSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rdataSq(rdataSqSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rdataAbs(rdataAbsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rdataInd(rdataIndSEXP);
    Rcpp::traits::input_parameter< double >::type tRange(tRangeSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type timeStampYears(timeStampYearsSEXP);
    Rcpp::traits::input_parameter< double >::type avgVol(avgVolSEXP);
    Rcpp::traits::input_parameter< double >::type referenceTime(referenceTimeSEXP);
    Rcpp::traits::input_parameter< bool >::type sepLR(sepLRSEXP);
    Rcpp::traits::input_parameter< double >::type timeDelta(timeDeltaSEXP);
    Rcpp::traits::input_parameter< double >::type yearLength(yearLengthSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernelType(kernelTypeSEXP);
    __result = Rcpp::wrap(spotVolBaseJump_cpp(spotPoint, rdataSq, rdataAbs, rdataInd, tRange, timeStampYears, avgVol, referenceTime, sepLR, timeDelta, yearLength, kernelType));
    return __result;
END_RCPP
}
