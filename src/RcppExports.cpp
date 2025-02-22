// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pfpop_interface
Rcpp::List pfpop_interface(const Rcpp::NumericVector degrees_vec, const double penalty, const Rcpp::NumericVector weight_vec);
RcppExport SEXP _pfpop_pfpop_interface(SEXP degrees_vecSEXP, SEXP penaltySEXP, SEXP weight_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type degrees_vec(degrees_vecSEXP);
    Rcpp::traits::input_parameter< const double >::type penalty(penaltySEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector >::type weight_vec(weight_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(pfpop_interface(degrees_vec, penalty, weight_vec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pfpop_pfpop_interface", (DL_FUNC) &_pfpop_pfpop_interface, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_pfpop(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
