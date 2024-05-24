// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// bpp
double bpp(int m, double alphat, double betat, double alphac, double betac, double eta, double delta);
RcppExport SEXP _bayespredp_bpp(SEXP mSEXP, SEXP alphatSEXP, SEXP betatSEXP, SEXP alphacSEXP, SEXP betacSEXP, SEXP etaSEXP, SEXP deltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type alphat(alphatSEXP);
    Rcpp::traits::input_parameter< double >::type betat(betatSEXP);
    Rcpp::traits::input_parameter< double >::type alphac(alphacSEXP);
    Rcpp::traits::input_parameter< double >::type betac(betacSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    rcpp_result_gen = Rcpp::wrap(bpp(m, alphat, betat, alphac, betac, eta, delta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bayespredp_bpp", (DL_FUNC) &_bayespredp_bpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_bayespredp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
