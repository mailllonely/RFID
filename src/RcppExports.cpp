// write by Jun Li 06-12-2019
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// log_likelihood_Y
double log_likelihood_Y(NumericMatrix Y, NumericMatrix PhiTheta, NumericVector alpha, NumericVector p, double lambda);
RcppExport SEXP _RFID_log_likelihood_Y(SEXP YSEXP, SEXP PhiThetaSEXP, SEXP alphaSEXP, SEXP pSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PhiTheta(PhiThetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(log_likelihood_Y(Y, PhiTheta, alpha, p, lambda));
    return rcpp_result_gen;
END_RCPP
}

double log_likelihood_Y_purity(NumericMatrix Y, NumericMatrix PhiTheta, NumericVector alpha, NumericVector p, double lambda);
RcppExport SEXP _RFID_log_likelihood_Y_purity(SEXP YSEXP, SEXP PhiThetaSEXP, SEXP alphaSEXP, SEXP pSEXP, SEXP lambdaSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PhiTheta(PhiThetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(log_likelihood_Y_purity(Y, PhiTheta, alpha, p, lambda));
    return rcpp_result_gen;
    END_RCPP
}


// likelihood_mat
NumericMatrix likelihood_mat(NumericMatrix Y, NumericMatrix PhiTheta, NumericVector alpha, NumericVector p, double lambda);
RcppExport SEXP _RFID_likelihood_mat(SEXP YSEXP, SEXP PhiThetaSEXP, SEXP alphaSEXP, SEXP pSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PhiTheta(PhiThetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(likelihood_mat(Y, PhiTheta, alpha, p, lambda));
    return rcpp_result_gen;
END_RCPP
}
// log_likelihood_alpha
double log_likelihood_alpha(NumericVector alpha, double g0, double h0);
RcppExport SEXP _RFID_log_likelihood_alpha(SEXP alphaSEXP, SEXP g0SEXP, SEXP h0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< double >::type h0(h0SEXP);
    rcpp_result_gen = Rcpp::wrap(log_likelihood_alpha(alpha, g0, h0));
    return rcpp_result_gen;
END_RCPP
}
// log_likelihood_Phi
double log_likelihood_Phi(NumericMatrix Phi, double eta);
RcppExport SEXP _RFID_log_likelihood_Phi(SEXP PhiSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(log_likelihood_Phi(Phi, eta));
    return rcpp_result_gen;
END_RCPP
}
// log_likelihood_Theta
double log_likelihood_Theta(NumericMatrix Theta, NumericVector r, NumericVector c);
RcppExport SEXP _RFID_log_likelihood_Theta(SEXP ThetaSEXP, SEXP rSEXP, SEXP cSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    rcpp_result_gen = Rcpp::wrap(log_likelihood_Theta(Theta, r, c));
    return rcpp_result_gen;
END_RCPP
}
// rmultinomial
Rcpp::IntegerVector rmultinomial(int N, Rcpp::NumericVector probs);
RcppExport SEXP _RFID_rmultinomial(SEXP NSEXP, SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(rmultinomial(N, probs));
    return rcpp_result_gen;
END_RCPP
}
// update_ell
arma::cube update_ell(Rcpp::NumericMatrix Phi, Rcpp::NumericMatrix Theta, Rcpp::IntegerMatrix ELL);
RcppExport SEXP _RFID_update_ell(SEXP PhiSEXP, SEXP ThetaSEXP, SEXP ELLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Phi(PhiSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type ELL(ELLSEXP);
    rcpp_result_gen = Rcpp::wrap(update_ell(Phi, Theta, ELL));
    return rcpp_result_gen;
END_RCPP
}
// rand_normal
double rand_normal(double mean, double stddev);
RcppExport SEXP _RFID_rand_normal(SEXP meanSEXP, SEXP stddevSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type stddev(stddevSEXP);
    rcpp_result_gen = Rcpp::wrap(rand_normal(mean, stddev));
    return rcpp_result_gen;
END_RCPP
}
// CrtRng
int CrtRng(int n, double r);
RcppExport SEXP _RFID_CrtRng(SEXP nSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(CrtRng(n, r));
    return rcpp_result_gen;
END_RCPP
}
// update_ELL
NumericMatrix update_ELL(NumericMatrix PhiTheta, IntegerMatrix x);
RcppExport SEXP _RFID_update_ELL(SEXP PhiThetaSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type PhiTheta(PhiThetaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(update_ELL(PhiTheta, x));
    return rcpp_result_gen;
END_RCPP
}
// update_z
IntegerMatrix update_z(NumericVector alpha, NumericMatrix PhiTheta, IntegerMatrix Y);
RcppExport SEXP _RFID_update_z(SEXP alphaSEXP, SEXP PhiThetaSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PhiTheta(PhiThetaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(update_z(alpha, PhiTheta, Y));
    return rcpp_result_gen;
END_RCPP
}

// update_z_purity
IntegerMatrix update_z_purity(NumericVector alpha, NumericMatrix PhiTheta, IntegerMatrix Y,double purity);
RcppExport SEXP _RFID_update_z(SEXP alphaSEXP, SEXP PhiThetaSEXP, SEXP YSEXP,SEXP puritySEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type PhiTheta(PhiThetaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Y(YSEXP);
     Rcpp::traits::input_parameter< IntegerMatrix >::type purity(puritySEXP);
    
    rcpp_result_gen = Rcpp::wrap(update_z_purity(alpha, PhiTheta, Y,purity));
    return rcpp_result_gen;
    END_RCPP
}


// update_m
IntegerMatrix update_m(NumericVector alpha, IntegerMatrix z, double lambda);
RcppExport SEXP _RFID_update_m(SEXP alphaSEXP, SEXP zSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(update_m(alpha, z, lambda));
    return rcpp_result_gen;
END_RCPP
}
// update_alpha
NumericVector update_alpha(double g0, double h0, NumericMatrix m, NumericVector p);
RcppExport SEXP _RFID_update_alpha(SEXP g0SEXP, SEXP h0SEXP, SEXP mSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type g0(g0SEXP);
    Rcpp::traits::input_parameter< double >::type h0(h0SEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type m(mSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(update_alpha(g0, h0, m, p));
    return rcpp_result_gen;
END_RCPP
}
// update_Phi
NumericMatrix update_Phi(double eta, NumericMatrix ell_ik, NumericMatrix Theta, NumericVector p);
RcppExport SEXP _RFID_update_Phi(SEXP etaSEXP, SEXP ell_ikSEXP, SEXP ThetaSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ell_ik(ell_ikSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(update_Phi(eta, ell_ik, Theta, p));
    return rcpp_result_gen;
END_RCPP
}
// update_Theta
NumericMatrix update_Theta(NumericVector r, NumericVector c, NumericVector p, NumericMatrix ell_jk, NumericMatrix Phi);
RcppExport SEXP _RFID_update_Theta(SEXP rSEXP, SEXP cSEXP, SEXP pSEXP, SEXP ell_jkSEXP, SEXP PhiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type ell_jk(ell_jkSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Phi(PhiSEXP);
    rcpp_result_gen = Rcpp::wrap(update_Theta(r, c, p, ell_jk, Phi));
    return rcpp_result_gen;
END_RCPP
}
// update_p
NumericVector update_p(double a0, double b0, NumericVector n, NumericMatrix Theta, NumericVector alpha);
RcppExport SEXP _RFID_update_p(SEXP a0SEXP, SEXP b0SEXP, SEXP nSEXP, SEXP ThetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Theta(ThetaSEXP);
    
    rcpp_result_gen = Rcpp::wrap(update_p(a0, b0, n, Theta));
    return rcpp_result_gen;
END_RCPP
}

// update_p_purity
NumericVector update_p_purity(double a0, double b0, NumericVector n, NumericMatrix Theta, NumericVector alpha);
RcppExport SEXP _RFID_update_p(SEXP a0SEXP, SEXP b0SEXP, SEXP nSEXP, SEXP ThetaSEXP) {
    BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Theta(ThetaSEXP);
    rcpp_result_gen = Rcpp::wrap(update_p(a0, b0, n, Theta));
    return rcpp_result_gen;
    END_RCPP
}


// update_c
NumericVector update_c(double e0, double f0, NumericVector r, NumericMatrix Theta);
RcppExport SEXP _RFID_update_c(SEXP e0SEXP, SEXP f0SEXP, SEXP rSEXP, SEXP ThetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< double >::type f0(f0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Theta(ThetaSEXP);
    rcpp_result_gen = Rcpp::wrap(update_c(e0, f0, r, Theta));
    return rcpp_result_gen;
END_RCPP
}
// update_gamma0
double update_gamma0(double a0, double b0, double p_tilde, double gamma0, NumericVector ell_k);
RcppExport SEXP _RFID_update_gamma0(SEXP a0SEXP, SEXP b0SEXP, SEXP p_tildeSEXP, SEXP gamma0SEXP, SEXP ell_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a0(a0SEXP);
    Rcpp::traits::input_parameter< double >::type b0(b0SEXP);
    Rcpp::traits::input_parameter< double >::type p_tilde(p_tildeSEXP);
    Rcpp::traits::input_parameter< double >::type gamma0(gamma0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ell_k(ell_kSEXP);
    rcpp_result_gen = Rcpp::wrap(update_gamma0(a0, b0, p_tilde, gamma0, ell_k));
    return rcpp_result_gen;
END_RCPP
}
// update_c0
double update_c0(double e0, double f0, double gamma0, NumericVector r);
RcppExport SEXP _RFID_update_c0(SEXP e0SEXP, SEXP f0SEXP, SEXP gamma0SEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e0(e0SEXP);
    Rcpp::traits::input_parameter< double >::type f0(f0SEXP);
    Rcpp::traits::input_parameter< double >::type gamma0(gamma0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(update_c0(e0, f0, gamma0, r));
    return rcpp_result_gen;
END_RCPP
}
// update_ell_tilde
NumericMatrix update_ell_tilde(NumericMatrix ell_jk, NumericVector r);
RcppExport SEXP _RFID_update_ell_tilde(SEXP ell_jkSEXP, SEXP rSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type ell_jk(ell_jkSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type r(rSEXP);
    rcpp_result_gen = Rcpp::wrap(update_ell_tilde(ell_jk, r));
    return rcpp_result_gen;
END_RCPP
}
// update_r
NumericVector update_r(double c0, double gamma0, NumericVector p_tilde_j, NumericVector ell_tilde_k);
RcppExport SEXP _RFID_update_r(SEXP c0SEXP, SEXP gamma0SEXP, SEXP p_tilde_jSEXP, SEXP ell_tilde_kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type c0(c0SEXP);
    Rcpp::traits::input_parameter< double >::type gamma0(gamma0SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p_tilde_j(p_tilde_jSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ell_tilde_k(ell_tilde_kSEXP);
    rcpp_result_gen = Rcpp::wrap(update_r(c0, gamma0, p_tilde_j, ell_tilde_k));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RFID_log_likelihood_Y", (DL_FUNC) &_RFID_log_likelihood_Y, 5},
    {"_RFID_likelihood_mat", (DL_FUNC) &_RFID_likelihood_mat, 5},
    {"_RFID_log_likelihood_alpha", (DL_FUNC) &_RFID_log_likelihood_alpha, 3},
    {"_RFID_log_likelihood_Phi", (DL_FUNC) &_RFID_log_likelihood_Phi, 2},
    {"_RFID_log_likelihood_Theta", (DL_FUNC) &_RFID_log_likelihood_Theta, 3},
    {"_RFID_rmultinomial", (DL_FUNC) &_RFID_rmultinomial, 2},
    {"_RFID_update_ell", (DL_FUNC) &_RFID_update_ell, 3},
    {"_RFID_rand_normal", (DL_FUNC) &_RFID_rand_normal, 2},
    {"_RFID_CrtRng", (DL_FUNC) &_RFID_CrtRng, 2},
    {"_RFID_update_ELL", (DL_FUNC) &_RFID_update_ELL, 2},
    {"_RFID_update_z", (DL_FUNC) &_RFID_update_z, 3},
    {"_RFID_update_m", (DL_FUNC) &_RFID_update_m, 3},
    {"_RFID_update_alpha", (DL_FUNC) &_RFID_update_alpha, 4},
    {"_RFID_update_Phi", (DL_FUNC) &_RFID_update_Phi, 4},
    {"_RFID_update_Theta", (DL_FUNC) &_RFID_update_Theta, 5},
    {"_RFID_update_p", (DL_FUNC) &_RFID_update_p, 5},
    {"_RFID_update_c", (DL_FUNC) &_RFID_update_c, 4},
    {"_RFID_update_gamma0", (DL_FUNC) &_RFID_update_gamma0, 5},
    {"_RFID_update_c0", (DL_FUNC) &_RFID_update_c0, 4},
    {"_RFID_update_ell_tilde", (DL_FUNC) &_RFID_update_ell_tilde, 2},
    {"_RFID_update_r", (DL_FUNC) &_RFID_update_r, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_RFID(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
