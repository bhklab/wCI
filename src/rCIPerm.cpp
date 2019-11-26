#include <Rcpp.h>

extern "C" SEXP permC(SEXP pin_x,
             SEXP pin_y,
             SEXP pobsCI,
             SEXP pR,
             SEXP pB,
             SEXP pn,
             SEXP pxties,
             SEXP pyties, 
             SEXP pseed);

// [[Rcpp::export]]
SEXP rCIPermC(SEXP pin_x, SEXP pin_y, SEXP pobsCI, SEXP pR, SEXP pB, SEXP pn, SEXP pxties, SEXP pyties,  SEXP pseed) {
    return permC(pin_x,
             pin_y,
             pobsCI,
             pR,
             pB,
             pn,
             pxties,
             pyties, 
             pseed);
}