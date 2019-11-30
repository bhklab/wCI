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

extern "C" SEXP newPairedConcIndex(SEXP pin_x,
             SEXP pin_y,
             SEXP pn,
             SEXP pxties,
             SEXP pyties,
             SEXP pdeltaX,
             SEXP pdeltaY,
             SEXP plogic);

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

// [[Rcpp::export]]
SEXP newPCI(SEXP pin_x,
             SEXP pin_y,
             SEXP pn,
             SEXP pxties,
             SEXP pyties,
             SEXP pdeltaX,
             SEXP pdeltaY,
             SEXP plogic) {
    return newPairedConcIndex(pin_x,
             pin_y,
             pn,
             pxties,
             pyties,
             pdeltaX,
             pdeltaY,
             plogic);
}
