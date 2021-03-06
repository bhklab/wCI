# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

usable <- function(x1, x2, delta) {
    .Call('_wCI_usable', PACKAGE = 'wCI', x1, x2, delta)
}

usableHard <- function(x1, x2, delta) {
    .Call('_wCI_usableHard', PACKAGE = 'wCI', x1, x2, delta)
}

concordanceIndex_modified_helper <- function(x, y, deltaX, deltaY, alpha, outx, alternative, logicOp) {
    .Call('_wCI_concordanceIndex_modified_helper', PACKAGE = 'wCI', x, y, deltaX, deltaY, alpha, outx, alternative, logicOp)
}

kernel_gaussian_C <- function(x, m, s) {
    .Call('_wCI_kernel_gaussian_C', PACKAGE = 'wCI', x, m, s)
}

kernel_laplace_C <- function(x, m, b) {
    .Call('_wCI_kernel_laplace_C', PACKAGE = 'wCI', x, m, b)
}

concordanceIndex_modified_helper_weighted <- function(x, y, deltaX, deltaY, weightingFun_pred, weightingFun_obs, alpha, outx, alternative, logicOp, max_weight, max_weight_obs) {
    .Call('_wCI_concordanceIndex_modified_helper_weighted', PACKAGE = 'wCI', x, y, deltaX, deltaY, weightingFun_pred, weightingFun_obs, alpha, outx, alternative, logicOp, max_weight, max_weight_obs)
}

do_swaps_c <- function(pvector, pvlen, pnswap, pseed) {
    .Call('_wCI_do_swaps_c', PACKAGE = 'wCI', pvector, pvlen, pnswap, pseed)
}

do_random_swaps_c <- function(pvector, pvlen, pnswap, pseed) {
    .Call('_wCI_do_random_swaps_c', PACKAGE = 'wCI', pvector, pvlen, pnswap, pseed)
}

rCIPermC <- function(pin_x, pin_y, pobsCI, pR, pB, pn, pxties, pyties, palternative, pseed) {
    .Call('_wCI_rCIPermC', PACKAGE = 'wCI', pin_x, pin_y, pobsCI, pR, pB, pn, pxties, pyties, palternative, pseed)
}

rCIBootC <- function(prcimat, pR, pn, pxties, pyties, pseed) {
    .Call('_wCI_rCIBootC', PACKAGE = 'wCI', prcimat, pR, pn, pxties, pyties, pseed)
}

newPCI <- function(pin_x, pin_y, pn, pxties, pyties, pdeltaX, pdeltaY, plogic) {
    .Call('_wCI_newPCI', PACKAGE = 'wCI', pin_x, pin_y, pn, pxties, pyties, pdeltaX, pdeltaY, plogic)
}

frCI <- function(pin_x, pin_y, pdelta_x, pdelta_y, pn, pdiscard_x_ties, pdiscard_y_ties) {
    .Call('_wCI_frCI', PACKAGE = 'wCI', pin_x, pin_y, pdelta_x, pdelta_y, pn, pdiscard_x_ties, pdiscard_y_ties)
}

KCI <- function(pin_x, pin_y, pn, psymmetric, pkern1, pkern2) {
    .Call('_wCI_KCI', PACKAGE = 'wCI', pin_x, pin_y, pn, psymmetric, pkern1, pkern2)
}

