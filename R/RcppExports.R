# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

GWAS_approx_pql_bed <- function(pA, PY, P, p, beg, end) {
    .Call('_milorGWAS_GWAS_approx_pql_bed', PACKAGE = 'milorGWAS', pA, PY, P, p, beg, end)
}

GWAS_logit_offset_bed <- function(pA, p, Y, Offset, Q, beg, end, tol, max_iter) {
    .Call('_milorGWAS_GWAS_logit_offset_bed', PACKAGE = 'milorGWAS', pA, p, Y, Offset, Q, beg, end, tol, max_iter)
}

min_ <- function(x) {
    .Call('_milorGWAS_min_', PACKAGE = 'milorGWAS', x)
}

max_ <- function(x) {
    .Call('_milorGWAS_max_', PACKAGE = 'milorGWAS', x)
}

manhattan_thinning <- function(x, y, mx, my) {
    .Call('_milorGWAS_manhattan_thinning', PACKAGE = 'milorGWAS', x, y, mx, my)
}

