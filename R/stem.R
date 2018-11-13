#' Stem height
#'
#' @description Computes the stem height from the total number of leaves emitted during
#' the palm life cycle
#'
#' @param X    Total number of leaves emmitted throughout the palm life
#' @param y0   Minimum height
#' @param coef Growth coefficient
#'
#' @details The stem height is computing following this equation:
#' \deqn{Stemheight=y_0\cdot e^{coef\cdot X}}{Stem_height= y0*exp(coef*X)}
#'
#' @return The stem height
#' @export
Stem_height=function(X,y0,coef){
  y0*exp(coef*X)
}
