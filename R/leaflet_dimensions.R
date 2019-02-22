
#' Leaflet dimensions
#'
#' @description Computes the leaflet length or maximum width using leaflet
#' position along the leaf and several parameters
#'
#' @param X      Relative position along the leaf
#' @param Ymax   Maximal Y
#' @param Y0     Intercept
#' @param Yfin   Value at the maximal X
#' @param X_Ymax X value at Ymax
#'
#' @aliases leaflet_max_width
#'
#' @return The leaflet length or maximum width
#' @export
leaflet_length=function(X,Ymax,Y0,Yfin,X_Ymax){
  ifelse(X<X_Ymax,
         Y0+X*2*(Ymax-Y0)/X_Ymax+X^2*(Y0-Ymax)/X_Ymax^2,
         (Yfin-Ymax)/((1-X_Ymax)^2)*(X^2-2*X_Ymax*X+X_Ymax^2)+Ymax)
}

#' @rdname leaflet_length
#' @export
leaflet_max_width=function(X,Ymax, Y0, Yfin, X_Ymax){
  ifelse(X<X_Ymax,
         Y0+X*(Ymax-Y0)/X_Ymax,
         X*(Yfin-Ymax)/(1-X_Ymax)+(Ymax-X_Ymax*Yfin)/(1-X_Ymax))
}


#' Fit nls on leaflet shape
#'
#' @description Find the coefficients of the beta function for leaflet shape.
#'
#'
#' @param width    Leaflet relative width (compared to maximum observed on the leaf, 0-1)
#' @param position Leaflet relative position on the leaf (0-1)
#'
#' @details The function determines the nonlinear (weighted) least-squares estimates of the
#' xm and ym parameters using the \code{\link[stats]{nls}} of the following model:
#' \code{width~shape_beta(x= position,xm= xm,ym= ym)}.
#'  See \code{\link{shape_beta}} for more informations.
#'
#' @return The estimation of the xm and ym parameters from the nls fit
#' @importFrom stats nls
#' @export
leaflet_shape_nls= function(width, position){
  ajust=
    try(stats::nls(width~shape_beta(x=position,xm= xm_estim,ym=ym_estim),
                   start=list(xm_estim=0.5,ym_estim=0.5),
                   control=list(maxiter=5000000,minFactor=0.000000000001,warnOnly=T),
                   trace=F),
        silent=T)
  if(!inherits(ajust,'try-error')){
    out= c(xm= stats::coef(ajust)[['xm_estim']],
           ym= stats::coef(ajust)[['ym_estim']])
  }else{
    out= c(NA,NA)
  }
  out
}

#' Beta function for leaflet shape
#'
#' @description The function computes the relative width of the leaflet at all x
#'  positions. The outputed width is relative to the leaflet max width.
#'
#' @param x  Relative position on the leaflet
#' @param xm Position of the maximum leaflet width (0-1)
#' @param ym Shape factor
#'
#' @return The relative width of the leaflet at x position
#' @export
#'
shape_beta=function(x,xm,ym){
  q= ((1-xm)*log(ym*xm*(1-xm))+(2*xm-1)*log(xm))/(xm*log(xm)+(1-xm)*log(1-xm))
  p= (1-2*xm+q*xm)/(1-xm)

  (1/ym)*(x**(p-1)*(1-x)**(q-1))
}
