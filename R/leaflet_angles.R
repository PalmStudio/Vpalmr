
#' Compute the leaflet angles
#'
#' @description Estimate the leaflet angles at point A and C using the tangents along the leaf curvature
#'
#' @param x      The values in the x axis
#' @param y      The values in the y axis
#' @param res    The re-sampling resolution, with the same unit as inputs.
#' @param method The smoothing method, either using [stats::loess()] or [stats::smooth.spline()]
#' @param ...    Further arguments to pass to the smoothing function
#'
#' @details The x and y coordinates must be ordered from first to last point on the object.
#'  The function is used to compute the tangents of points for a curved object. The function
#'  re-sample the object using a smoothing function to give the closest tangeant to the point.
#'  The higher the resolution, the higher the sub-sampling, and the closer the tangent will be.
#'  Here is how the function works:
#'  1. The user sample some points from a curved object and register their X and Y positions
#'     in a euclidean coordinate system, and give the points coordinates as input of the function
#'  1. An approximation of the object curvature is made using a smoothing function on all
#'  measured points
#'  1. The object is re-sampled using the resulting smoothing function with a
#'  sampling resolution of `res`.
#'  1. The tangent is computed for each input point using the closer re-sampled points
#'
#' @return The tangent for each point of a curved object
#' @seealso [object_lenght()]
#' @importFrom dplyr lead "%>%"
#' @importFrom stats smooth.spline loess
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' object_tans(x= 1:10,y= rep(1,10))
#'}
#' @export
find_angle= function(x,y){
  angles= object_tans(x,y,method= "smooth.spline")$angle
  angles[angles<0]= pi + angles[angles<0]
  angles= angles*180/pi
  angle_C= angles[1]
  angle_A= tail(angles,1)
  return(c(angle_C,angle_A))
}
