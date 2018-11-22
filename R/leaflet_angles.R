
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


#' Axial insertion angle
#'
#' @description Computes the leaflet axial insertion angle
#'
#' @param position Relative position of the leaflet
#' @param angle_C  Angle of the leaflet compared to C point
#' @param slope_C  Slope of the leaflet to C point relationship
#' @param angle_A  Angle of the leaflet compared to A point
#'
#' @return The leaflet axial angle
#' @export
leaflet_axial_angle=function(position,angle_C,slope_C,angle_A){
  a= angle_C**2
  b= slope_C*2*sqrt(a)
  c= angle_A**2-a-b
  sqrt(a+b*position+c*(position**3))
}



#' Radial insertion angle
#'
#' @description Computes the leaflet radial insertion angle
#'
#' @param position Relative position of the leaflet
#' @param A0   Intercept
#' @param Amax Maximum angle
#' @param Xm   Maximum position value
#'
#' @return The leaflet axial angle
#' @export
leaflet_radial_angle=function(position,A0,Amax,Xm){
  c1=(A0-Amax)/(Xm**2)
  b1=-2*c1*Xm
  a1=A0
  c2=-Amax/((Xm-1)**2)
  b2=-2*c2*Xm
  a2=-b2-c2

  ifelse(position < Xm,
         a1+b1*position+c1*(position**2),
         a2+b2*position+c2*(position**2))
}
