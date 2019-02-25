#' Sigmoidal function
#'
#' @description Estimate a sigmoidal function from parameters.
#'
#' @param X       The x values
#' @param max     The y maximum value
#' @param slope   The slope at the inflection point
#' @param infl    The X position of the inflection point
#'
#' @return A sigmoid
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' sigmoid(X= 1/10,max= 3,slope= 1,infl=5)
#'}
#' @export
sigmoid=function(X,max,slope,infl){
  max/(1+exp(4* slope*(infl-X)))
}


#' Distance bewteen 2D points
#'
#' @description Estimate the linear distance between points in 2D coordinates
#'
#' @param x The values in the x axis
#' @param y The values in the y axis
#'
#' @details The function returns a 0 value for the first two points.
#'
#' @return The total distance (lenght) of the series of objects
#'
#' @importFrom dplyr lag "%>%"
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' distance_2D(x= 1:10,y= 1:10)
#'}
#' @export
distance_2D= function(x,y){
  sqrt((x-dplyr::lag(x,default = 0))^2+
         (y-dplyr::lag(y,default = 0))^2)
}


#' Lenght of non-linear objects
#'
#' @description Estimate the smoothed distance between points sampled
#' in an object
#'
#' @param x      The values in the x axis
#' @param y      The values in the y axis
#' @param res    The re-sampling resolution, with the same unit as inputs.
#' @param method The smoothing method, either using [stats::loess()] or [stats::smooth.spline()]
#' @param ...    Further arguments to pass to the smoothing function
#'
#' @details The x and y coordinates must be ordered from first to last point on the object.
#'  The function uses the c(0,0) point coordinate as the reference for the first point length.
#'  The function is used to compute the length of an object from subsampled points
#'  coordinates. Here is how the function works:
#'  1. The user sample some points from an object and register their X and Y positions
#'     in a euclidean coordinate system, and give the points coordinates as input to the function
#'  1. An approximation of the object curvature is made using a smoothing function on all
#'  measured points
#'  1. The object is re-sampled using the resulting smoothing function with a
#'  sampling resolution of `res`.
#'  1. The linear distances between each re-sampled points is computed
#'  1. The distances of each points bewteen the input coordinates are summed, so the estimation of
#'  the distances between the input coordinates are returned
#'
#' @note The distance from the 0 reference for each points can be computed using [base::cumsum()],
#'  and the total length of the object using [base::sum()].
#'  [stats::smooth.spline()] needs at least four unique values, so if the number of points is unknown,
#'  consider using [stats::loess()] instead.
#'
#' @return The distance between points from a curved object
#'
#' @importFrom dplyr lag "%>%"
#' @importFrom stats smooth.spline loess
#' @importFrom rlang .data
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' Segment_lengths= object_lenght(x= 1:10,y= rep(1,10))
#' Length_from_0= cumsum(Segment_lengths)
#' Total_length= sum(Segment_lengths)
#'}
#' @export
object_lenght= function(x,y,res = 1,method = c("loess","smooth.spline"),...){

  if(length(x)!=length(y)){stop("x and y must have the same length")}
  if(length(x)<4 & method=="smooth.spline"){
    warning(paste("x must have at least four distinct values for smooth.spline,",
                  "using linear interpolation instead"))
  }

  # Re-sampling with the given resolution
  x_pred= seq(0,max(x),res)%>%dplyr::union(x)

  if(method=="smooth.spline"){
    df_model=
      predict(object= smooth.spline(x,y,...),
              x= seq(0,max(x),res)%>%dplyr::union(x))%>%
      data.frame()%>%dplyr::rename(x_res= .data$x, y_res= .data$y)
  }else if(method=="loess"){
    df_model=
      stats::predict(object= stats::loess(y~x,...), x_pred)%>%
      data.frame(x_res= x_pred, y_res= .)
  }

  df_model%>%
    mutate(Distance= distance_2D(x = .data$x_res,y = .data$y_res),
           x= cut(.data$x_res, x, labels= x[-1]))%>%
    dplyr::slice(-1)%>%
    group_by(.data$x)%>%
    summarise(Distance= sum(.data$Distance))%>%
    merge(data.frame(x,y),., by="x", all.x= TRUE)%>%
    mutate(Distance= replace(.data$Distance,1,0))
}


#' Compute the tangants along non-linear objects
#'
#' @description Estimate the tangants for points from a curved object using smoothing
#'
#' @param x      The values in the x axis
#' @param y      The values in the y axis
#' @param res    The re-sampling resolution, with the same unit as inputs.
#' @param method The smoothing method, either using [stats::loess()] or [stats::smooth.spline()]
#' @param ...    Further arguments to pass to the smoothing function
#'
#' @details The x and y coordinates must be ordered from first to last point on the object.
#'  The function is used to compute the tangents of points for a curved object. The function
#'  re-sample the object using a smoothing function to give the closest tangant to the point.
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
#' @importFrom dplyr lead "%>%" n
#' @importFrom stats smooth.spline loess
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' object_tans(x= 1:10,y= rep(1,10))
#'}
#' @export
object_tans= function(x,y,res = 1, method = c("loess","smooth.spline"),...){

  if(length(x)!=length(y)){stop("x and y must have the same length")}
  if(length(x)<4 & method=="smooth.spline"){
    warning(paste("x must have at least four distinct values for smooth.spline,",
                  "using linear interpolation instead"))
  }

  # Re-sampling with the given resolution
  x_pred= seq(0,max(x),res)%>%dplyr::union(x)

  if(method=="smooth.spline"){
    df_model=
      stats::predict(object= smooth.spline(x,y,...),
                     x= seq(0,max(x),res)%>%dplyr::union(x))%>%
      data.frame()%>%
      dplyr::rename(x_res= .data$x, y_res= .data$y)
  }else if(method=="loess"){
    df_model=
      stats::predict(object= stats::loess(y~x,...),x_pred)%>%
      data.frame(x_res= x_pred, y_res= .)
  }

  df_model%>%
    dplyr::mutate(angle= atan((lead(.data$x_res)-.data$x_res)/
                                (lead(.data$y_res)-.data$y_res)))%>%
    dplyr::mutate(angle= replace(.data$angle,n(),lag(.data$angle)%>%.[n()]))%>%
    dplyr::filter(.data$x_res%in%x)
}



#' AlKashi
#'
#' @description Get the length of the 3 side of a triangle knowing the 2 other sides using
#' the law of cosines, also known as the generalized Pythagorean theorem.
#'
#' @param x lenght of one side
#' @param l lenght of the second side
#'
#' @return The length of the third side
#' @export
AlKashi= function(x,l){
  acos((l^2+l^2-x^2)/(2*l*l))*180/pi
}



#' warning construction
#'
#' @description Helper to add warnings and return only one
#'
#' @param x A previous warning message as a character string
#'
#' @return One warning from two: the warn object and x
#'
#' @keywords internal
warn_inc= function(warn,x){
  warn= paste(warn,"\n",x)
}


#' Parallel mean
#'
#' @description Returns the parallel mean, in the same fashion as [base::Extremes()]
#'
#' @param ...   Numeric arguments
#' @param na.rm A logical indicating whether missing values should be removed.
#'
#' @details The function uses a combination of [base::cbind()] and [rowMeans][base::colSums()].
#' The function can take one or more vectors or matrices as arguments and returns a single
#' 'parallel' mean of the vectors. The first element of the result is the mean of the first
#' elements of all the arguments and so on.
#'
#' @return A vector of the parallel mean of same length as the arguments.
#' @export
#'
#' @examples
#' pmean(1:10, 21:30)
pmean= function(..., na.rm = T){
  x= cbind(...)
  if(!is.null(x)){
    x= rowMeans(x, na.rm= na.rm)
  }else{
    x= NULL
  }
  x
}



#' Update progress bar
#'
#' @param prog A progress function
#' @param name The progress info (e.g. name of the function executed)
#'
#' @details This function is a helper that is called to increment a progress function
#' passed as an argument. It is used for example for a Shiny progress bar using
#'  [progress](https://shiny.rstudio.com/articles/progress.html).
#'
#' @return Nothing. Just update the prog function.
#' @export
#'
up_progress= function(prog,name){
  if(is.function(prog)){
    prog(name)
  }
}





