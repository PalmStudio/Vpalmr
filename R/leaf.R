
#' Curvature
#'
#' @description Compute the angles of points from relative positions and curvature
#'
#' @param position The positions along the axis relative to the tip (A point)
#' @param coefCurv Curvature coefficient
#'
#' @return The curvature angles for each position
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' curvature(position= 0.2, coefCurv= 0.5)
#'}
#' @export
curvature=function(position,coefCurv){
  ((1+ coefCurv)*(position**2))/(1+ coefCurv*(position**2))
}


#' Leaf curvature
#'
#' @description Compute the leaf curvature as VPalm do
#'
#' @param angC     The angle of the leaf at C point
#' @param angA     The angle of the leaf at A point (tip)
#' @param position The positions along the rachis relative to the tip (A point)
#' @param Length   The total rachis length
#' @param coefCurv Curvature coefficient
#'
#' @return The leaf curvature as a data.frame of the X and Y coordinates of the input relative positions
#'
#' @importFrom dplyr mutate transmute "%>%"
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' leaf_curvature(position= c(0,0.43,0.69,1),
#'                angC= c(27.14,27.14,27.14,27.14),
#'                angA= c(66.59,66.59,66.59,66.59),
#'                coefCurv= 3.58,
#'                Lenght= c(320.11,320.11,320.11,320.11))
#'}
#' @export
leaf_curvature=function(position,angC,angA,coefCurv,Length){
  nbSegment= length(position)-1
  angleSimu=
    angC*pi/180 +
    curvature(position= position, coefCurv= coefCurv)*
    (angA*pi/180-angC*pi/180)

  data.frame(position,angleSimu=angleSimu,segment=Length/nbSegment)%>%
    mutate(xSimu= segment*sin(lag(angleSimu,default=0)),
           ySimu= segment*cos(lag(angleSimu)))%>%
    mutate(xSimu= replace(xSimu,1,0),
           ySimu= replace(ySimu,1,0))%>%
    transmute(X= cumsum(xSimu),
              Y= cumsum(ySimu))
}


#' Leaf curvature error
#'
#' @description Compute the leaf curvature as VPalm do, compare
#' the output coordinates to observations, and return the error. This function is
#' generally used for parameter optimization
#'
#' @param angC     The angle of the leaf at C point
#' @param position The position of each point along the rachis relative to the tip (A point)
#' @param Length   The total rachis length
#' @param X_pos    The observed X positions in centimeters
#' @param Y_pos    The observed Y positions in centimeters
#' @param coefCurv Curvature coefficient
#' @param decAInfl The inflection point of the sigmoidal curve that predict point A
#'
#' @details The error is computed as the sum of the summed squarred root errors
#' of X and Y coordinates
#'
#' @note The function must be applied to each leaf independently
#'
#' @return The error of prediction of the leaf curvature
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' curvature_error(angC,position,Length,coefCurv,decAInfl)
#'}
#' @export
curvature_error= function(angC,position,Length,X_pos,Y_pos,coefCurv,decAInfl){
  angA= sigmoid(X= angC, max= 160, slope= 0.02,infl= decAInfl)

  leaf_curvature(position= position,
                  angC= angC, angA= angA, coefCurv= coefCurv,
                  Length= Length)%>%
    mutate(dmin= sqrt((X_pos-X)**2 +
                        (Y_pos-Y)**2))%>%
    summarise(dmin= sum(dmin))%>%.$dmin
}



#' Optimize leaf curvature parameters
#'
#' @description Optimize coefCurv and decAInfl, two leaf curvature parameters from VPalm
#'  using field data
#'
#' @param data     The input data.frame
#'
#' @details The input data.frame must have the following columns:
#' * Rank: the leaf rank. If not available, at least an index or a grouping factor for each leaf.
#' * angC: the leaf C angle
#' * RelativePositionRachisEstim: the relative position of the given points in the leaf compared
#' to the tip.
#' * rachisLength: total length of the leaf
#' * X_distance_cm: X coordinates in centimeters
#' * Y_distance_cm: Y coordinates in centimeters
#'
#' @return A data.frame with coefCurv and decAInfl values, the final error value returned by
#' [curvature_error()] and the convergence value from [stats::optim()].
#'
#' @importFrom dplyr group_by summarise "%>%" ungroup
#' @importFrom stats optim
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' optim_Leaf_Curv_Infl(data)
#'}
#' @export
# Function to apply to each leaf independently:
optim_Leaf_Curv_Infl= function(data){
  nrows= nrow(data)
  optimDist=
    optim(par= list(coefCurv=0.5,decAInfl=30),
          fn= function(par){
            data%>%
              group_by(Rank)%>%
              summarise(curv= curvature_error(angC= angC, position = RelativePositionRachisEstim,
                                              Length = rachisLength, X_pos= X_distance_cm,
                                              Y_pos= Y_distance_cm, coefCurv=par[1],
                                              decAInfl=par[2]))%>%
              ungroup()%>%
              summarise(curv= sum(curv))%>%as.data.frame()
          },
          method="L-BFGS-B")
  data.frame(decAInfl= optimDist$par[["decAInfl"]],
             coefCurv= optimDist$par[["coefCurv"]],
             value= optimDist$value/nrows,
             conv=optimDist$convergence)
}



#' Leaf deviation
#'
#' @description Compute the leaf deviation on point A
#'
#' @param devC   Deviation on each point compared to C point
#' @param Length The total rachis length
#'
#' @return The leaf deviation in A point
#'
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' deviation_A(devC= 86, Lenght= 320)
#'}
#' @export
deviation_A=function(devC,Length){
  # NB: leaf deviation in Z is like leaf curvature in X
  angles= 0:90
  devCm= data.frame(devA= angles, devC= devC)
  devCm$devCm=
    sapply(devCm$devA, function(x){
      leaf_curvature(position = seq(0,1,0.01), angC = 0,
                     angA = x,coefCurv = 2,
                     Length = Length)$X%>%tail(.,1)
    })

  devCm$dif= abs(devCm$devCm-devCm$devC)
  devA= devCm[devCm$dif==min(devCm$dif),]$devA

  return(devA)
}

#' Rachis section height
#'
#' @description Computes the leaf section height from its
#' relative position on the leaf and the \code{a} parameter.
#'
#' @param Position Relative position on the rachis (0-1)
#' @param a Parameter
#'
#' @return The section height
#' @export
section_height= function(Position,a){
  1 + a*(Position**3)
}
