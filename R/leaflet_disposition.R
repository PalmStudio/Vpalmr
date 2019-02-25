#' Leaflet position
#'
#' @description Predict the leaflet position along the rachis
#' using the leaflet relative rank on rachis and a parameter
#'
#' @param Rank  Leaflet rank on rachis
#' @param param Distribution parameter
#'
#' @return The leaflet position on the rachis
#' @export
leaflet_position=function(Rank,param){
  (param*(Rank**2))/(1+(param-1)*(Rank**2))
}


#' Positions on leaflet
#'
#' @description Compute positions along the leaflet length
#'
#' @param x The (maximum) length of the leaflet
#'
#' @return A vector of 8 positions along the leaflet
#' @export
#'
#' @examples
#' position_on_leaflet(25)
position_on_leaflet= function(x){
  (c(0,1/x,seq(0.1,1,0.2),1)*x)%>%round()
}
