####leaflets dicposition on rachis



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
