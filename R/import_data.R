


#' Sigmoidal function
#'
#' @description Estimate a sigmoidal function from parameters.
#'
#' @param file_archi Path to the architectural file
#' @param file_la    Path to the leaf area file
#' @param file_phylo Path to the phylotaxy file
#' @param file_angles Path to the leaf declination angles file
#'
#' @return A sigmoid
#'
#' @importFrom dplyr ungroup group_by summarise "%>%"
#' @importFrom data.table fread
#'
#' @export
#'
import_data= function(file_archi,file_la,file_phylo,file_angles){

  # File input with date and number of leaves emitted
  Parameter=
    fread('1-Data/Archi/ParameterSimu.csv', data.table = F)%>%
    mutate(Date= as.Date(Date,format='%d/%m/%Y'))



  return(list(Parameter,Dec,Curve))
}
