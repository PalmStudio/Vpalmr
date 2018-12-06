#' Compute Architectural parameters
#'
#' @description Imports the data and fit the models to compute the architectural
#' parameters for VPalm.
#'
#' @param map        Palm trees age desired in months after planting
#' @param data_path  Path to the folder were all data files resides
#' @param write_path Path were to write the function outputs. If `NULL` (the default),
#' the data is not written.
#'
#' @note Uses [import_data()] and [mod_all()] under the hood
#'
#' @return Returns invisibly the input data and the models fit as a list for each progeny.
#' @export
compute_archi= function(map, data_path, write_path= NULL){
  # Importing all data used to fit the models:
  Inputs= import_data(parameter= file.path(data_path,'ParameterSimu.csv'),
                      development= file.path(data_path,'Development_Rep4_SMSE.csv'),
                      phylotaxy= file.path(data_path,'Stem_SMSE14.csv'),
                      declination= file.path(data_path,'AnglesC&A_SMSE_Nov14.csv'),
                      curvature= file.path(data_path,'LeafCurvature_SMSE14.csv'),
                      toricity= file.path(data_path,'AnglesC&A_SMSE_Nov14.csv'),
                      leaf_area= file.path(data_path,'LeafArea_monitoring_SMSE.csv'),
                      axial_angle= file.path(data_path,'LeafDispositionComp_SMSE14.csv'),
                      petiole_width= file.path(data_path,'Petiole_SMSE14.csv'),
                      twist= file.path(data_path,'Torsion_SMSE14.csv'),
                      map = map)

  # Fit the models on data:
  mods= mod_all(x= Inputs)
  if(!is.null(write_path)){
    write_models(data = list(input= Inputs, model= mods), path = write_path)
    message("Data and models were successfully written in: ", write_path)
  }
  invisible(list(input= Inputs, model= mods))
}
