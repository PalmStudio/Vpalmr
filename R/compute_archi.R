#' Compute Architectural parameters
#'
#' @description Imports the data and fit the models to compute the architectural
#' parameters for VPalm.
#'
#' @param map        Palm trees age desired in months after planting
#' @param data_path  The folder path, or the files path (see details)
#' @param write_path Path were to write the function outputs. If `NULL` (the default),
#' the data is not written.
#'
#' @details If `data_path` is a character vector of length one, it is used as the path to the folder where all data files are,
#' and all files are read from this folder using the default file names. If it is a named list (the defaukt), it is used as
#' the path for each file.
#'
#' @note Uses [import_data()] and [mod_all()] under the hood
#'
#' @return Returns invisibly the input data and the models fit as a list for each progeny.
#'
#' @examples
#' \dontrun{
#' # Using the data_path as a folder path with default file names
#' Palm_Param= compute_archi(map = 47, data_path = "1-Data/Archi", write_path = "2-Outputs")
#'
#' # Or using it with custom file names:
#' Palm_Param= compute_archi(map = 47, data_path = list(parameter= '1-Data/Parameter.csv',
#'                                                      development= '1-Data/Development.csv',
#'                                                      phylotaxy= '1-Data/Stem.csv',
#'                                                      declination= '1-Data/AnglesC&A.csv',
#'                                                      curvature= '1-Data/LeafCurvature.csv',
#'                                                      torsion= '1-Data/AnglesC&A.csv',
#'                                                      leaf_area= '1-Data/LeafArea.csv',
#'                                                      axial_angle= '1-Data/LeafDisposition.csv',
#'                                                      petiole_width= '1-Data/Petiole.csv',
#'                                                      twist= '1-Data/Torsion.csv'),
#'                           write_path = "2-Outputs")
#' }
#'
#' @export
compute_archi= function(map, data_path= list(parameter= 'ParameterSimu.csv',
                                             development= 'Development_Rep4_SMSE.csv',
                                             phylotaxy= 'Stem_SMSE14.csv',
                                             declination= 'AnglesC&A_SMSE_Nov14.csv',
                                             curvature= 'LeafCurvature_SMSE14.csv',
                                             torsion= 'AnglesC&A_SMSE_Nov14.csv',
                                             leaf_area= 'LeafArea_monitoring_SMSE.csv',
                                             axial_angle= 'LeafDispositionComp_SMSE14.csv',
                                             petiole_width= 'Petiole_SMSE14.csv',
                                             twist= 'Torsion_SMSE14.csv'),
                        write_path= NULL){


  if(length(data_path)==1){
    data_path=
      list(
        parameter= file.path(data_path,'ParameterSimu.csv'),
        development= file.path(data_path,'Development_Rep4_SMSE.csv'),
        phylotaxy= file.path(data_path,'Stem_SMSE14.csv'),
        declination= file.path(data_path,'AnglesC&A_SMSE_Nov14.csv'),
        curvature= file.path(data_path,'LeafCurvature_SMSE14.csv'),
        torsion= file.path(data_path,'AnglesC&A_SMSE_Nov14.csv'),
        leaf_area= file.path(data_path,'LeafArea_monitoring_SMSE.csv'),
        axial_angle= file.path(data_path,'LeafDispositionComp_SMSE14.csv'),
        petiole_width= file.path(data_path,'Petiole_SMSE14.csv'),
        twist= file.path(data_path,'Torsion_SMSE14.csv')
      )
  }
  data_path$map= map

  # Importing all data used to fit the models:
  Inputs= do.call(import_data, data_path)

  # Fit the models on data:
  mods= mod_all(x= Inputs)
  if(!is.null(write_path)){
    if(!dir.exists(write_path)){
      dir.create(write_path, recursive = TRUE)
    }
    write_models(data = list(input= Inputs, model= mods), path = write_path)
    message("Data and models were successfully written in: ", write_path)
  }
  invisible(list(input= Inputs, model= mods))
}
