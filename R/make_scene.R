#' Make a full 3D scene
#'
#' @description Sample trees for each progenies from the models provided, write parameters, call
#' VPalm to create OPF files, and write the OPS files to disk for each progeny.
#'
#' @param data    The output from [compute_archi()]
#' @param nleaves    The numbers of leaves desired for the 3D palms
#' @param Progeny    The progeny desired. If `NULL` (the default), uses all progenies from `data`.
#' @param path       The path for writing the scene
#' @param AMAPStudio The path to AMAPStudio were VPalm lives
#' @param plot_design The design of the plot if custom. If `NULL`, calls [design_plot()]
#' with a quincunx disposition
#' @param plant_dist The distance between palm trees, used as y_dist in [design_plot()]
#' @param seed  The seed for random parameter generation (must be length ntrees)
#' @param overwrite Boolean. Should pre-existing scene files overwriten ?
#' @param ntrees  The number of trees to be sampled for each progeny if genetic variability
#' is required. If `NULL` or `0`, uses only average palm trees.
#'
#' @note To extract only average trees from progenies, simply set `ntrees= 0`.
#'
#' @return Writes a full 3D scene with a list of OPS for each progeny and a list of OPF for each plant
#' of each progeny. Returns a list of two:
#' * A ggplot of the design
#' * The OPF values for each progeny and for an average plant, as [extract_progeny()] do.
#' @export
#'
#' @examples
#' \dontrun{
#' # Initializing:
#' ntrees= 0 ; nleaves= 45 ; map= 47
#'
#' # Importing the inputs, and fitting the models:
#' Palm_Param= compute_archi(map = map, data_path = "1-Data/Archi", write_path = "3-Outputs")
#'
#' # Writing the scene:
#' scene= make_scene(data = Palm_Param, ntrees = ntrees, nleaves = nleaves, path = "3-Outputs",
#'                   AMAPStudio = "../../AMAPSTUDIO/VPalm", plant_dist = 9.2)
#' }
make_scene= function(data, nleaves= 45, Progeny= NULL,
                     path, AMAPStudio,
                     plant_dist= 9.2, plot_design= NULL,
                     seed= NULL,overwrite= T, ntrees= NULL){

  if(is.null(ntrees)){ntrees= 0}
  # Formating the VPalm outputs
  VPalm_list= Vpalmr::extract_progeny(data= data$input, model= data$model,
                                      n= ntrees, leaves= nleaves, seed= seed,
                                      average = T)
  if(!is.null(Progeny)){
    progenies= names(VPalm_list)
    Progeny= match.arg(Progeny, progenies, several.ok = TRUE)
    VPalm_list= VPalm_list[match(Progeny,names(VPalm_list))]
  }

  VPalm_in= format_progeny(data = VPalm_list)

  # Write the plants architectural parameters
  params= write_progeny(data = VPalm_in, path = file.path(path, "VPalm_inputs"),
                        verbose = F, overwrite = overwrite)
  if(all(params)){
    message("All VPalm parameters files successfully written in: ", file.path(path, "1-VPalm_inputs"))
  }else{
    stop("Error during VPalm parameter files writing")
  }


  # Make the OPFs (carefull, OPFs must be in a subfolder of OPS for AMAPStudio) :
  OPFs= make_opf_all(parameter = file.path(path, "VPalm_inputs"),
                     opf = file.path(path, "scenes/opf"),
                     AMAPStudio = AMAPStudio, overwrite = overwrite)

  # Design the planting pattern:
  if(is.null(plot_design)){
    if(ntrees>0){
      plot_design=
        design_plot(ntrees = ntrees, x0 = 0, y_dist = plant_dist)$design
    }else{
      plot_design=
        design_plot(ntrees = ntrees, x0 = 0, y_dist = plant_dist)$design
    }
  }

  # plot:
  design_ggplot=
    plot_design%>%
    ggplot2::ggplot(ggplot2::aes(x= x, y= y, color= Border))+
    ggplot2::geom_point()+
    ggplot2::ylim(low= unique(plot_design$ymin),
                  high= unique(plot_design$ymax))+
    ggplot2::xlim(low= unique(plot_design$xmin),
                  high= unique(plot_design$xmax))

  # Make the OPS:
  OPSs= make_ops_all(Progeny = as.list(names(VPalm_in[-grep("Average",names(VPalm_in))])),
                     design = list(plot_design), map = list(map),
                     path = file.path(path,"scenes"), overwrite = overwrite,
                     average = ifelse(ntrees>0,FALSE,TRUE))

  message("Scene successfully created")

  invisible(list(ggplot_design= design_ggplot, VPalm_parameters= VPalm_list))
}
