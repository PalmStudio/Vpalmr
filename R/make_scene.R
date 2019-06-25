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
#' @param plant_dist The distance between palm trees, used as y_dist in [design_plot()]. Used only
#' if `plot_design` is `NULL`.
#' @param seed  The seed for random parameter generation (see details)
#' @param overwrite Boolean. Should pre-existing scene files overwriten ?
#' @param ntrees  The number of trees to be sampled for each progeny if genetic variability
#' is required. If `NULL` or `0`, uses only average palm trees.
#' @param progress A progress function (see details) Shiny progress bar if ran using Shiny application.
#'
#' @note To extract only average trees from progenies, simply set `ntrees= 0`. To make a scene with the average
#' tree from all progenies, set `Progeny= "Average"`.
#'
#' @details `seed` should be a named list of length equal to the number of progenies, and each
#' object with n seeds each. The names of the list objects should match the ones of the
#' progenies. Alternatively, `seed` can be a vector of seeds of length n that will be recycled
#' for each progeny. If `NULL`, a random seed is generated for each tree of each progeny
#' using [base::sample()].
#' @section Progress bar:
#' As this function can take some time to run, it is possible to pass a
#' progress function to the `progress` argument as for Shiny progress bar using
#' [progress](https://shiny.rstudio.com/articles/progress.html). This function calls
#' [up_progress()] under the hood. There are 7 calls to the progress function in total.
#'
#' @return Writes a full 3D scene with a list of OPS for each progeny and a list of OPF for each plant
#' of each progeny. Returns a list of three:
#' * The OPF values for each progeny and for an average plant, as [extract_progeny()] do (VPalm_parameters).
#' * The data.frame of the plot design (plot_design)
#' * A ggplot of the design (ggplot_design)
#'
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
                     seed= NULL,overwrite= T, ntrees= NULL,
                     progress= NULL){

  path= normalizePath(path, winslash = "/", mustWork = FALSE)
  AMAPStudio= normalizePath(AMAPStudio, winslash = "/", mustWork = TRUE)

  # If the user input a jar file, only use the directory name:
  if(grepl('.jar',AMAPStudio)){
    AMAPStudio= dirname(AMAPStudio)
  }

  if(is.null(ntrees)){ntrees= 0}
  # Formating the VPalm outputs
  VPalm_list= Vpalmr::extract_progeny(data= data$input, model= data$model,
                                      n= ntrees, leaves= nleaves, seed= seed,
                                      average = T)
  up_progress(progress,'extract_progeny')

  map= unlist(VPalm_list)%>%.[grep("MAP_requested",names(.))]

  progenies= names(VPalm_list)
  if(is.null(Progeny)){
    Progeny= progenies
  }else{
    Progeny= match.arg(Progeny, progenies, several.ok = TRUE)
    VPalm_list= VPalm_list[match(Progeny,names(VPalm_list))]
    map= map[match(Progeny,names(VPalm_list))]
  }

  VPalm_in= format_progeny(data = VPalm_list)
  up_progress(progress,'format_progeny')

  # Write the plants architectural parameters
  params= write_progeny(data = VPalm_in, path = file.path(path, "VPalm_inputs"),
                        verbose = F, overwrite = overwrite)
  up_progress(progress,'write_progeny')

  if(all(params)){
    message("All VPalm parameters files successfully written in: ", file.path(path, "1-VPalm_inputs"))
  }else{
    stop("Error during VPalm parameter files writing")
  }


  # Make the OPFs (carefull, OPFs must be in a subfolder of OPS for AMAPStudio) :
  OPFs= make_opf_all(parameter = file.path(path, "VPalm_inputs"),
                     opf = file.path(path, "scenes/opf"),
                     AMAPStudio = AMAPStudio, overwrite = overwrite)
  up_progress(progress,'make_opf_all')

  # Design the planting pattern:
  if(is.null(plot_design)){
    if(ntrees>0){
      ro= floor(sqrt(ntrees/2))
      co= ceiling(sqrt(ntrees/2))
      if((ro*co*2)<ntrees){
        co= co+1
      }
      plot_design=
        design_plot(rows= ro,
                    cols= co, x0 = 0, x_dist = plant_dist)$design
      plot_design= plot_design[1:ntrees,]
    }else{
      plot_design=
        design_plot(rows=1, cols= 1, x0 = 0, x_dist = plant_dist)$design
    }
    up_progress(progress,'design_plot')
  }

  # plot:
  design_ggplot=
    plot_design%>%
    ggplot2::ggplot(ggplot2::aes(x= .data$x, y= .data$y, color= .data$Border))+
    ggplot2::geom_point()+
    ggplot2::ylim(low= min(unique(plot_design$ymin)),
                  high= max(unique(plot_design$ymax)))+
    ggplot2::xlim(low= min(unique(plot_design$xmin)),
                  high= max(unique(plot_design$xmax)))
  up_progress(progress,'ggplot of the design')

  # Make the OPS:
  OPSs= make_ops_all(Progeny = ifelse(Progeny=="Average","Average",
                                      as.list(names(VPalm_in)[
                                        names(VPalm_in)!="Average"])),
                     design = list(plot_design), map = as.list(map),
                     path = file.path(path,"scenes"), overwrite = overwrite,
                     average = ifelse(ntrees>0,FALSE,TRUE))
  up_progress(progress,'make_ops_all')

  message("Scene successfully created")

  invisible(list(VPalm_parameters = VPalm_list, plot_design = plot_design,
                 ggplot_design = design_ggplot))
}
