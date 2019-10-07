#' Design the plot
#'
#' @description Help designing a quincunx planting pattern of palm stand.
#'
#' @param cols   How many times the voronoi sub-plot is repeated in columns?
#' @param rows   How many times the voronoi sub-plot is repeated in rows?
#' @param x_dist The inter-row distance (m). See details.
#' @param y_dist The intra-row distance (m). See details.
#' @param x0     The minimum X coordinates
#'
#' @details If only one of \eqn{x_{dist}}{x_dist} or \eqn{y_{dist}}{y_dist} is given, the function
#'  will compute the other distance using the following equation:
#'  \deqn{\sqrt{z_{dist}^2-\left(\frac{z_{dist}}{2}\right)^2}}{sqrt(z_dist*z_dist-((z_dist/2)**2))}
#' with \eqn{z_{dist}}{z_dist} being \eqn{y_{dist}}{y_dist} if only \eqn{x_{dist}}{x_dist} is provided,
#' and reciprocally.
#'
#' @section Voronoi:
#' The design of the plot is based on a quincunx planting pattern. The whole stand can be thought
#' as a matrix with each cell being a Voronoi sub-plot of two trees. The `cols` argument represent
#' the number of cells repeated in x, and the `rows` argument in y.
#'
#' @section Torricity in ARCHIMED:
#' If the user wants to use the ARCHIMED model for further computations, and if the torricty option
#' is activated, only the Voronoï sample of the design is needed (*i.e.* the smaller possible
#' representation of the system) because the design will be virtually duplicated to infinity.
#'
#' @note  If \eqn{y_{dist}}{y_dist} and \eqn{x_{dist}}{x_dist} are given for the distance
#' between trees and not intra/inter-row distance, a transformation has to be done beforehand (see example)
#'
#' @return A list of two:
#' * the data.frame with all information to make an Open Plant Scene file.
#' * the ggplot object associated to the design
#'
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' # design a plot with a distance of 9.2 m between each palm trees:
#'   design_plot(rows=2, cols= 2, y_dist = 9.2)
#'
design_plot= function(rows=1, cols= 1, x_dist= NULL, y_dist= NULL, x0= 0){

  if(is.null(y_dist)){
    if(is.null(x_dist)){stop('At least x_dist or y_dist are required')}
    y_dist= sqrt(x_dist*x_dist-((x_dist/2)**2))
  }else if(is.null(x_dist)){
    x_dist= sqrt(y_dist*y_dist-((y_dist/2)**2))
  }

  # Voronoi of the quincunx design:
  voronoi_plot= data.frame(x= c(x_dist/4,x_dist/4+x_dist/2),
                           y= c(y_dist/2,y_dist/2+y_dist),
                           xmin= c(0,0), xmax= rep(x_dist,2),
                           ymin= c(0,0), ymax= rep(y_dist*2,2))

  # Matrix of the design (each cell is a Voronoi):
  mat_plot= expand.grid(Row= 1:rows, Col= 1:cols)

  # Full design:
  design=
    mapply(function(Row,Col){
      voronoi_plot%>%
        dplyr::select(.data$x,.data$y,.data$xmax,.data$ymax,
                      .data$xmin,.data$ymin)%>%
        dplyr::mutate(xmin= .data$xmax*(Col-1), ymin= .data$ymax*(Row-1),
                      x= .data$x+.data$xmin, y= .data$y+.data$ymin,
                      xmax= .data$xmax*Col, ymax= .data$ymax*Row,
                      Col= Col, Row= Row)
    }, Row= mat_plot$Row, Col= mat_plot$Col)%>%t()%>%
    dplyr::as_tibble()%>%
    tidyr::unnest()%>%
    dplyr::mutate(xmax= max(.data$xmax), ymax= max(.data$ymax),
                  xmin= min(.data$xmin), ymin= min(.data$ymin))


  result=
    design%>%
    dplyr::mutate(z= 0.0, scale= 1.0,
                  inclinationAzimut= 0.0, inclinationAngle= 0.0,
                  stemTwist= 0.0)%>%
    dplyr::group_by(.data$Col)%>%
    dplyr::mutate(Border_x= ifelse(.data$x==min(.data$x)|
                                     .data$x==max(.data$x),"out","in"))%>%
    dplyr::group_by(.data$Row)%>%
    dplyr::mutate(Border_y= ifelse(.data$y==min(.data$y)|
                                     .data$y==max(.data$y),"out","in"))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(Border= ifelse(.data$Border_x=="out"|.data$Border_y=="out","out","in"))

  # plot:
  plot_bounds=
    result%>%
    ggplot2::ggplot(ggplot2::aes(x= .data$x, y= .data$y, color= .data$Border))+
    ggplot2::geom_point()+
    ggplot2::ylim(low= min(result$ymin), high= max(result$ymax))+
    ggplot2::xlim(low= min(result$xmin), high= max(result$xmax))

  list(design= result, plot= plot_bounds)
}





#' Format ops information
#'
#' @description Make the planting design to the Open Plant Scene format
#' to prepare for OPS writing
#'
#' @param design  The planting design, generally computed using [design_plot()]
#' @param Progeny The progeny name
#' @param map     The tree age in month after planting
#' @param id      The scene ID
#' @param bounds  Are the plot bounds written in the OPS file ? Used for backward
#' compatibility.
#' @param pavement Pavement file path. Only used if the pavement file (`.gwa`) has to
#'  be linked in the OPS. Used for backward compatibility.
#' @param average  Boolean. Use average tree instead of the sampled ones.
#'
#' @return A pre-formatted OPS
#' @export
#'
format_ops=function(design, Progeny, map, id= 1, bounds= TRUE,
                    pavement= NULL, average= FALSE){

  nbTree= nrow(design)

  # generate a line of config file for each tree
  opf_table=
    design%>%
    dplyr::transmute(
      plantId= 1:dplyr::n(),
      out=
        paste(id,.data$plantId,
              if(average){
                paste0('opf/',Progeny,'_Average_MAP_',map,'.opf')
              }else{
                paste0('opf/',Progeny,'_Tree_',.data$plantId,'_MAP_',map,'.opf')
              },
              .data$x,.data$y,.data$z,.data$scale,.data$inclinationAzimut ,
              .data$inclinationAngle, .data$stemTwist,sep='\t'))%>%
    dplyr::select(-.data$plantId)%>%as.matrix

  if(bounds){
    plot_box=
      design%>%
      dplyr::mutate(zmin= 0.0, zmax= 0.0)%>%
      dplyr::select(.data$xmin, .data$ymin, .data$zmin,
                    .data$xmax, .data$ymax, .data$zmax)%>%
      dplyr::summarise_all(function(x){round(unique(x),3)})%>%
      dplyr::transmute("# T xOrigin yOrigin zOrigin xSize ySize flat"=
                         glue::glue("T {xmin} {ymin} {zmin} {xmax} {ymax} {zmax} flat"))%>%
      as.matrix()
    plot_box= rbind(colnames(plot_box),plot_box)

    plot_groups= c("# Archimed metadata format :",
                   '# start comment with "[Archimed]" stuck to # sign',
                   "# field 1 : functional group name",
                   "#[Archimed] two")

  }else{
    plot_box= NULL
    plot_groups= NULL
  }

  if(!is.null(pavement)){
    pav=
      paste(1,21,pavement,unique(design$xmax/2),unique(design$ymax/2),0,1,0,0,0,sep='	')
  }else{
    pav= NULL
  }

  c(plot_box,
    paste('# Part 1: one line per plant in the scene'),
    paste(plot_groups),
    paste('#sceneId plantId plantFileName x y z scale inclinationAzimut inclinationAngle stemTwist'),
    paste(opf_table),
    pav,
    paste('# [Optional] Part 2, chaining: only if scenario or project, one line per sceneId in part1'),
    paste('#motherId sceneId date'),
    paste(-1,id,1, sep='\t')
  )
}


#' Write OPS
#'
#' @param data A pre-formated OPS file, generally the output from [format_ops()]
#' @param path File path and name
#' @param overwrite Boolean. Should pre-existing OPS files overwriten ?
#'
#' @return Writes an OPS to disk, and returns `TRUE` if successfull or `FALSE` otherwise.
#' @export
#'
write_ops= function(data, path, overwrite= T){

  file_time= NULL
  path= normalizePath(path, winslash= "/", mustWork = FALSE)


  if(!dir.exists(file.path(dirname(path)))){
    dir.create(file.path(dirname(path)))
  }else{
    if(file.exists(path)){
      if(!overwrite){
        warning("OPS file ", basename(path), " already exists",
                " and overwrite is set to FALSE. Please set overwrite= TRUE or change the",
                " file name or directory")
        return(basename(path))
      }else{
        file_time= file.mtime(path)
      }
    }
  }

  write(data, file= path)

  # Return TRUE if written and FALSE if not, or not replaced:
  is_written= file.exists(path)
  if(is.null(file_time)){
    return(is_written)
  }else{
    if(file_time<file.mtime(path)){
      return(is_written)
    }else{
      return(FALSE)
    }
  }
}


#' Format and write all OPS
#'
#' @description Format and write OPS from a set of design experiment(s),
#' map and progenies, by applying [format_ops()] and [write_ops()] sequentially
#'
#' @param Progeny The progenies
#' @param design  The design of experiment, generally computed using [design_plot()]
#' @param map     The age of the plantation in month after planting
#' @param path    The path of the target folder
#' @param overwrite Boolean. Should pre-existing OPS files overwriten ?
#' @param average Boolean. Use average tree instead of the sampled ones.
#' @param ...     Further arguments to pass to [format_ops()]
#'
#' @details The function uses [base::mapply()] to apply both [format_ops()] and
#'  [write_ops()] to any number of progeny, design, or map. So if these arguments
#'  are provided with equal length, they will be applied in parallel (*i.e.* in a
#'  multivariate mode) as [base::mapply()] do. Each can also have length one. in this
#'  case it will be re-cycled for each combination.
#'
#' @return An OPS for each progeny, each design and each map.
#' @export
#'
make_ops_all= function(Progeny, design, map, path, overwrite,
                       average= FALSE,...){
  out=
    mapply(function(x,y,z){
      format_ops(design = z, Progeny = x, map = y,
                 average = average,...)%>%
        write_ops(file.path(path,paste0(x,'_MAP_',y,'.ops')),
                  overwrite = overwrite)
    },
    Progeny, map, design)

  if(all(unlist(out))){
    message("All OPS files were successfully written to disk")
    return(TRUE)
  }else{
    stop("Error during OPS file writting for file\n", out[out==FALSE])
  }
}



