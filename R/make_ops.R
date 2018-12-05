#' Design the plot
#'
#' @description Help designing a quincunx planting pattern of palm stand.
#'
#' @param ntrees The numbers of trees required in the scene
#' @param x_dist The inter-row distance (m). See details.
#' @param y_dist The intra-row distance (m). See details.
#' @param x0     The minimum X coordinates
#'
#' @details If only one of \eqn{x_{dist}}{x_dist} or \eqn{y_{dist}}{y_dist} is given, the function
#'  will compute the other distance using the following equation:
#'  \deqn{\sqrt{z_{dist}^2-\left(\frac{z_{dist}}{2}\right)^2}}{sqrt(x_dist*x_dist-((x_dist/2)**2))}
#' with \eqn{z_{dist}}{z_dist} being \eqn{y_{dist}}{y_dist} if only \eqn{x_{dist}}{x_dist} is provided,
#' and reciprocally.
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
#'   design_plot(ntrees = 5, y_dist = 9.2)
#'
design_plot= function(ntrees= 2, x_dist= NULL, y_dist= NULL, x0= 0){

  if(is.null(y_dist)){
    if(is.null(x_dist)){stop('At least x_dist or y_dist are required')}
    y_dist= sqrt(x_dist*x_dist-((x_dist/2)**2))
  }else if(is.null(x_dist)){
    x_dist= sqrt(y_dist*y_dist-((y_dist/2)**2))
  }

  if(ntrees<9){
    ntrees_tot= 9
  }else{
    ntrees_tot= ntrees
  }

  plan=
    data.frame(x= rep(NA,ntrees_tot), y= rep(NA,ntrees_tot),
               Row= rep(seq_len(floor(sqrt(ntrees_tot))),length.out= ntrees_tot))%>%
    dplyr::mutate(Col= rep(seq_len(ceiling(sqrt(ntrees_tot))),
                           each= length(unique(.data$Row)))[1:ntrees_tot])%>%
    dplyr::mutate(odd= .data$Row%%2,
                  x= (.data$Row-1/2)*x_dist+x0,
                  y= ifelse(.data$odd==1,(.data$Col-1/2)*y_dist, .data$Col*y_dist))

  # Particular cases:
  if(ntrees!=2){
    plan=
      plan%>%
      dplyr::arrange(.data$y,.data$x)
  }
  if(ntrees==4){
    plan=
      plan%>%
      dplyr::filter(.data$x!=max(.data$x))
  }

  if(ntrees<9){
    plan= plan[1:ntrees,]
  }

  # # Test for eventual issues:
  # if (round(plan[1,]$y-plan[nCol,]$y,4)!= round(y_dist,4)){
  #   warning('Toricity issue for intra-row spacing')
  # }
  # if (round(plan[1,]$x-plan[nRow*nCol,]$x,4)!=round(x_dist,4)){
  #   warning('Toricity issue for inter-row spacing')
  # }

  result=
    plan%>%
    dplyr::mutate(z= 0.0, scale= 1.0,
                  inclinationAzimut= 0.0, inclinationAngle= 0.0,
                  stemTwist= 0.0,
                  ymin= min(.data$y)-y_dist/2,
                  ymax= max(.data$y)+y_dist/2,
                  xmin= min(.data$x)-x_dist/2,
                  xmax= max(.data$x)+x_dist/2)%>%
    dplyr::group_by(.data$Col)%>%
    dplyr::mutate(Border_x= ifelse(.data$x==min(.data$x)|
                                     .data$x==max(.data$x),"out","in"))%>%
    dplyr::group_by(.data$Row)%>%
    dplyr::mutate(Border_y= ifelse(.data$y==min(.data$y)|
                                     .data$y==max(.data$y),"out","in"))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(Border= ifelse(Border_x=="out"|Border_y=="out","out","in"))


  # plot:
  plot_bounds=
    result%>%
    ggplot2::ggplot(ggplot2::aes(x= x, y= y, color= Border))+
    ggplot2::geom_point()+
    ggplot2::ylim(low= min(result$y), high= max(result$y))

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
#' compatibility
#' @param pavement Pavement file path. Only used if the pavement file (`.gwa`) has to
#'  be linked in the OPS. Used for backward compatibility.
#'
#' @return A pre-formatted OPS
#' @export
#'
format_ops=function(design, Progeny, map, id= 1, bounds= FALSE, pavement= NULL){

  nbTree= nrow(design)

  # generate a line of config file for each tree
  opf_table=
    design%>%
    dplyr::transmute(
      plantId= 1:dplyr::n(),
      out=
        paste(id,plantId,paste0('opf/',Progeny,'_Tree_',plantId,'_MAP_',map,'.opf'),
              .data$x,.data$y,.data$z,.data$scale,.data$inclinationAzimut ,
              .data$inclinationAngle, .data$stemTwist,sep='\t')
    )%>%dplyr::select(-.data$plantId)%>%as.matrix

  if(bounds){
    plot_box=
      design%>%
      dplyr::mutate(zmin= 0, zmax= 0)%>%
      dplyr::select(.data$xmin, .data$ymin, .data$zmin,
                    .data$xmax, .data$ymax, .data$zmax)%>%
      dplyr::summarise_all(function(x){round(unique(x),3)})%>%
      dplyr::transmute("# T xOrigin yOrigin zOrigin xSize ySize flat"=
                         glue::glue("T {xmin} {ymin} {zmin} {xmax} {ymax} {zmax} flat"))%>%
      as.matrix()
  }else{
    plot_box= NULL
  }

  if(!is.null(pavement)){
    pav=
      paste(1,21,pavement,unique(design$xmax/2),unique(design$ymax/2),0,1,0,0,0,sep='	')
  }else{
    pav= NULL
  }

  c(plot_box,
    paste('# Part 1: one line per plant in the scene'),
    paste('#sceneId plantId plantFileName x y z scale inclinationAzimut inclinationAngle stemTwist'),
    paste(opf_table),
    pav,
    paste('# [Optional] Part 2, chaining: only if scenario or project, one line per sceneId in part1'),
    paste('#motherId sceneId date'),
    paste(-1,id,1)
  )
}


#' Write OPS
#'
#' @param data A pre-formated OPS file, generally the output from [format_ops()]
#' @param path File path and name
#' @param overwrite Boolean. Should pre-existing OPS files overwriten ?
#'
#' @return Writes an OPS to disk
#' @export
#'
write_ops= function(data, path, overwrite= T){
  if(!dir.exists(file.path(dirname(path)))){
    dir.create(file.path(dirname(path)))
  }else{
    if(file.exists(path) & !overwrite){
      warning("OPS file ", basename(path), " already exists",
              " and overwrite is set to FALSE. Please set overwrite= TRUE or change the",
              " file name or directory")
      return(basename(file_name))
    }
  }
  write(data, file= path)
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
#'
#' @details The function uses [base::mapply()] to apply both [format_ops()] and
#'  [write_ops()] to any number of progeny, design, or map. So if these arguments
#'  are provided with equal length, they will be applied in parallel (*i.e.* in a
#'  multivariate mode) as [base::mapply()] do. Each can also have length one, so it
#'  will be re-cycled for each combination.
#'
#' @return An OPS for each progeny, each design and each map.
#' @export
#'
make_ops_all= function(Progeny, design, map, path){
  mapply(function(x,y,z){
    format_ops(z,x,y)%>%
      write_ops(file.path(path,paste0(x,'_',y,'MAP.ops')))
  },
  Progeny, map, design)
}



