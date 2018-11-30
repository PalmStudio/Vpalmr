#' Design the plot
#'
#' @description Help designing the planting pattern of the plot
#'
#' @param nRow Number of tree rows
#' @param nCol Number of tree columns
#' @param x0   The minimum X coordinates
#' @param x_distance The inter-row distance
#' @param y_distance The intra-row distance
#'
#' @details For the moment the function only creates triangular planting pattern,
#' but more options will be added with time
#'
#' @return A data.frame with all information to make an Open Plant Scene file.
#' @export
#'
#' @examples
#' design_plot(nRow = 4, nCol = 5, x0 = 0, y_distance = 9.2,
#'             x_distance = 8)
#'
design_plot= function(nRow,nCol,x0,y_distance,x_distance){
  plan= data.frame(x=rep(NA,nRow*nCol),y=rep(NA,nRow*nCol),
                   Row=rep(c(1:nRow),each=nCol),Col=rep(c(1:nCol)))

  plan$x= (plan$Row-1/2)*x_distance

  for (i in 1:nRow){
    if (i%%2==0){
      plan[plan$Row==i,]$y=
        seq(y_distance-y_distance/4,y_distance*(nCol),y_distance)
    }
    if (i%%2==1){
      plan[plan$Row==i,]$y=
        seq(-y_distance/4,y_distance*(nCol-1),y_distance)+y_distance/2
    }
  }

  # inner palms:
  plan$Border='out'
  plan[plan$Row!=max(plan$Row) &
         plan$Row!=min(plan$Row) &
         plan$Col!=max(plan$Col) &
         plan$Col!=min(plan$Col),]$Border='In'

  # Test for eventual issues:
  if (round(plan[1,]$y-plan[nCol,]$y,4)!= round(y_distance,4)){
    warning('Toricity issue for intra-row spacing')
  }
  if (round(plan[1,]$x-plan[nRow*nCol,]$x,4)!=round(x_distance,4)){
    warning('Toricity issue for inter-row spacing')
  }

  # borders:
  ymin= min(plan$y)-y_distance/4
  ymax= max(plan$y)+y_distance/4
  xmin= min(plan$x)-x_distance/2
  xmax= max(plan$x)+x_distance/2

  result= data.frame(sceneId=1,
                     plantId=paste(rownames(plan)),
                     plantFileName= paste('opf/DA1_Tree',rownames(plan),
                                          '_47MAP.opf',sep=''),
                     x= plan$x, y= plan$y, z= 0.0, scale= 1.0,
                     inclinationAzimut= 0.0, inclinationAngle= 0.0,
                     stemTwist= 0.0, Border= plan$Border)

  plot_bounds=
    result[1:2,]%>%
    dplyr::mutate_all(function(x)x=NA)%>%
    dplyr::mutate(plantId= c('min','max'),
                  x= c(xmin,xmax),y= c(ymin,ymax),
                  Border= 'Border')

  rbind(result,plot_bounds)
}





#' Format ops information
#'
#' @description Make the planting design to the Open Plant Scene format
#' to prepare for OPS writing
#'
#' @param design  The planting design, generally computed using [design_plot()]
#' @param Progeny The progeny name
#' @param map     The tree age in month after planting
#'
#' @return A pre-formatted OPS
#' @export
#'
format_ops=function(design,Progeny,map){

  # map = month after planting (~palm age)

  # number of opf in the ops (if change need to change design file as well!!!!)
  nbTree= nrow(head(design,-2))
  xmax= max(design$x)
  ymax= max(design$y)
  ###generate ligne of config file for each tree
  opf_table=NULL
  for (t in 1:nbTree){
    tableSub=
      paste(1,t,paste('opf/',Progeny,'_Tree_',t,'_MAP_',map,'.opf',sep=''),
            design$x[t],design$y[t],design$z[t],design$scale[t],
            design$inclinationAzimut [t],design$inclinationAngle[t],
            design$stemTwist[t],sep='\t')
    opf_table= rbind(opf_table,tableSub)
  }

  c(
    # paste('# T xOrigin yOrigin zOrigin xSize ySize flat'),
    # paste('T 0 0 0 ',xmax,' ',ymax,' flat',sep=''),
    paste('# Part 1: one line per plant in the scene'),
    paste('#sceneId plantId plantFileName x y z scale inclinationAzimut inclinationAngle stemTwist'),
    paste(opf_table[1,]),
    paste(opf_table[2,],collapse=' '),
    paste(opf_table[3,],collapse=' '),
    paste(opf_table[4,],collapse=' '),
    paste(opf_table[5,],collapse=' '),
    paste(opf_table[6,],collapse=' '),
    paste(opf_table[7,],collapse=' '),
    paste(opf_table[8,],collapse=' '),
    paste(opf_table[9,],collapse=' '),
    paste(opf_table[10,],collapse=' '),
    paste(opf_table[11,],collapse=' '),
    paste(opf_table[12,],collapse=' '),
    paste(opf_table[13,],collapse=' '),
    paste(opf_table[14,],collapse=' '),
    paste(opf_table[15,],collapse=' '),
    paste(opf_table[16,],collapse=' '),
    paste(opf_table[17,],collapse=' '),
    paste(opf_table[18,],collapse=' '),
    paste(opf_table[19,],collapse=' '),
    paste(opf_table[20,],collapse=' '),
    paste(1,21, paste('opf/pavement4x5.gwa',sep=''),xmax/2,ymax/2,0,1,0,0,0,sep='	'),
    paste('# [Optional] Part 2, chaining: only if scenario or project, one line per sceneId in part1'), #####pavement to be added in opfFile!!!
    paste('#motherId sceneId date'),
    paste(-1,1,1,sep='	')
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



