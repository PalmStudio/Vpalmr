#' Make OPF
#'
#' @description Use VPalm from AMAPStudio to make an Open Plant Format file using
#' a parameter file generally created by [write_progeny()] or [write_tree()].
#'
#' @param parameter The VPalm parameter file path and name (see [write_progeny()]
#' or [write_tree()])
#' @param opf The target OPF file path and name
#' @param AMAPStudio The root path to AMAPStudio
#'
#' @return Writes an OPF file
#' @export
#'
make_opf= function(parameter,opf,AMAPStudio){

  if(!dir.exists(file.path(dirname(opf)))){
    dir.create(file.path(dirname(opf)))
  }

  this_wd= getwd()

  parameter=
    tryCatch(
      normalizePath(parameter, winslash= "/", mustWork = T),
      error=function(cond) {
        message(paste("The following input parameter file do not exist:"))
        message(sub('.*?\\"(.*?)\\".*', '\\1', cond))
        return(FALSE)
      }
    )

  opf= normalizePath(opf, winslash= "/", mustWork = F)
  setwd(AMAPStudio)
  exportFile=
    paste('java -cp bin;ext/* jeeb.workspace.palms.elaeisRaphael.ElaeisArchiTree',
          parameter,opf,sep=' ')
  # add trycatch here
  out= system(exportFile, intern= TRUE)
  setwd(this_wd)
  out
}





#' Make a set of opf file
#'
#' @description Make a call to [make_opf()] successively on all VPalm parameter
#' files present on the `parameter` folder to create an Open Plant Format from
#' each.
#'
#' @param parameter The folder of the parameter files
#' @param opf The target folder for resulting OPF files
#' @param AMAPStudio The root path to AMAPStudio where VPalm lives
#'
#' @details The parameter folder should only contain parameter files. Subfolders
#' are tolerated though.
#'
#' @return Creates one OPF file from each VPalm parameter file
#' @export
#'
make_opf_all= function(parameter,opf,AMAPStudio){
  param_files=
    list.files(parameter,full.names = T)%>%
    .[!file.info(.)$isdir]%>%
    normalizePath(., winslash= "/", mustWork = T)

  if(is.null(param_files)|length(param_files)<1){
    stop("No parameter files found in: ",parameter)
  }

  opf_path=
    basename(param_files)%>%gsub("\\.txt","\\.opf",.)%>%
    file.path(opf,.)%>%
    normalizePath(., winslash= "/", mustWork = F)

  AMAPStudio=
    tryCatch(normalizePath(AMAPStudio, winslash= "/", mustWork = T),
             error=function(cond) {
               message(paste("AMAPStudio not found at", AMAPStudio))
               return(FALSE)
             })

  mapply(function(x,y){
    make_opf(parameter = x, opf = y, AMAPStudio = AMAPStudio)
  },x= param_files, y= opf_path)
}