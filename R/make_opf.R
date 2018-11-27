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
  parameter= file.path(this_wd,parameter)
  opf= file.path(this_wd,opf)
  setwd(AMAPStudio)
  exportFile=
    paste('java -cp bin;ext/* jeeb.workspace.palms.elaeisRaphael.ElaeisArchiTree',
          parameter,opf,sep=' ')
  # add trycatch here
  out= system(exportFile, intern= TRUE)
  setwd(this_wd)
  out
}
