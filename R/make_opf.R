#' Make OPF
#'
#' @description Use VPalm from AMAPStudio to make an Open Plant Format file using
#' a parameter file generally created by [write_progeny()] or [write_tree()].
#'
#' @param parameter The VPalm parameter file path and name (see [write_progeny()]
#' or [write_tree()])
#' @param opf The target OPF file path and name
#' @param AMAPStudio The root path to AMAPStudio
#' @param overwrite Boolean. Should pre-existing OPF files overwriten ?
#'
#' @return Writes an OPF file
#' @export
#'
make_opf= function(parameter,opf,AMAPStudio,overwrite=T){

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

  if(file.exists(opf)&!overwrite){
    stop("OPF file ", basename(opf)," already exist in ", dirname(opf),
         "\nPlease set overwrite= TRUE, or change the destination or file name.")
  }

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
#' @param overwrite Boolean. Should pre-existing OPF files overwriten ?
#' @param parallel  Boolean. Is the OPF making to be distributed on available machine cores?
#' @param NbCores   The number of cores to use for parallel making. If `NULL` (the default)
#' uses all cores minus one.
#'
#' @details The parameter folder should only contain parameter files. Subfolders
#' are tolerated though. The function uses [parallel::detectCores()] to find how many
#' cores are available on the machine. This function has known issues, see help for more
#' details.
#'
#' @return Creates one OPF file from each VPalm parameter file
#' @export
#'
make_opf_all= function(parameter,opf,AMAPStudio,overwrite=T,parallel=T,NbCores=NULL){
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

  existing_opfs= sapply(opf_path, file.exists)
  existing_opfs= existing_opfs[existing_opfs]
  if(length(existing_opfs)>0 & !overwrite){
    stop("One or more OPF file already exist:\n",
         paste(basename(names(existing_opfs)), collapse= ", "),
         "\nPlease set overwrite= TRUE, or change the destination folder.")
  }

  AMAPStudio=
    tryCatch(normalizePath(AMAPStudio, winslash= "/", mustWork = T),
             error=function(cond) {
               message(paste("AMAPStudio not found at", AMAPStudio))
               return(FALSE)
             })

  if(parallel){
    if(is.null(NbCores)){
      NbCores= parallel::detectCores()-1
    }
    cl= parallel::makeCluster(min(NbCores,length(param_files)))
    parallel::makeCluster(cl)
    parallel::clusterExport(cl=cl,
                            varlist=c("AMAPStudio","make_opf","overwrite"),
                            envir=environment())

    parallel::clusterMap(
      cl = cl,
      fun = function(x,y){
        make_opf(parameter = x, opf = y, AMAPStudio = AMAPStudio,
                 overwrite = overwrite)
      }, x= param_files, y= opf_path)

    parallel::stopCluster(cl)

  }else{
    mapply(function(x,y){
      make_opf(parameter = x, opf = y, AMAPStudio = AMAPStudio,
               overwrite = overwrite)
    },x= param_files, y= opf_path)
  }


}
