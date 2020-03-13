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
#' @param verbose Should the VPalm writting informations printed to the console ?
#' @param java    Java path (optionnal, see details).
#'
#' @details The `java` argument can be a path to the java
#' executable if the user needs a particular version (for example if the default Java used is the Open JDK,
#' because ARCHIMED is only compatible with the Oracle version).
#'
#' @return Writes an OPF file, and return `TRUE` if the file was successfully written.
#' @export
#'
make_opf= function(parameter,opf,AMAPStudio,overwrite=T,verbose=F, java=NULL){

  # Normalize all paths:
  AMAPStudio= normalizePath(AMAPStudio, winslash = "/", mustWork = TRUE)
  opf= normalizePath(opf, winslash = "/", mustWork = FALSE)

  if(!dir.exists(file.path(dirname(opf)))){
    dir.create(file.path(dirname(opf)), recursive = T)
  }

  # If the user input a jar file, only use the directory name:
  if(grepl('.jar',AMAPStudio)){
    AMAPStudio= dirname(AMAPStudio)
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

  if(file.exists(opf)){
    if(!overwrite){
      stop("OPF file ", basename(opf)," already exist in ", dirname(opf),
           "\nPlease set overwrite= TRUE, or change the destination or file name.")
    }else{
      file_time= file.mtime(opf)
    }
  }else{
    file_time= NULL
  }

  # Test if there are some white spaces in the config path, and if so double quote it:
  if(grepl(" ",parameter)){
    parameter= paste0("\"",parameter,"\"")
  }
  # Test if there are some white spaces in the config path, and if so double quote it:
  if(grepl(" ",opf)){
    opf_call= paste0("\"",opf,"\"")
  }else{
    opf_call= opf
  }

  on.exit(setwd(this_wd))
  setwd(AMAPStudio)
  if(is.null(java)){
    java= "java"
  }

  exportFile=
    paste(java,'-jar', "vpalm.jar",parameter,opf_call,sep=' ')
    # paste('java -cp bin;ext/* jeeb.workspace.palms.elaeisRaphael.ElaeisArchiTree',
    #       parameter,opf,sep=' ')
  out=
    tryCatch(system(exportFile, intern= TRUE),
             error=function(cond) {
               message("VPalm encountered an issue")
               return(FALSE)
             },
             finally= setwd(this_wd))
  if(verbose){
    message(out)
  }

  # Return TRUE if written and FALSE if not, or not replaced:
  is_written= file.exists(opf)
  if(is.null(file_time)){
    return(is_written)
  }else{
    if(file_time<file.mtime(opf)){
      return(is_written)
    }else{
      return(FALSE)
    }
  }
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
#' @param java    Java path (optionnal, see details)
#'
#' @details The parameter folder should only contain parameter files. Subfolders
#' are tolerated though. The function uses [parallel::detectCores()] to find how many
#' cores are available on the machine. This function has known issues, see help for more
#' details. The `java` argument can be a path to the java
#' executable if the user needs a particular version (for example if the default Java used is the Open JDK,
#' because ARCHIMED is only compatible with the Oracle version).
#'
#' @return Creates one OPF file from each VPalm parameter file, and returns `TRUE` if all OPFs
#' were successfully written.
#' @export
#'
make_opf_all= function(parameter,opf,AMAPStudio,overwrite=T,parallel=T,NbCores=NULL,java=NULL){

  # normalize all paths:
  AMAPStudio= normalizePath(AMAPStudio, winslash = "/", mustWork = TRUE)
  opf= normalizePath(opf, winslash = "/", mustWork = FALSE)

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
    normalizePath(., winslash= "/", mustWork = FALSE)

  existing_opfs= sapply(opf_path, file.exists)
  existing_opfs= existing_opfs[existing_opfs]
  if(length(existing_opfs)>0 & !overwrite){
    stop("One or more OPF file already exist:\n",
         paste(basename(names(existing_opfs)), collapse= ", "),
         "\nPlease set overwrite= TRUE, or change the destination folder.")
  }

  AMAPStudio=
    tryCatch(AMAPStudio,
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
    out=
      parallel::clusterMap(
        cl = cl,
        fun = function(x,y){
          out_tmp= make_opf(parameter = x, opf = y, AMAPStudio = AMAPStudio,
                            overwrite = overwrite, java = java)
          out_tmp
        }, x= param_files, y= opf_path)

    parallel::stopCluster(cl)

  }else{
    out=
      mapply(function(x,y){
        make_opf(parameter = x, opf = y, AMAPStudio = AMAPStudio,
                 overwrite = overwrite)
      },x= param_files, y= opf_path)
  }
  if(all(unlist(out))){
    message("All OPF files were successfully written to disk")
    return(TRUE)
  }else{
    stop("Error during OPF file writting for file\n", opf_path[opf_path==FALSE])
  }

}
