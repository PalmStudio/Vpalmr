#' Write VPalm parameter lists
#'
#' @description Writes all the lists from [format_list()] as VPalm input files.
#'
#' @details This function applies [write_params()] sequentially to all progenies and to
#' the average list.
#'
#' @param data A [format_list()] output
#' @param path The path of the directory to write the files to.
#' @param name (Optionnal) The names of the files. If `NULL`, the names of the progenies
#' will be used
#'
#' @export
#'
#' @examples
#' \dontrun{
#' out= format_list(data = VPalm_out)
#' write_list(data = out)
#' }
#'
write_list= function(data, path= getwd(), name= NULL){

  if(is.null(name)){
    name= names(data)
  }else{
    if(length(name)!=length(data)){
      stop("The name parameter length do not match the number of lists in data")
    }
  }


  filenamePro= file.path(path,name,"MAP.txt")

  write_param(data = data, )

  # Average values from all progenies:
  filenamePro= file.path(path,paste(pro,'_AverageTree_', map,"MAP.txt",sep = ""))
  capture.output(ArchiOutputMean(Progeny=pro, nbLeafEmitted = nbLeafEmitted),file = filenamePro)

  ###individuals
  for (i in (1:20)){
    filename= file.path("Outputs/ParametersFiles",paste(pro,'_','Tree',i,'_', map,"MAP.txt",sep = ""))
    capture.output(ArchiOutput(Progeny=pro, nbLeafEmitted = nbLeafEmitted,seed=i),file = filename)
  }
}


#' Write VPalm parameter
#'
#' @description Writes a parameter list from [format_param()] (or one from
#'  [format_list()]) as a VPalm input files.
#'
#'
#' @param data A [format_param()] output (or one of the [format_list()])
#' @param path The path of the VPalm input file
#'
#' @export
#'
#' @examples
#' \dontrun{
#' out= format_list(data = VPalm_out)
#' write_param(data = out$Average)
#' }
#'
write_param= function(data,path= getwd()){
  capture.output(data,file = path)
}
