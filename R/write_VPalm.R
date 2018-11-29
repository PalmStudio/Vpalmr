
#' Writes all progenies VPalm inputs
#'
#' @description Applies [write_tree()] successively to a [format_progeny()] output to
#' write all trees from all progenies to a VPalm input file.
#'
#' @param data A pre-formated [format_progeny()] output list of matrices
#' @param path The path of the directory to write in
#' @param wforce  Boolean. If `TRUE` (the default), continue to write files even if one
#' returned an error, and return warnings instead. If `FALSE`, stop at first error.
#' @param overwrite Boolean. Overwrite a pre-existing VPalm input file if set to `TRUE`
#' @param verbose Prints a messages during writing. Especially usefull for debugging.
#'
#' @return A matrix informing about the success of the writing of the files. All values
#' should be `TRUE` if all file were successfully written. If any error was encountered, and if
#' `wforce= TRUE` all file writing the file name if any error occured.
#'
#' @export
#'
write_progeny= function(data, path= getwd(), wforce= FALSE,
                        overwrite = TRUE, verbose= TRUE){
  prog= names(data)
  name= mapply(function(x,y)paste(y,names(x),sep='_'),x=data, y=prog)

  if(!is.null(data$Average)){
    success= write_tree(data$Average, path = path, name = "All_Progenies_Average",
                        overwrite = overwrite, verbose= verbose)
    if(success!=TRUE){
      if(wforce){
        warning("write_tree returned an error during writing for file: ",success)
      }else{
        stop("write_tree returned an error during writing for file: ",success)
      }
    }
    data= data[-grep('Average',names(data))]
    name= name[-grep("Average",names(name))]
  }

  successes=
    mapply(FUN = function(x,y){
      mapply(function(z,a){
        success_tmp= write_tree(data = z, path = path, name= a,
                                overwrite = overwrite, verbose= verbose)
        if(success_tmp!=TRUE){
          if(wforce){
            warning("write_tree returned an error during writing for file: ",success_tmp)
          }else{
            stop("write_tree returned an error during writing for file: ",success_tmp)
          }
        }
        success_tmp
      },z= x, a= y)
    }, x= data, y= name)

  return(cbind(successes,Average= success))
}



#' Write VPalm tree parameter file
#'
#' @description Writes the input parameter file for VPalm from a pre-formated
#' matrix, generally from [format_tree()]
#'
#' @param data A formated matrix input
#' @param path The target directory
#' @param name (Optionnal) The progeny/tree name
#' @param age  (Optionnal) The age at which the progeny parameters were
#'  computed from (see details)
#' @param overwrite Boolean. Overwrite a pre-existing VPalm input file if set to `TRUE`
#' @param verbose Prints a message if the writing is successfull
#'
#' @details If the `age` argument is not provided, the function will try to extract it from
#' the matrix (from the `Modelled Months After Planting` row name).
#'
#' @return `TRUE` if the file was successfully written, `FALSE` if any error
#' was encountered during the process.
#'
#' @export
#'
write_tree= function(data, path= getwd(), name= "vpalm_input", age=NULL,
                     overwrite= TRUE, verbose= TRUE){

  MAP= grep('Modelled Months After Planting',rownames(data))
  if(is.null(age)&length(MAP)>0){
    age= paste0('_',data[MAP,])
    data= data[-MAP,]
  }else{
    age= '_unknown'
  }

  file_name= file.path(path,paste0(name,age,'.txt'))

  if(!is.data.frame(data)){
    warning("data should be a matrix")
    return(basename(file_name))
  }


  if(!dir.exists(file.path(path))){
    dir.create(file.path(path))
  }else{
    if(file.exists(file_name) & !overwrite){
      warning("VPalm parameter file ", basename(file_name), " already exists",
              " and overwrite is set to TRUE. Please set overwrite= FALSE or change the",
              " file name or directory")
      return(basename(file_name))
    }
  }
  tryCatch({
    writeLines(paste(names(data),data),file_name)
    if(verbose){message("VPalm input file written in ",file_name)}
    return(TRUE)
  },
  error=function(cond){
    if(verbose){
      message(paste("Couldn't write VPalm parameters to", basename(file_name)))
      message("Error message from writeLines:")
      message(cond,"\n")
    }
    return(basename(file_name))
  },
  warning=function(cond){
    if(verbose){
      message(paste("warning during VPalm parameters writing to", basename(file_name)))
      message("Warning message from writeLines:")
      message(cond,"\n")
    }
    return(TRUE)
  }
  )
}


#' Write model
#'
#' @description Writes the models outputs to disk.
#' @param data Model outputs, generally from [mod_all()]
#' @param path Path to the folder where to write it
#' @param name Name of the file to write
#'
#' @details The file can be read back using [base::readRDS()]
#' @export
#'
write_models= function(data,path,name= 'models'){
  saveRDS(data, file = file.path(path,paste0(name,".RData")))
}
