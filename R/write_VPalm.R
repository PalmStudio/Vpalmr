
#' Writes all progenies VPalm inputs
#'
#' @description Applies [write_tree()] successively to a [format_progeny()] output to
#' write all trees from all progenies to a VPalm input file.
#'
#' @param data A pre-formated [format_progeny()] output list of matrices
#' @param path The path of the directory to write in
#'
#' @export
#'
write_progeny= function(data,path){
  prog= names(data)
  name= mapply(function(x,y)paste(y,names(x),sep='_'),x=data, y=prog)

  if(!is.null(data$Average)){
    write_tree(data$Average, path = path, name = "All_Progenies_Average")
    data= data[-grep('Average',names(data))]
    name= name[-grep("Average",names(name))]
  }

  mapply(FUN = function(x,y){
    mapply(function(z,a){
      write_tree(data = z, path = path, name= a)
    },z= x, a= y)
  }, x= data, y= name)

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
#'
#' @details If the `age` argument is not provided, the function will try to extract it from
#' the matrix (from the `Modelled Months After Planting` row name).
#'
#' @export
#'
write_tree= function(data, path= getwd(), name= "vpalm_input", age=NULL){
  if(!dir.exists(file.path(path))){
    dir.create(file.path(path))
  }

  MAP= grep('Modelled Months After Planting',rownames(data))
  if(is.null(age)&length(MAP)>0){
    age= paste0('_',data[MAP,])
    data= data[-MAP,]
  }else{
    age= '_unknown'
  }
  writeLines(paste(names(data),data),file.path(path,paste0(name,age,'.txt')))
  message("VPalm input file written in ",
          file.path(getwd(),path,paste0(name,age,'.txt')), "\n")
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
