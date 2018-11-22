#' VPalm parameter list for all progenies
#'
#' @description Make a VPalm parameter list for each progeny provided in the
#' data from the [models][mod_stem_diameter()] outputs using [VPalm_list()] and replace
#' missing value by the average of all progenies.
#'
#' @param data  A list of all data (generally from [import_data()])
#' @param model A list of the models fitted on the data (generally from [mod_all()])
#' @param nb_leaves The number of leaves needed
#' @param seed      The seed for random parameter generation
#'
#' @return A list of two:
#' * Progeny, a list of VPalm parameters for each progeny
#' * Average, A list of VPalm parameters with average parameter values
#'
#' @export
make_list= function(data, model, nb_leaves= 45, seed= sample.int(1000,1)){
  # Testing if any model has missing Progenies:
  Progenies=
    lapply(model, function(x){
      unique(x$Progeny)
    })
  Prog= unique(unlist(Progenies))
  nprog=
    lapply(Progenies, function(x){
      setdiff(Prog,x)
    })

  mod_missing= nprog[which(lengths(nprog)>0)]
  mod_missing= mod_missing[-grep("Rep",names(mod_missing))]
  missing_string=
    sapply(mod_missing, function(x){
      paste(x,collapse = ", ")
    })
  missing_string= paste(names(missing_string),":",missing_string)

  warning("Missing Progeny data to model:\n","* ",missing_string,
          "\nUsing average data for all progenies")

  # Make the VPalm list for each Progeny, then take average data for
  # parameters that are missing for the progeny (should return NA first):
  Prog_list= list()
  for(i in Prog){
    mod_pro=
      lapply(model, function(x){
        try({x= x%>%filter(Progeny==i)}, silent = T)
        x
      })
    data_pro=
      lapply(Inputs, function(x){
        try({x= x%>%filter(Progeny==i)}, silent = T)
        x
      })
    Prog_list[[i]]=
      suppressWarnings(
        Vpalmr::VPalm_list(data= data_pro, model= mod_pro,
                           nb_leaves= nb_leaves, seed = seed)
      )
  }

  # Computing the average of all progenies:
  keys= unique(unlist(lapply(Prog_list, names)))
  # grouped_list=
  #   setNames(do.call(mapply, c(FUN=c, lapply(Prog_list, `[`, keys))), keys)
  pmean= function(...){
    x= cbind(...)
    rowMeans(x, na.rm = T)
  }
  Average_list= setNames(do.call(mapply, c(FUN=pmean, lapply(Prog_list, `[`, keys))), keys)

  # Replacing missing values (if any) for each Progeny by the average value:

  Prog_list=
    lapply(Prog_list,
           function(z){
             Map(function(x,y){
               if(length(x)<1||any(is.na(x))){x=y}; x
             },z,Average_list)
           }
    )

  return(list(Progeny= Prog_list, Average= Average_list))
}
