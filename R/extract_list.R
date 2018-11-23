#' VPalm parameter list for all progenies
#'
#' @description Make a VPalm parameter list for n trees sampled in the parameters
#' distributions of each progeny, and also for an average tree from the n trees distribution.
#'  provided in the data from the [models][mod_stem_diameter()]
#' outputs using [extract_params()] and replace missing value by the average of all progenies.
#'
#' @param data   A list of all data (generally from [import_data()])
#' @param model  A list of the models fitted on the data (generally from [mod_all()])
#' @param n      The number of samplings for each progeny, or in other words, the
#'  number of trees to be generated for each progeny
#' @param leaves The number of leaves needed
#' @param seed   The seeds for random parameter generation. Should be a named list of length equal
#' to the number of progenies (names = names of the progenies), with n seeds each. If it is
#' a length n vector only, it will be recycled through each progeny. If `NULL`, take a random seed.
#' @param init   Initialisation values (see [extract_trees()])
#'
#' @return A list of several VPalm parameters lists, one for each tree sampled per progeny
#' and the average tree of the progeny, and one with the average (non-NA) values for all
#' progenies. Each list in the list has as many n lists that represents each sampled tree
#' for each progeny.
#'
#' @export
extract_progenies= function(data, model, n, leaves= 45, seed= NULL,
                            init= init_list()){

  # Testing if any model has missing Progenies:
  Progenies=
    lapply(model, function(x){
      unique(x$Progeny)
    })
  Prog= unique(unlist(Progenies))

  # Testing the number of seeds:
  if(length(seed)==n & !is.list(seed)){
    # Make a list with repeated seed value for each progeny:
    seed= setNames(rep(list(seed),length(Prog)),Prog)
    warning("Number of seeds match n, recycle seeds through progenies")
  }
  if(length(unlist(seed))!=(n*length(Prog))&is.null(seed)){
    stop("Number of seeds should match n or total number of trees (n * number of progenies)")
  }
  if(all(Prog%in%names(seed))){
    stop("Missing seed for progeny",Prog[!(Prog%in%names(seed))])
  }

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
        try({x= x%>%dplyr::filter(Progeny==i)}, silent = T)
        x
      })
    data_pro=
      lapply(Inputs, function(x){
        try({x= x%>%dplyr::filter(Progeny==i)}, silent = T)
        x
      })
    Prog_list[[i]]=
      suppressWarnings(
        Vpalmr::extract_trees(data= data_pro, model= mod_pro,
                              leaves= leaves, seed = seed[[i]],
                              init = init)
      )
  }


  return(c(Prog_list, Average= list(Average_list)))
}
