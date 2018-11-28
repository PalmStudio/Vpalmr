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
#' @param seed   The seeds for random parameter generation, see details.
#' @param init   Initialisation values (see [extract_trees()])
#' @param average Boolean, is mean of all average trees to be used to replace missing values ?
#'
#' @details `seed` should be a named list of length equal to the number of progenies, and each
#' object with n seeds each. The names of the list objects should match the ones of the
#' progenies. Alternatively, `seed` can be a vector of seeds of length n that will be recycled
#' for each progeny. If `NULL`, a random seed is generated for each tree of each progeny
#' using [base::sample()].
#'
#' @return A nested list, containing a list of VPalm parameters for each tree sampled
#' for each Progeny, and the average tree of each progeny, and one with the average (non-NA) values for all
#' progenies. Each list in the list has as many n lists that represents each sampled tree
#' for each progeny.
#'
#' @importFrom rlang .data
#'
#' @export
extract_progeny= function(data, model, n, leaves= 45, seed= NULL,
                          init= init_list(),average= T){

  # Testing if any model has missing Progenies:
  Progenies=
    lapply(model, function(x){
      unique(x$Progeny)
    })
  Prog= unique(unlist(Progenies))


  # Control/Making of seeds  ------------------------------------------------

  if(is.null(seed)){
    warning("No seeds provided, generating random seeds")
    seed= sample.int(1000,n*length(Prog))%>%
      split(., rep(Prog, each= n))
  }
  if(length(seed)==n & !is.list(seed)){
    # Make a list with repeated seed value for each progeny:
    seed= stats::setNames(rep(list(seed),length(Prog)),Prog)
    warning("Number of seeds match n, recycle seeds through progenies")
  }
  if(length(unlist(seed))!=(n*length(Prog))){
    stop("Number of seeds should match n or total number of trees (n * number of progenies)")
  }
  if(!all(Prog%in%names(seed))){
    stop("Missing seed for progeny: ",paste(Prog[!(Prog%in%names(seed))], ' '))
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

  if(length(missing_string)>0){
    missing_string= paste0(names(missing_string),", progeny: ",missing_string)
    if(average){
      warning("Missing Progeny data to model:",paste("\n* ",missing_string),
              "\nUsing all progenies average data for these missing values")
    }else{
      warning("Missing Progeny data to model:",paste("\n* ",missing_string),
              "\nset 'average' argument to TRUE to use average data for",
              " these missing values")
    }
  }

  # Make the VPalm list for each Progeny, then take average data for
  # parameters that are missing for the progeny (should return NA first):
  Prog_list= list()
  for(i in Prog){
    mod_pro=
      lapply(model, function(x){
        try({x= x%>%dplyr::filter(.data$Progeny==i)}, silent = T)
        x
      })
    data_pro=
      lapply(data, function(x){
        try({x= x%>%dplyr::filter(.data$Progeny==i)}, silent = T)
        x
      })
    Prog_list[[i]]=
      suppressWarnings(
        extract_trees(data= data_pro, model= mod_pro, n = n,
                      leaves= leaves, seed = seed[[i]],
                      init = init)
      )
  }

  # Extract average tree from each progeny in Prog_list
  Averages= lapply(Prog_list, function(x)x$Average)

  # Identifying all parameter names:
  keys= unique(unlist(lapply(Prog_list[[1]], names)))
  # Computing the averages from the average Tree of each progeny:
  Average=
    stats::setNames(
      do.call(mapply, c(FUN= pmean, lapply(Averages, `[`, keys))),
      keys)
  # Order the average by the key just to be sure:
  Average_list= Average[keys]
  # nbLeafEmitted should remain an integer:
  Average_list$nbLeafEmitted= round(Average_list$nbLeafEmitted)

  if(average){
    # Replacing missing values (if any) for each tree by the average value:
    Prog_list=
      lapply(Prog_list, function(y){
        lapply(y,
               function(x){
                 # making sure that x is rightly ordered:
                 x= x[keys]
                 # Finding the x with missing values:
                 cond= lengths(x)<1 | is.na(x)
                 # replacing those with average value:
                 x[cond]= Average_list[cond]
                 x
               })
      })
  }

  return(c(Prog_list, Average= list(Average_list)))
}
