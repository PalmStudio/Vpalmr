#' Extract several VPalm parameters
#'
#' @description Extract Vpalm parameters for n trees and for the average tree using
#' [extract_params()] recursively
#'
#' @param data  A list of all data (generally from [import_data()])
#' @param model A list of the models fitted on the data (generally from
#'  [mod_all()])
#' @param n     The number of trees desired (how many samples are made)
#' @param leaves The number of leaves required
#' @param seed  The seed for random parameter generation (must be length n)
#' @param init  The initialisation parameters (see details)
#'
#' @details n controls the number of trees sampled in the parameter value
#' distribution. Parameters sampled include:
#' * The rachis curvature
#' * The rachis relative height
#' * All lme models (see [models][mod_stem_diameter()])
#' * All nlme models (see [models][mod_stem_diameter()])
#'
#' The init parameter should be a named list of initialization
#' parameters as in [init_list()]. The user can use [init_list()] as a
#' template to change the values of the init argument (see examples).
#'
#' Missing values are replaced by values from the average tree
#'
#' @return n+1 list of the parameters used as VPalm input, one for each n tree, and
#' one for the average tree
#'
#' @export
#'
extract_trees= function(data, model, n, leaves= 45, seed= sample.int(1000,n),
                        init= init_list()){
  if(length(seed)!=n){
    stop("Number of seeds do not match number of samples (n)")
  }

  Trees_params=
    lapply(seed, function(x){
      extract_params(data = data, model = model, leaves = leaves, seed = x,
                     init= init, type= 'sample')
    })
  names(Trees_params)= paste0('Tree_',1:n)

  Average_params=
    extract_params(data = data, model = model, leaves = leaves, seed = NULL,
                   init= init, type= 'mean')

  # Identifying the progenies parameter values:
  keys= unique(unlist(lapply(Trees_params, names)))

  # Order the average by this value just to be sure:
  Average_list= Average_params[keys]

  # Replacing missing values (if any) for each tree by the average value:
  Trees_params=
    lapply(Trees_params,
           function(x){
             # making sure that x is rightly ordered:
             x= x[keys]
             # Finding the x with missing values:
             cond= lengths(x)<1 | is.na(x)
             # replacing those with average value:
             x[cond]= Average_list[cond]
             x
           }
    )

  return(c(Trees_params, Average= list(Average_list)))
}



#' Extract VPalm parameters
#'
#' @description Format the model outputs as a list for VPalm inputs
#'
#' @param data  A list of all data (generally from [import_data()])
#' @param model A list of the models fitted on the data (generally from
#'  [mod_all()])
#' @param leaves The number of leaves needed
#' @param seed  The seed for random parameter generation
#' @param init  The initialisation parameters (see details)
#' @param type The type of sampling performed on the parameter distribution (see
#'  [coef_sample()])
#'
#' @details Parameters sampled include:
#' * The rachis curvature
#' * The rachis relative height
#' * All lme models (see [models][mod_stem_diameter()])
#' * All nlme models (see [models][mod_stem_diameter()])
#' The type parameter controls whether a sample or an average tree is needed.
#'
#' The init parameter should be a named list of initialization
#' parameters as in [init_list()]. The user can use [init_list()] as a
#' template to change the values of the init argument (see examples).
#'
#' @return A list of the parameters used as VPalm input
#' @export
#'
#' @examples
#' # Change the initial values of lamina angle :
#' ini= init_list()
#' ini$laminaAngle= 145
#' \dontrun{
#'  # call the function to create a VPalm input:
#'  extract_params(data,model, leaves= 45, seed= sample.int(1000,1),
#'  initialisation= ini)
#' }
extract_params= function(data, model, leaves= 45, seed= sample.int(1000,1),
                         init= init_list(), type= c('sample','mean')){


  # Sampling the parameter distribution: ------------------------------------

  set.seed(seed)

  # nlme:
  model$leafleftLength.nlme%<>%coef_sample(epsilon= 5e-06, type= type)
  model$leafletWidth.nlme%<>%coef_sample(epsilon= 10^-6, type= type)
  model$axialAngle.nlme%<>%coef_sample(epsilon= 10^-2, type= type)

  # lme:
  model$rachisLength.lme%<>%pull_lme(epsilon= 10^-6, type= type)
  model$decliC.lme%<>%pull_lme(epsilon= 0, type= type)
  model$leaflet_length_B.lme%<>%pull_lme(epsilon= 2e-06, type= type)
  model$leaflet_width_B.lme%<>%pull_lme(epsilon= 2e-06, type= type)
  model$petioleWidthC.lme%<>%pull_lme(epsilon= 10^-6, type= type)


  VP_list=
    list(
      seed= seed,

      nbLeafEmitted= mean(data$DataAll$Physio_age),
      # FIXED PARAMETERS PER PROGENY --------------------------------------------
      # Stem phylotaxis
      frondPhyllo_M= mean(data$Phylo$Phylo,na.rm=T),
      frondPhyllo_SD= sd(data$Phylo$Phylo,na.rm=T),

      # stem height
      residStemHeight= model$model.stemHeight$sigma,
      coefStemHeight= model$model.stemHeight$coef,
      # stem diameter
      stemDiamMax= model$StemDiam.nls$finalStemDiam,
      stemDiamSlope= model$StemDiam.nls$slope,
      stemDiamInfl= model$StemDiam.nls$Infl,
      residStemDiam= model$StemDiam.nls$sigma,

      # number of frond
      nbFronds_M= leaves,
      nbFronds_SD= 0,

      # frond twist
      rachisTwistFinalAngle_M= mean(data$Tor$TwistA,na.rm=T),
      rachisTwistFinalAngle_SD= sd(data$Tor$TwistA,na.rm=T),

      # frond deviation
      rachisDevFinalAngle_M= model$Dev$DevA_deg,
      rachisDevFinalAngle_SD= model$Dev$DevA_deg,

      # parameters ratio petiole/rachis final
      petioleRachisRatio_M= model$Pet$RatioPetiole,
      petioleRachisRatio_SD= model$Pet$RatioPetiole_sd,

      # B point relative position on rachis
      pointBrelativePosition_M= model$Bpos$PosB,
      pointBrelativePosition_SD= model$Bpos$PosB_sd,

      # leaflet number
      nbMax= model$nbLeaflets.nls$nbMax,
      nbSlope = model$nbLeaflets.nls$nbSlope,
      nbInfl= model$nbLeaflets.nls$nbInfl,
      nbLeaflets_SDP= round(model$nbLeaflets.nls$sigma),

      # parameters Leaflet radial angle
      leafletRadialHighA0Sup=
        model$radial.nls%>%dplyr::filter(Type=="High"&Mode=="sup")%>%
        dplyr::pull(A0),
      leafletRadialHighAmaxSup=
        model$radial.nls%>%dplyr::filter(Type=="High"&Mode=="sup")%>%
        dplyr::pull(Amax),
      leafletRadialHighA0Inf=
        model$radial.nls%>%dplyr::filter(Type=="High"&Mode=="inf")%>%
        dplyr::pull(A0),
      leafletRadialHighAmaxInf=
        model$radial.nls%>%dplyr::filter(Type=="High"&Mode=="inf")%>%
        dplyr::pull(Amax),
      leafletRadialLowA0Sup=
        model$radial.nls%>%dplyr::filter(Type=="Low"&Mode=="sup")%>%
        dplyr::pull(A0),
      leafletRadialLowAmaxSup=
        model$radial.nls%>%dplyr::filter(Type=="Low"&Mode=="sup")%>%
        dplyr::pull(Amax),
      leafletRadialLowA0Inf=
        model$radial.nls%>%dplyr::filter(Type=="Low"&Mode=="inf")%>%
        dplyr::pull(A0),
      leafletRadialLowAmaxInf=
        model$radial.nls%>%dplyr::filter(Type=="Low"&Mode=="inf")%>%
        dplyr::pull(Amax),

      # Frequency of leaflets type
      leafletFrequencyHigh= model$Rep[model$Rep$Type==1,]$Prop,
      leafletFrequencyLow= model$Rep[model$Rep$Type==-1,]$Prop,

      # leaflet shape
      xm_intercept= model$Shape$xm_intercept,
      xm_slope= model$Shape$xm_slope,
      ym_intercept= model$Shape$ym_intercept,
      ym_slope=  model$Shape$ym_slope,

      # ratio height/width sectionC
      heightWidthCRatio=
        mean(data$PetioleSectionC$Petiole_height_C_cm/
               data$PetioleSectionC$Petiole_width_C_cm),

      #------------------------------------PARAMETERS_PER_TREE---------------------------------------------

      # rachis length
      rachisLengthIntercept= model$rachisLength.lme$coef_simu[[1]][1],
      rachisLengthSlope= model$rachisLength.lme$coef_simu[[1]][2],
      rachisLength_SDP= model$rachisLength.lme$sigma,

      # C angle declination#
      decliCintercept= model$decliC.lme$coef_simu[[1]][1],
      decliCslope= model$decliC.lme$coef_simu[[1]][2],

      cPointAngle_SDP = model$decliC.lme$sigma,

      # Declination frond tip
      decMaxA= 180,
      decSlopeA= 0.01,
      decInflA= model$decliA_nls$decInflA,
      APointAngle_SDP = model$decliA_nls$sigma,

      # Rachis curvature
      coefCurve= rnorm(n=1,mean=mean(model$df_optim$coefCurv),
                       sd=sd(model$df_optim$coefCurv)),

      # Leaflets position
      coefDispo= model$dispo_nls$coefDispo,
      Dispo_SDP= model$dispo_nls$sigma,

      # Leaflet length b point
      lenfletLengthAtBIntercept= model$leaflet_length_B.lme$coef_simu[[1]][1],
      leafletLengthAtBSlope= model$leaflet_length_B.lme$coef_simu[[1]][2],

      # Parameters leaflet width b point
      bWidthIntercept= model$leaflet_width_B.lme$coef_simu[[1]][1],
      bWidthSlope= model$leaflet_width_B.lme$coef_simu[[1]][2],

      # Relative leaflet length
      lengthFirst= max(0,model$leafleftLength.nlme$coef_simu[[1]]['L0']),
      lengthLast= max(0,model$leafleftLength.nlme$coef_simu[[1]]['Lfin']),
      posLengthMax= model$leafleftLength.nlme$coef_simu[[1]]['Pos_Lmax'],

      # Relative leaflet width
      widthFirst= max(0,model$leafletWidth.nlme$coef_simu[[1]]['W0']),
      widthLast= max(0,model$leafletWidth.nlme$coef_simu[[1]]['Wfin']),
      posWidthMax= model$leafletWidth.nlme$coef_simu[[1]]['Pos_Wmax'],

      # Parameters Leaflet axial angle
      leafletAxialAngleC= model$axialAngle.nlme$coef_simu[[1]]['angleC'],
      leafletAxialAngleA= model$axialAngle.nlme$coef_simu[[1]]['angleA'],
      leafletAxialAngleSlope= model$axialAngle.nlme$coef_simu[[1]]['slopeC'],
      leafletAxialAngle_SDP= model$axialAngle.nlme$sigma,

      # Parameters width sectionC
      frondCpointWidthIntercept= model$petioleWidthC.lme$coef_simu[[1]]['(Intercept)'],
      frondCpointWidthSlope= model$petioleWidthC.lme$coef_simu[[1]]['CtoA'],

      # Parameters rachis relative height (mean + sd)
      rachisHeightTappering=
        model$rachisRelativeHeight.nlme$coef_mean +
        rnorm(n= 1, mean= 0, sd= model$rachisRelativeHeight.nlme$SdG1)
    )

  c(init,VP_list)
}

