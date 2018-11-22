#' VPalm input
#'
#' @description Format the model outputs as a list for VPalm inputs
#'
#' @param data  A list of all data (generally from [import_data()])
#' @param model A list of the models fitted on the data (generally from
#'  [mod_all()])
#' @param nb_leaves The number of leaves needed
#' @param seed      The seed for random parameter generation
#'
#' @return A list of the parameters used as VPalm input
#' @export
#'
VPalm_list= function(data, model, nb_leaves= 45, seed= sample.int(1000,1)){
  list(
    seed= seed,

    # FIXED PARAMETERS PER PROGENY --------------------------------------------
    # Stem phylotaxis
    frondPhyllo_M= mean(data$Phylo$Phylo,na.rm=T),
    frondPhyllo_SD= sd(data$Phylo$Phylo,na.rm=T),

    # stem height
    residStemHeight= model$model.stemHeight$residStemHeight,

    # stem diameter
    stemDiamMax= model$StemDiam.nls$finalStemDiam,
    stemDiamSlope= model$StemDiam.nls$StemDiamSlope,
    stemDiamInfl= model$StemDiam.nls$StemDiamInfl,

    # number of frond
    nbFronds_M= nb_leaves,
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
      model$radial.nls%>%filter(Type=="High"&Mode=="sup")%>%
      pull(A0),
    leafletRadialHighAmaxSup=
      model$radial.nls%>%filter(Type=="High"&Mode=="sup")%>%
      pull(Amax),
    leafletRadialHighA0Inf=
      model$radial.nls%>%filter(Type=="High"&Mode=="inf")%>%
      pull(A0),
    leafletRadialHighAmaxInf=
      model$radial.nls%>%filter(Type=="High"&Mode=="inf")%>%
      pull(Amax),
    leafletRadialLowA0Sup=
      model$radial.nls%>%filter(Type=="Low"&Mode=="sup")%>%
      pull(A0),
    leafletRadialLowAmaxSup=
      model$radial.nls%>%filter(Type=="Low"&Mode=="sup")%>%
      pull(Amax),
    leafletRadialLowA0Inf=
      model$radial.nls%>%filter(Type=="Low"&Mode=="inf")%>%
      pull(A0),
    leafletRadialLowAmaxInf=
      model$radial.nls%>%filter(Type=="Low"&Mode=="inf")%>%
      pull(Amax),

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
    lengthFirst= max(0,model$leafleftLength.nlme$coef_simu[[1]][,'L0']),
    lengthLast= max(0,model$leafleftLength.nlme$coef_simu[[1]][,'Lfin']),
    posLengthMax= model$leafleftLength.nlme$coef_simu[[1]][,'Pos_Lmax'],

    # Relative leaflet width
    widthFirst= max(0,model$leafletWidth.nlme$coef_simu[[1]][,'W0']),
    widthLast= max(0,model$leafletWidth.nlme$coef_simu[[1]][,'Wfin']),
    posWidthMax= model$leafletWidth.nlme$coef_simu[[1]][,'Pos_Wmax'],

    # Parameters Leaflet axial angle
    leafletAxialAngleC= model$axialAngle.nlme$coef_simu[[1]][,'angleC'],
    leafletAxialAngleA= model$axialAngle.nlme$coef_simu[[1]][,'angleA'],
    leafletAxialAngleSlope= model$axialAngle.nlme$coef_simu[[1]][,'slopeC'],
    leafletAxialAngle_SDP= model$axialAngle.nlme$sigma,

    # Parameters width sectionC
    frondCpointWidthIntercept= model$petioleWidthC.lme$coef_simu[[1]][,'(Intercept)'],
    frondCpointWidthSlope= model$petioleWidthC.lme$coef_simu[[1]][,'CtoA'],

    # Parameters rachis relative height (mean + sd)
    rachisHeightTappering=
      model$rachisRelativeHeight.nlme$coef_mean +
      rnorm(n= 1, mean= 0, sd= model$rachisRelativeHeight.nlme$SdG1)
  )
}


write_VPalm=function(data){

  paramNames=c(
    paste('long seed = '),
    paste('int nbLeafEmitted = '),

    paste('double frondPhyllo_M = '),
    paste('double frondPhyllo_SD = '),

    paste('double H0 = '),
    paste('double coefStemHeight = '),
    paste('double residStemHeight= '),

    paste('double trunkBending_M = '),
    paste('double trunkBending_SD = '),

    paste('int nbFronds_M = '),
    paste('int nbFronds_SD = '),

    paste('double stemDiamMax = '),
    paste('double stemDiamSlope = '),
    paste('double stemDiamInfl = '),
    paste('double residStemDiam = '),

    paste('double decliCintercept ='),
    paste('double decliCslope ='),
    paste('double cPointAngle_SDP ='),

    paste('double decMaxA ='),
    paste('double decSlopeA ='),
    paste('double decInflA ='),
    paste('double APointAngle_SDP ='),

    paste('double rachisTwistFinalAngle_M ='),
    paste('double rachisTwistFinalAngle_SD ='),
    paste('double rachisTwistCoef ='),

    paste('double coefCurve ='),
    #paste('double coefCurve_SD ='),

    paste('double rachisDevFinalAngle_M ='),
    paste('double rachisDevFinalAngle_SD ='),
    paste('double rachisDevCoef ='),

    paste('double petioleRachisRatio_M ='),
    paste('double petioleRachisRatio_SD ='),

    paste('double rachisLengthIntercept ='),
    paste('double rachisLengthSlope ='),
    paste('double rachisLength_SDP ='),

    paste('double laminaAngle = '),

    paste('double nbMax ='),
    paste('double nbSlope ='),
    paste('double nbInfl ='),
    paste('double nbLeaflets_SDP ='),

    paste('double coefDispo= '),
    paste('double Dispo_SDP= '),
    paste('double pointBrelativePosition_M= '),
    paste('double pointBrelativePosition_SD= '),

    paste('double lenfletLengthAtBIntercept = '),
    paste('double leafletLengthAtBSlope = '),

    paste('double lengthFirst = '),
    paste('double lengthLast = '),
    paste('double posLengthMax = '),

    paste('double widthFirst = '),
    paste('double widthLast = '),
    paste('double posWidthMax = '),

    paste('double bWidthIntercept = '),
    paste('double bWidthSlope ='),

    paste('double xm_intercept= '),
    paste('double xm_slope= '),
    paste('double ym_intercept= '),
    paste('double ym_slope= '),

    paste('double leafletAxialAngleC = '),
    paste('double leafletAxialAngleA = '),
    paste('double leafletAxialAngleSlope = '),
    paste('double leafletAxialAngle_SDP = '),

    paste('double leafletStiffness = ' ),
    paste('double leafletStiffness_SD ='),

    paste('double[] leafletFrequencyHigh ='),
    paste('double[] leafletFrequencyLow =' ),

    paste('int nbInflorescences ='),

    paste('double frondBaseWidth = '),
    paste('double frondCpointWidthIntercept = '),
    paste('double frondCpointWidthSlope = '),
    paste('double frondtipWidth = '),
    paste('double frondBaseHeight = '),
    paste('double frondTipHeight = '),
    paste('double heightWidthCRatio = '),
    paste('double rachisHeightTappering = '),

    paste('double leafletRadialHighA0Sup= '),
    paste('double leafletRadialHighAmaxSup= '),
    paste('double leafletRadialHighA0Inf= '),
    paste('double leafletRadialHighAmaxInf= '),
    paste('double leafletRadialLowA0Sup= '),
    paste('double leafletRadialLowAmaxSup= '),
    paste('double leafletRadialLowA0Inf= '),
    paste('double leafletRadialLowAmaxInf= ')
  )

  paramValue=c(
    data$seed,
    data$nbLeafEmitted,
    data$frondPhyllo_M,
    data$frondPhyllo_SD,
    data$H0,
    data$coefStemHeight,
    data$residStemHeight,
    data$trunkBending_M,
    data$trunkBending_SD,
    data$nbFronds_M,
    data$nbFronds_SD,
    data$stemDiamMax,
    data$stemDiamSlope,
    data$stemDiamInfl,
    data$residStemDiam,
    data$decliCintercept,
    data$decliCslope,
    data$cPointAngle_SDP,
    data$decMaxA,
    data$decSlopeA,
    data$decInflA,
    data$APointAngle_SDP,
    data$rachisTwistFinalAngle_M,
    data$rachisTwistFinalAngle_SD,
    data$rachisTwistCoef,
    data$coefCurve,
    #coefCurve_SD,
    data$rachisDevFinalAngle_M,
    data$rachisDevFinalAngle_SD,
    data$rachisDevCoef,
    data$petioleRachisRatio_M,
    data$petioleRachisRatio_SD,
    data$rachisLengthIntercept,
    data$rachisLengthSlope,
    data$rachisLength_SDP,
    data$laminaAngle,
    data$nbMax,
    data$nbSlope,
    data$nbInfl,
    data$nbLeaflets_SDP,
    data$coefDispo,
    data$Dispo_SDP,
    data$pointBrelativePosition_M,
    data$pointBrelativePosition_SD,
    data$lenfletLengthAtBIntercept,
    data$leafletLengthAtBSlope,
    data$lengthFirst,
    data$lengthLast,
    data$posLengthMax,
    data$widthFirst,
    data$widthLast,
    data$posWidthMax,
    data$bWidthIntercept,
    data$bWidthSlope,
    data$xm_intercept,
    data$xm_slope,
    data$ym_intercept,
    data$ym_slope,
    data$leafletAxialAngleC,
    data$leafletAxialAngleA,
    data$leafletAxialAngleSlope,
    data$leafletAxialAngle_SDP,
    data$leafletStiffness,
    data$leafletStiffness_SD,
    paste(list(data$leafletFrequencyHigh)),
    paste(list(data$leafletFrequencyLow)),
    data$nbInflorescences,
    data$frondBaseWidth,
    data$frondCpointWidthIntercept,
    data$frondCpointWidthSlope,
    data$frondtipWidth,
    data$frondBaseHeight,
    data$frondTipHeight,
    data$heightWidthCRatio,
    data$rachisHeightTappering,
    data$leafletRadialHighA0Sup,
    data$leafletRadialHighAmaxSup,
    data$leafletRadialHighA0Inf,
    data$leafletRadialHighAmaxInf,
    data$leafletRadialLowA0Sup,
    data$leafletRadialLowAmaxSup,
    data$leafletRadialLowA0Inf,
    data$leafletRadialLowAmaxInf
  )


  output= as.matrix(paramValue)
  rownames(output)= paramNames

  return(output)
}
