
#' Format VPalm inputs
#'
#' @description Format the [VPalm_list()] output to be written in a VPalm
#' input file
#'
#' @param data A [VPalm_list()] output
#'
#' @return A string in the VPalm input format
#' @export
#'
format_VPalm=function(data){

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
