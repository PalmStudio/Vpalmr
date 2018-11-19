


# RV: attention, nbLeafEmitted est maintenant un data.frame ! Adapter l'input !
write_VPalm=function(Progeny,nbLeafEmitted,seed){
  #-----------------------------------FIXED_PARAMETERS-----------------------------------------------


  # -----PARAMETERS_FILE-----

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
    List$seed,
    List$nbLeafEmitted,
    List$frondPhyllo_M,
    List$frondPhyllo_SD,
    List$H0,
    List$coefStemHeight,
    List$residStemHeight,
    List$trunkBending_M,
    List$trunkBending_SD,
    List$nbFronds_M,
    List$nbFronds_SD,
    List$stemDiamMax,
    List$stemDiamSlope,
    List$stemDiamInfl,
    List$residStemDiam,
    List$decliCintercept,
    List$decliCslope,
    List$cPointAngle_SDP,
    List$decMaxA,
    List$decSlopeA,
    List$decInflA,
    List$APointAngle_SDP,
    List$rachisTwistFinalAngle_M,
    List$rachisTwistFinalAngle_SD,
    List$rachisTwistCoef,
    List$coefCurve,
    #coefCurve_SD,
    List$rachisDevFinalAngle_M,
    List$rachisDevFinalAngle_SD,
    List$rachisDevCoef,
    List$petioleRachisRatio_M,
    List$petioleRachisRatio_SD,
    List$rachisLengthIntercept,
    List$rachisLengthSlope,
    List$rachisLength_SDP,
    List$laminaAngle,
    List$nbMax,
    List$nbSlope,
    List$nbInfl,
    List$nbLeaflets_SDP,
    List$coefDispo,
    List$Dispo_SDP,
    List$pointBrelativePosition_M,
    List$pointBrelativePosition_SD,
    List$lenfletLengthAtBIntercept,
    List$leafletLengthAtBSlope,
    List$lengthFirst,
    List$lengthLast,
    List$posLengthMax,
    List$widthFirst,
    List$widthLast,
    List$posWidthMax,
    List$bWidthIntercept,
    List$bWidthSlope,
    List$xm_intercept,
    List$xm_slope,
    List$ym_intercept,
    List$ym_slope,
    List$leafletAxialAngleC,
    List$leafletAxialAngleA,
    List$leafletAxialAngleSlope,
    List$leafletAxialAngle_SDP,
    List$leafletStiffness,
    List$leafletStiffness_SD,
    paste(list(List$leafletFrequencyHigh)),
    paste(list(List$leafletFrequencyLow)),
    List$nbInflorescences,
    List$frondBaseWidth,
    List$frondCpointWidthIntercept,
    List$frondCpointWidthSlope,
    List$frondtipWidth,
    List$frondBaseHeight,
    List$frondTipHeight,
    List$heightWidthCRatio,
    List$rachisHeightTappering,
    List$leafletRadialHighA0Sup,
    List$leafletRadialHighAmaxSup,
    List$leafletRadialHighA0Inf,
    List$leafletRadialHighAmaxInf,
    List$leafletRadialLowA0Sup,
    List$leafletRadialLowAmaxSup,
    List$leafletRadialLowA0Inf,
    List$leafletRadialLowAmaxInf
  )


  output=as.matrix(paramValue)
  rownames(output)= paramNames

  return(output)

}
