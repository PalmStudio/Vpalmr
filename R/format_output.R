#' Format all parameters lists
#'
#' @description Format all the lists from [extract_progeny()] output to prepare for writting.
#'
#' @details This function applies [format_tree()] sequentially to all trees from all
#' progenies and to the average list.
#'
#' @param data A [extract_progeny()] output
#'
#' @return A list of VPalm inputs in the [format_tree()] format.
#' @export
#'
#' @examples
#' \dontrun{
#' VPalm_out= Vpalmr::extract_progeny(data= Inputs, model= models, nb_leaves= 45, seed= 10)
#' out= format_list(data = VPalm_out)
#' }
#'
format_progeny= function(data){
  if(!is.null(data$Average)){
    out= lapply(data[-grep('Average',names(data))], function(x){lapply(x, format_tree)})
    out$Average= format_tree(data$Average)
  }else{
    out= lapply(data, function(x){lapply(x, format_tree)})
  }
  out
}

#' Format a VPalm input list
#'
#' @description Format an [extract_params()] list or one of the [extract_progeny()] list
#' to prepare for writting
#'
#' @param data An [extract_params()] list, or one from the [extract_progeny()] lists
#'
#' @return A [tibble::tibble()] with two columns:
#'
#' * name: the Java-formatted variable initialization for VPalm
#' * value: the variable value
#'
#' @export
#' @examples
#' \dontrun{
#' VPalm_out= Vpalmr::extract_progeny(data= Inputs, model= models, nb_leaves= 45, seed= 10)
#' average_params= format_param(data = VPalm_out$Average)
#' }
format_tree=function(data){

  paramNames=c(
    paste("Modelled Months After Planting = "),
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

  paramValue=
    list(
      data$MAP_requested,
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
      ifelse(!is.null(data$leafletFrequencyHigh)||
               !is.na(data$leafletFrequencyHigh),
             paste(list(data$leafletFrequencyHigh)),NULL),
      ifelse(!is.null(data$leafletFrequencyLow)||
               !is.na(data$leafletFrequencyLow),
             paste(list(data$leafletFrequencyLow)),NULL),
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

  # Error handling:
  no_val= unlist(lapply(paramValue,function(x)is.null(x)||is.na(x)))
  if(any(no_val)){
    stop("Missing value for VPalm parameter:",
         paste(gsub("=|double|long|int|\\[|\\]","",paramNames[no_val]),
               collapse = ", "))
  }

  tibble::tibble(name= paramNames, value= unlist(paramValue))

}
