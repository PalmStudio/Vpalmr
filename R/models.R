#' Fit the models
#'
#' @description These functions fit a given model from which the user can
#'    extract parameters used as VPalm inputs
#'
#' @param data      Input data, generally one of the objects from [import_data()]
#' @param decMaxA   The maximum allowed declination angle at the A point
#' @param decSlopeA The declination slope at point A
#'
#' @section Purpose of each function:
#'
#' * [mod_stem_diameter()]: model the stem diameter from the rachis length using a
#'  [sigmoid()]
#'
#' * [mod_stem_height()]: model the stem height using the total emitted
#' leaves throughout the palm life using [stem_height()]
#'
#' * [mod_rachis_length()]: model the leaf length from the leaf number (index of the
#'  leaves starting from the first leaf ever emitted)
#'
#' * [mod_petiole_ratio()]: model the petiole/rachis ratio by taking the average value for
#'  the progeny at the given physiological age +/- 60. If there is no value, return the
#'  average value by progeny for all leaves older than the reference leaf rank (>17).
#'
#' * [mod_nb_leaflet()]: model the number of leaflets from the rachis length using a
#'   [sigmoid()]
#'
#' * [mod_C_declination()]: model the leaf declination angle at C point from leaf rank
#'
#' * [mod_leaf_curvature()]: find the leaf inclination parameters (coefCurv and decAInfl)
#'  by minimizing the difference between a computed leaf using these parameters and the
#'  measurements. Uses [optim_Leaf_Curv_Infl()] under the hood.
#'
#' * [mod_A_declination()]: model the leaf declination angle at A point from C angle. This
#'  function needs the decAInfl parameter computed by [mod_leaf_curvature()], so it first
#'  checks if the parameter is present in the input data, and computes it if absent. The user
#'  should avoid this step by providing the parameter when possible because the
#'  [mod_leaf_curvature()] call can be computationally intensive.
#'
#' @section The models:
#'
#' * [stats::nls()] is used for [mod_stem_diameter()], [mod_stem_height()], [mod_nb_leaflet()]
#'
#' * [nlme::lme()] is used for [mod_rachis_length()], [mod_C_declination()]
#'
#' * [nlme::nlme()] is used for
#'
#' * [base::mean()] is used for [mod_petiole_ratio()]
#'
#' * [stats::optim()] is used by [mod_leaf_curvature()]
#'
#' @return A [tibble::tibble()] of the model outputs, eventually grouped by Progeny
#'  if several in input
#'
#'  @importFrom dplyr group_by do mutate filter
#'
#' @export
#'
mod_stem_diameter= function(data){
  data%>%
    dplyr::group_by(Progeny)%>%
    dplyr::do(mod=
                nls(data = .,
                    formula =
                      StemDiameter~sigmoid(X= RachisLength,max= finalStemDiam,
                                           slope= StemDiamSlope,infl= StemDiamInfl),
                    start= c(finalStemDiam= 300, StemDiamSlope= 0.01, StemDiamInfl= 200)))%>%
    dplyr::mutate(finalStemDiam= coef(mod)['finalStemDiam'],
                  StemDiamSlope= coef(mod)['StemDiamSlope'],
                  StemDiamInfl= coef(mod)['StemDiamInfl'],
                  residStemDiam=summary(mod)$sigma)
}


#' @rdname mod_stem_diameter
#' @export
mod_stem_height= function(data){
  data%>%
    dplyr::filter(!is.na(StemHeight17))%>%
    dplyr::group_by(Progeny)%>%
    dplyr::do(mod=
                nls(data = .,
                    formula =
                      StemHeight17 ~ Stem_height(X= TotalEmitted, y0= 5, coef= coefStem),
                    start= c(coefStem= 0.1)))%>%
    dplyr::mutate(coefStemHeight= coef(mod),
                  residStemHeight= summary(mod)$sigma)
}

#' @rdname mod_stem_diameter
#' @export
mod_rachis_length= function(data){
  data%>%
    merge(nbLeafEmitted)%>%
    filter(TotalEmitted<= Physio_age+30 & TotalEmitted>= Physio_age-30 &
             !is.na(RachisLength))%>%
    group_by(Progeny)%>%
    do(mod=
         nlme::lme(RachisLength ~ LeafNumber,
                   data= .,
                   random=~1+ LeafNumber |TreeNumber,method='ML',
                   control= param_controle))%>%
    mutate(rachisLength_intercept= summary(mod)$coefficients$fixed[1],
           rachisLength_slope= summary(mod)$coefficients$fixed[2],
           rachisLength_cov= list(vcov(mod)),
           coef.rachisLength_mean= list(c(rachisLength_intercept, rachisLength_slope)),
           SigmaR_rachisLength= summary(mod)$sigma,
           label.rachisLength.lme= list(colnames(summary(mod)$coefficients$random$TreeNumber)),
           SdG1_rachisLength= as.numeric(VarCorr(mod)['(Intercept)','StdDev']),
           SdG2_rachisLength= as.numeric(VarCorr(mod)['LeafNumber','StdDev']),
           corG_rachisLength= as.numeric(VarCorr(mod)['LeafNumber','Corr']),
           MatG_rachisLength= list(
             matrix(data=c(SdG1_rachisLength^2,
                           SdG1_rachisLength* SdG2_rachisLength* corG_rachisLength,
                           SdG1_rachisLength* SdG2_rachisLength* corG_rachisLength,
                           SdG2_rachisLength ^2),
                    nrow= length(coef.rachisLength_mean), ncol=length(coef.rachisLength_mean),
                    dimnames=list(label.rachisLength.lme, label.rachisLength.lme))))
}


#' @rdname mod_stem_diameter
#' @export
mod_petiole_ratio= function(data,Physio_age){

  Pet=
    data%>%
    group_by(Progeny)%>%
    filter(TotalEmitted<= Physio_age+60 &
             TotalEmitted>= Physio_age-60 &
             FrondRank>17)%>%
    summarise(RatioPetiole_age= mean(RatioPetiole, na.rm= T),
              RatioPetiole_age_sd= sd(RatioPetiole, na.rm= T))

  All_data=
    data%>%
    group_by(Progeny)%>%
    filter(FrondRank>17)%>%
    summarise(RatioPetiole_all= mean(RatioPetiole, na.rm= T),
              RatioPetiole_all_sd= sd(RatioPetiole, na.rm= T))

  out=
    merge(Pet,All_data,by= "Progeny")%>%
    mutate(RatioPetiole=
             ifelse(is.na(RatioPetiole_age),
                    RatioPetiole_all,
                    RatioPetiole_age),
           RatioPetiole_sd=
             ifelse(is.na(RatioPetiole_age),
                    RatioPetiole_all_sd,
                    RatioPetiole_age_sd),
           na_ratio=
             ifelse(is.na(RatioPetiole_age),
                    FALSE,TRUE))

  if(length(out$Progeny[!out$na_ratio])){
    warning("No petiol/rachis ratio data available for the physiological age requested for",
            " progenies: ", paste(out$Progeny[!out$na_ratio], collapse=", "),
            "using average from all leaf ranks > 17")
  }

  out%>%select(Progeny,RatioPetiole,RatioPetiole_sd)
}



#' @rdname mod_stem_diameter
#' @export
mod_B_position= function(data,Physio_age){

  Pet=
    data%>%
    group_by(Progeny)%>%
    filter(TotalEmitted<= Physio_age+60 &
             TotalEmitted>= Physio_age-60 &
             FrondRank>17 & !is.na(PosB))%>%
    summarise(PosB_age= mean(PosB, na.rm= T),
              PosB_age_sd= sd(PosB, na.rm= T))

  All_data=
    data%>%
    group_by(Progeny)%>%
    filter(FrondRank>17)%>%
    summarise(PosB_all= mean(PosB, na.rm= T),
              PosB_all_sd= sd(PosB, na.rm= T))

  out=
    merge(Pet,All_data,by= "Progeny")%>%
    mutate(PosB=
             ifelse(is.na(PosB_age),
                    PosB_all,
                    PosB_age),
           PosB_sd=
             ifelse(is.na(PosB_age),
                    PosB_all_sd,
                    PosB_age_sd),
           na_ratio=
             ifelse(is.na(PosB_age),
                    FALSE,TRUE))

  if(length(out$Progeny[!out$na_ratio])){
    warning("No relative position of B point data for the physiological age requested for",
            " progenies: ", paste(out$Progeny[!out$na_ratio], collapse=", "),
            "using average from all leaf ranks > 17")
  }

  out%>%select(Progeny,PosB,PosB_sd)
}


#' @rdname mod_stem_diameter
#' @export
mod_nb_leaflet= function(data){
  data%>%
    filter(!is.na(Nb_leaflets))%>%
    group_by(Progeny)%>%
    do(mod=
         nls(data = .,
             formula =
               Nb_leaflets ~ sigmoid(X= RachisLength,max= nbMax,slope= nbSlope,infl= nbInfl),
             start=list(nbMax =100,nbSlope=0.01,nbInfl= 140),
             upper=c(nbMax=180,nbSlope=10,nbInfl=300),
             lower= list(nbMax=0,nbSlope=-10,nbInfl=0),
             algorithm='port'))%>%
    mutate(nbMax= coef(mod)['nbMax'],
           nbSlope= coef(mod)['nbSlope'],
           nbInfl= coef(mod)['nbInfl'],
           SigmaR_nbLeaflets= summary(mod)$sigma,
           coef.nbLeaflets_mean= list(c(nbMax, nbSlope, nbInfl)),
           nbLeaflets_cov= list(vcov(mod)))
}


#' @rdname mod_stem_diameter
#' @export
mod_C_declination= function(data){
  data%>%
    filter(!is.na(Decli_C))%>%
    group_by(Progeny)%>%
    do(mod=
         lme(Decli_C ~ Rank,
             data= .,
             random= ~1 + Rank |TreeNumber,method='ML',
             control= param_controle))%>%
    mutate(decliC_intercept= summary(mod)$coefficients$fixed[1],
           decliC_slope= summary(mod)$coefficients$fixed[2],
           SigmaR_decliC= summary(mod)$sigma,
           SdG1_decli= as.numeric(VarCorr(mod)['(Intercept)','StdDev']),
           SdG2_decli= as.numeric(VarCorr(mod)['Rank','StdDev']),
           corG_decli= as.numeric(VarCorr(mod)['Rank','Corr']),
           decliC_cov= list(vcov(mod)),
           coef.decliC_mean= list(c(decliC_intercept, decliC_slope)),
           label.decliC.lme= list(colnames(summary(mod)$coefficients$random$TreeNumber)),
           MatG_decliC= list(
             matrix(data= c(SdG1_decli^2,
                            SdG1_decli*SdG2_decli* corG_decli,
                            SdG1_decli*SdG2_decli* corG_decli,SdG2_decli^2),
                    nrow=length(coef.decliC_mean),
                    ncol=length(coef.decliC_mean),
                    dimnames=list(label.decliC.lme, label.decliC.lme))))
}



#' @rdname mod_stem_diameter
#' @export
mod_leaf_curvature= function(data){
  data%>%
    group_by(Progeny, TreeNumber)%>%
    do(optim_Leaf_Curv_Infl(data= .))
}

#' @rdname mod_stem_diameter
#' @export
mod_A_declination= function(data, decMaxA= 180, decSlopeA= 0.01){

  # Checking if decAInfl was already computed, and if not, compute it:
  if(is.null(data$decAInfl)){
    data=
      merge(data,mod_leaf_curvature(data)%>%select(-value,-conv),
            by = c('Progeny','TreeNumber'),all.x = T, sort = F)%>%
      arrange(Progeny, TreeNumber, Rank, RelativePositionRachisEstim)
  }

  data%>%
    group_by(Progeny)%>%
    do(decliA_nls=
         nls(data = .,
             formula =
               angA~sigmoid(X= angC, max= decMaxA, slope= decSlopeA, infl= decAInfl),
             start=list(decAInfl= mean(.$decAInfl)),upper=c(decAInfl =120),
             lower= list(infl=0),algorithm='port'))%>%
    mutate(decInflA= coef(decliA_nls),
           SigmaR_decliA= summary(decliA_nls)$sigma)
}
