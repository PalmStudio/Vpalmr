#' Fit the models
#'
#' @description These functions fit a given model from which the user can
#'    extract parameters used as VPalm inputs
#'
#' @param data      Input data, generally one of the objects from [import_data()]
#' @param decMaxA   The maximum allowed declination angle at the A point
#' @param decSlopeA The declination slope at point A
#' @param control   The control parameters for each kind of model
#' @param Area      The area data from [import_data()]
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
#' * [mod_A_deviation()]: model the leaf deviation by using the distance from C point
#'  and the length of the leaf between A and C point as inputs for [deviation_A()], and
#'  return the average value and standard deviation per progeny.
#'
#' * [mod_leaflet_position()]: model the leaflet position along the leaf rachis.
#'
#' * [mod_leaflet_length_B()]: model the leaflet length at B point
#'
#' * [mod_leaflet_width_B()]: model the leaflet width at B point
#'
#' * [mod_leaflet_length()]: model the leaflet relative length from its relative position
#'   on the rachis
#'
#' * [mod_leaflet_width()]:  model the leaflet relative width from its relative position
#'   on the rachis
#'
#' * [mod_leaflet_axial_angle()]: model the leaflet axial angle from its relative position
#'   on the rachis
#'
#' * [mod_leaflet_radial_angle()]: model the leaflet radial angle from its relative position
#'   on the rachis
#'
#' * [mod_leaflet_type_freq()]: model the leaflet radial angle frequency along the rachis.
#'  Leaflets can be positionned up, medium or down.
#'
#' * [mod_leaflet_shape()]: model the leaflet shape
#'
#' * [mod_petiole_width_C()]: model the petiole width at C point
#'
#' * [mod_petiole_height()]: model the petiole section height using [section_height()]
#'
#' @section The models:
#'
#' * [stats::nls()] is used for [mod_stem_diameter()], [mod_stem_height()],
#'   [mod_nb_leaflet()], [mod_leaflet_position()]
#'
#' * [nlme::lme()] is used for [mod_rachis_length()], [mod_C_declination()],
#'   [mod_leaflet_length_B()], [mod_petiole_width_C()], [mod_petiole_height()]
#'
#' * [nlme::nlme()] is used for [mod_leaflet_length()], [mod_leaflet_width()],
#'   [mod_leaflet_axial_angle()]
#'
#' * [base::mean()] is used for [mod_petiole_ratio()]
#'
#' * [stats::optim()] is used by [mod_leaf_curvature()]
#'
#' @return A [tibble::tibble()] of the model outputs, eventually grouped by Progeny
#'  if several in input
#'
#' @importFrom dplyr group_by do mutate filter rename ungroup
#' @importFrom nlme VarCorr
#' @importFrom rlang .data
#'
#' @export
#'
mod_stem_diameter= function(data){
  data%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nls(data = .,
                    formula =
                      StemDiameter~sigmoid(X= RachisLength,max= finalStemDiam,
                                           slope= StemDiamSlope,infl= StemDiamInfl),
                    start= c(finalStemDiam= 300, StemDiamSlope= 0.01, StemDiamInfl= 200)))%>%
    dplyr::mutate(finalStemDiam= stats::coef(.data$mod)['finalStemDiam'],
                  slope= stats::coef(.data$mod)['StemDiamSlope'],
                  Infl= stats::coef(.data$mod)['StemDiamInfl'],
                  sigma= summary(.data$mod)$sigma)
}


#' @rdname mod_stem_diameter
#' @export
mod_stem_height= function(data){
  data%>%
    dplyr::filter(!is.na(.data$StemHeight17))%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nls(data = .,
                    formula =
                      StemHeight17 ~ stem_height(X= TotalEmitted, y0= 5, coef= coefStem),
                    start= c(coefStem= 0.1)))%>%
    dplyr::mutate(coef= stats::coef(.data$mod),
                  sigma= summary(.data$mod)$sigma)
}

#' @rdname mod_stem_diameter
#' @export
mod_rachis_length= function(data,
                            control= nlme::lmeControl(maxIter=500000,
                                                      niterEM=25000)){
  data%>%
    filter(.data$TotalEmitted<= .data$Physio_age+30 &
             .data$TotalEmitted>= .data$Physio_age-30 &
             !is.na(.data$RachisLength))%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nlme::lme(RachisLength ~ LeafNumber,
                          data= .,
                          random=~1+ LeafNumber |TreeNumber,method='ML',
                          control= control))
}


#' @rdname mod_stem_diameter
#' @export
mod_petiole_ratio= function(data){

  Pet=
    data%>%
    group_by(.data$Progeny)%>%
    filter(.data$TotalEmitted<= .data$Physio_age+60 &
             .data$TotalEmitted>= .data$Physio_age-60 &
             .data$FrondRank>17)%>%
    summarise(RatioPetiole_age= mean(.data$RatioPetiole, na.rm= T),
              RatioPetiole_age_sd= stats::sd(.data$RatioPetiole, na.rm= T))

  All_data=
    data%>%
    group_by(.data$Progeny)%>%
    filter(.data$FrondRank>17)%>%
    summarise(RatioPetiole_all= mean(.data$RatioPetiole, na.rm= T),
              RatioPetiole_all_sd= stats::sd(.data$RatioPetiole, na.rm= T))

  out=
    merge(Pet,All_data,by= "Progeny")%>%
    mutate(RatioPetiole=
             ifelse(is.na(.data$RatioPetiole_age),
                    .data$RatioPetiole_all,
                    .data$RatioPetiole_age),
           RatioPetiole_sd=
             ifelse(is.na(.data$RatioPetiole_age),
                    .data$RatioPetiole_all_sd,
                    .data$RatioPetiole_age_sd),
           na_ratio=
             ifelse(is.na(.data$RatioPetiole_age),
                    FALSE,TRUE))

  if(length(out$Progeny[!out$na_ratio])){
    warning("No petiol/rachis ratio data available for the physiological age requested for",
            " progenies: ", paste(out$Progeny[!out$na_ratio], collapse=", "),
            "using average from all leaf ranks > 17")
  }

  out%>%select(.data$Progeny,.data$RatioPetiole,.data$RatioPetiole_sd)
}



#' @rdname mod_stem_diameter
#' @export
mod_B_position= function(data){

  Pet=
    data%>%
    group_by(.data$Progeny)%>%
    filter(.data$TotalEmitted<= .data$Physio_age+60 &
             .data$TotalEmitted>= .data$Physio_age-60 &
             .data$FrondRank>17 & !is.na(.data$PosB))%>%
    summarise(PosB_age= mean(.data$PosB, na.rm= T),
              PosB_age_sd= stats::sd(.data$PosB, na.rm= T))

  All_data=
    data%>%
    group_by(.data$Progeny)%>%
    filter(.data$FrondRank>17)%>%
    summarise(PosB_all= mean(.data$PosB, na.rm= T),
              PosB_all_sd= stats::sd(.data$PosB, na.rm= T))

  out=
    merge(Pet,All_data,by= "Progeny")%>%
    mutate(PosB=
             ifelse(is.na(.data$PosB_age),
                    .data$PosB_all,
                    .data$PosB_age),
           PosB_sd=
             ifelse(is.na(.data$PosB_age),
                    .data$PosB_all_sd,
                    .data$PosB_age_sd),
           na_ratio=
             ifelse(is.na(.data$PosB_age),
                    FALSE,TRUE))

  if(length(out$Progeny[!out$na_ratio])){
    warning("No relative position of B point data for the physiological age requested for",
            " progenies: ", paste(out$Progeny[!out$na_ratio], collapse=", "),
            "using average from all leaf ranks > 17")
  }

  out%>%select(.data$Progeny,.data$PosB,.data$PosB_sd)
}


#' @rdname mod_stem_diameter
#' @export
mod_nb_leaflet= function(data){
  data%>%
    filter(!is.na(.data$Nb_leaflets))%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nls(data = .,
                    formula =
                      Nb_leaflets ~ sigmoid(X= RachisLength,max= nbMax,slope= nbSlope,infl= nbInfl),
                    start=list(nbMax =100,nbSlope=0.01,nbInfl= 140),
                    upper=c(nbMax=180,nbSlope=10,nbInfl=300),
                    lower= list(nbMax=0,nbSlope=-10,nbInfl=0),
                    algorithm='port'))%>%
    mutate(nbMax= stats::coef(.data$mod)['nbMax'],
           nbSlope= stats::coef(.data$mod)['nbSlope'],
           nbInfl= stats::coef(.data$mod)['nbInfl'],
           sigma= summary(.data$mod)$sigma,
           coef_mean= list(stats::coef(.data$mod)),
           cov= list(stats::vcov(.data$mod)))
}


#' @rdname mod_stem_diameter
#' @export
mod_C_declination= function(data, control= nlme::lmeControl(maxIter=500000,
                                                            niterEM=25000)){
  data%>%
    filter(!is.na(.data$Decli_C))%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nlme::lme(Decli_C ~ Rank,
                          data= .,
                          random= ~1 + Rank |TreeNumber,method='ML',
                          control= control))
}



#' @rdname mod_stem_diameter
#' @export
mod_leaf_curvature= function(data){
  data%>%
    group_by(.data$Progeny, .data$TreeNumber)%>%
    dplyr::do(optim_Leaf_Curv_Infl(data= .))%>%
    dplyr::ungroup()
}

#' @rdname mod_stem_diameter
#' @export
mod_A_declination= function(data, decMaxA= 180, decSlopeA= 0.01){

  # Checking if decAInfl was already computed, and if not, compute it:
  if(is.null(data$decAInfl)){
    data=
      merge(data,mod_leaf_curvature(data)%>%
              select(-.data$value,-.data$conv),
            by = c('Progeny','TreeNumber'),all.x = T, sort = F)%>%
      arrange(.data$Progeny, .data$TreeNumber, .data$Rank,
              .data$RelativePositionRachisEstim)
  }

  data%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nls(data = .,
                    formula =
                      angA~sigmoid(X= angC, max= decMaxA, slope= decSlopeA, infl= decAInfl),
                    start=list(decAInfl= mean(.$decAInfl)),upper=c(decAInfl =120),
                    lower= list(infl=0),algorithm='port'))%>%
    mutate(decInflA= stats::coef(.data$mod),
           sigma= summary(.data$mod)$sigma)
}


#' @rdname mod_stem_diameter
#' @export
mod_A_deviation= function(data){
  data%>%
    filter(.data$Point=='A' & !is.na(.data$Z_distance_cm))%>%
    group_by(.data$Progeny,.data$TreeNumber,.data$Rank)%>%
    mutate(DevA_deg= deviation_A(devC = .data$Z_distance_cm,
                                 Length= .data$rachisLength))%>%
    ungroup()%>%
    group_by(.data$Progeny)%>%
    summarise(DevA_deg_sd= stats::sd(.data$DevA_deg),
              DevA_deg= mean(.data$DevA_deg))
}

#' @rdname mod_stem_diameter
#' @export
mod_leaflet_position= function(data, Area,
                               control= stats::nls.control(maxiter=500000)){
  data%>%
    group_by(.data$TreeNumber,.data$LeafIndex)%>%
    summarise(LeafNumber= min(.data$LeafNumber,na.rm = T))%>%
    merge(Area,., by=c("TreeNumber","LeafIndex"), sort= F)%>%
    filter(.data$PositionOnLeaflet==0)%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nls(data = .,
                    formula =
                      RelativePositionRachis~
                      leaflet_position(Rank= RelativeLeafletRank, param= coefDispo),
                    start= list(coefDispo=2),
                    control=control))%>%
    mutate(coefDispo= stats::coef(.data$mod),
           coefDispo_SD= summary(.data$mod)$parameters[,'Std. Error'],
           sigma= summary(.data$mod)$sigma)
}


#' @rdname mod_stem_diameter
#' @export
mod_leaflet_length_B= function(data,control= nlme::lmeControl(maxIter=500000,
                                                              niterEM=25000)){
  data%>%
    filter(!is.na(.data$LeafletBLength))%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nlme::lme(LeafletBLength ~ RachisLength,
                          data= .,
                          random= ~1 + RachisLength |TreeNumber,method='ML',
                          control= control))
}


#' @rdname mod_stem_diameter
#' @export
mod_leaflet_width_B= function(data, control= nlme::lmeControl(maxIter=500000,
                                                              niterEM=25000)){
  data%>%
    filter(!is.na(.data$LeafletBLength))%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nlme::lme(LeafletBWidth ~ RachisLength,
                          data= .,
                          random= ~1 + RachisLength |TreeNumber,method='ML',
                          control= control))
}


#' @rdname mod_stem_diameter
#' @export
mod_leaflet_length= function(data, control= nlme::nlmeControl(maxIter=500000,
                                                              niterEM=25000)){
  data%>%
    filter(.data$Width==0 & .data$PositionOnLeaflet!=0)%>%
    mutate(Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength)%>%
    group_by(.data$Progeny,.data$TreeNumber,.data$LeafIndex)%>%
    mutate(Max_length= max(.data$PositionOnLeaflet))%>%
    ungroup()%>%
    mutate(Relative_length= .data$PositionOnLeaflet/.data$Max_length)%>%
    group_by(.data$Progeny)%>%
    mutate(L0_obs=
             mean(.data$Relative_length[.data$Position_rachis_rel==
                                          min(.data$Position_rachis_rel, na.rm=T)]),
           Pos_Lmax_obs=
             mean(.data$Position_rachis_rel[.data$Relative_length==
                                              max(.data$Relative_length,na.rm=T)]))%>%
    dplyr::do(mod=
                nlme::nlme(Relative_length~
                             leaflet_length(X=Position_rachis_rel,
                                            Ymax=1,Y0= L0,Yfin=Lfin,X_Ymax=Pos_Lmax),
                           data= .,
                           start=
                             list(fixed= c(L0= mean(.$L0_obs),
                                           Lfin= mean(.$L0_obs),
                                           Pos_Lmax= mean(.$Pos_Lmax_obs))),
                           fixed= L0+Lfin+Pos_Lmax~1,
                           random= L0+Lfin+Pos_Lmax~1|TreeNumber,
                           control= control))
}



#' @rdname mod_stem_diameter
#' @export
mod_leaflet_width= function(data, control= nlme::nlmeControl(maxIter=500000,
                                                             niterEM=25000)){
  data%>%
    filter(.data$PositionOnLeaflet!=0)%>%
    group_by(.data$Progeny,.data$TreeNumber,.data$LeafIndex,.data$Section)%>%
    mutate(Leaflet_length= max(.data$PositionOnLeaflet),
           Leaflet_max_width= max(.data$Width),
           Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength,
           Width_rel= .data$Width/.data$Leaflet_max_width)%>%
    filter(.data$Width==0)%>%ungroup()%>%
    group_by(.data$TreeNumber,.data$LeafIndex)%>%
    mutate(Max_max_width= max(.data$Leaflet_max_width),
           Relative_max_width= .data$Leaflet_max_width/.data$Max_max_width)%>%
    group_by(.data$Progeny)%>%
    mutate(Wfin_obs=
             mean(.data$Relative_max_width[.data$Position_rachis_rel==
                                             max(.data$Position_rachis_rel,na.rm=T)]),
           Pos_Wmax_obs=
             mean(.data$Position_rachis_rel[.data$Relative_max_width==
                                              max(.data$Relative_max_width,na.rm=T)])
    )%>%
    dplyr::do(mod=
                nlme::nlme(Relative_max_width~
                             leaflet_max_width(X=Position_rachis_rel,Ymax=1,Y0=W0,
                                               Yfin=Wfin,X_Ymax=Pos_Wmax),
                           data= .,
                           start=
                             list(fixed= c(W0=0.2,
                                           Wfin= mean(.$Wfin_obs),
                                           Pos_Wmax= mean(.$Pos_Wmax_obs))),
                           fixed=W0+Wfin+Pos_Wmax~1,
                           random=W0+Wfin+Pos_Wmax~1|TreeNumber,
                           control= control))
}


#' @rdname mod_stem_diameter
#' @export
mod_leaflet_axial_angle= function(data,control= nlme::nlmeControl(maxIter=500000,
                                                                  niterEM=25000)){
  data%>%
    filter(!is.na(.data$Axial))%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nlme::nlme(Axial~
                             leaflet_axial_angle(position= Position_rel, angle_C=angleC,
                                                 slope_C= slopeC, angle_A= angleA),
                           data= .,
                           start= list(fixed=c(angleC=101,slopeC=-2,angleA=10)),
                           fixed= angleC+slopeC+angleA~1,
                           random= angleC+slopeC+angleA~1|TreeNumber,
                           control= control))
}



#' @rdname mod_stem_diameter
#' @export
mod_leaflet_radial_angle= function(data){
  dataAngle=
    data%>%
    group_by(.data$Progeny,.data$Type,.data$Section)%>%
    summarise(position= unique(.data$Section/10),
              mean= mean(.data$Radial_deg,na.rm=T),
              sd= stats::sd(.data$Radial_deg,na.rm=T))%>%
    mutate(up= .data$mean+2*.data$sd,
           dwn= .data$mean-2*.data$sd)%>%
    ungroup()

  radialHigh.nls=
    dataAngle%>%
    filter(.data$Type!=-1)%>%stats::na.omit()%>%
    group_by(.data$Progeny,.data$Type)%>%
    dplyr::do(mod=
                nls(formula =
                      up~
                      leaflet_radial_angle(position= position, A0= A0, Amax= Amax, Xm= 0.5),
                    data = .,start= list(A0= 10, Amax= max(.$up)), trace= F,
                    control= list(maxiter= 5000000, minFactor= 0.000000000001,
                                  warnOnly= T)))%>%
    mutate(Mode= ifelse(.data$Type==1,"sup","inf"), Type= "High")

  radialLow.nls=
    dataAngle%>%
    filter(.data$Type!=1)%>%stats::na.omit()%>%
    group_by(.data$Progeny,.data$Type)%>%
    dplyr::do(mod=
                nls(formula =
                      dwn~
                      leaflet_radial_angle(position= position, A0= A0, Amax= Amax, Xm= 0.5),
                    data = .,start= list(A0= -10, Amax= max(.$up)), trace= F,
                    control= list(maxiter= 5000000, minFactor= 0.000000000001, warnOnly= T)))%>%
    mutate(Mode= ifelse(.data$Type==0,"sup","inf"), Type= "Low")

  radial=
    rbind(radialHigh.nls,radialLow.nls)%>%
    arrange(.data$Progeny,.data$Type,.data$Mode)
  radial$A0= sapply(radial%>%.$mod,function(x){stats::coef(x)['A0']})
  radial$Amax= sapply(radial%>%.$mod,function(x){stats::coef(x)['Amax']})
  radial
}



#' @rdname mod_stem_diameter
#' @export
mod_leaflet_type_freq= function(data){
  Rep=
    table(data$Type, data$Section)%>%
    as.data.frame()
  colnames(Rep)= c('Type','Section','Freq')

  RepSection= as.data.frame(table(data$Section))
  colnames(RepSection)=c('Section','FreqTotal')

  Rep= cbind(Rep,FreqSection= rep(RepSection$FreqTotal,each=3))
  Rep$Position_rel= as.numeric(Rep$Section)/10
  Rep$Prop= Rep$Freq/Rep$FreqSection

  Rep$Observation= Rep$Type
  levels(Rep$Observation)= c("down",'medium',"up")
  Rep
}

#' @rdname mod_stem_diameter
#' @export
mod_leaflet_shape= function(data){
  Shape=
    data%>%
    group_by(.data$TreeNumber,.data$Section)%>%
    mutate(Leaflet_length= max(.data$PositionOnLeaflet),
           Leaflet_max_width= max(.data$Width),
           Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength,
           Position_leaflet_rel= .data$PositionOnLeaflet/.data$Leaflet_length,
           Width_rel= .data$Width/.data$Leaflet_max_width,
           PositionRelative= .data$Section/10)%>%
    summarise(
      Progeny= unique(.data$Progeny),
      PositionRelative= unique(.data$PositionRelative),
      shape=
        list(suppressWarnings(
          leaflet_shape_nls(.data$Width_rel,.data$Position_leaflet_rel))))

  Shape=
    data.frame(Shape%>%select(-.data$shape),
               unlist(Shape$shape)%>%matrix(ncol = 2, byrow = T))%>%
    rename(xm= .data$X1, ym= .data$X2)%>%
    group_by(.data$Progeny)%>%
    mutate(xm.lm= list(stats::lm(data=.,xm~PositionRelative)),
           ym.lm= list(stats::lm(data=.,ym~PositionRelative)))%>%
    ungroup()

  Shape$xm_intercept= sapply(Shape$xm.lm, function(x){stats::coef(x)[1]})
  Shape$xm_slope= sapply(Shape$xm.lm, function(x){stats::coef(x)[2]})
  Shape$ym_intercept= sapply(Shape$ym.lm, function(x){stats::coef(x)[1]})
  Shape$ym_slope= sapply(Shape$ym.lm, function(x){stats::coef(x)[2]})
  Shape%>%
    group_by(.data$Progeny)%>%
    summarise(xm_intercept= mean(.data$xm_intercept, na.rm= T),
              xm_slope= mean(.data$xm_slope, na.rm= T),
              ym_intercept= mean(.data$ym_intercept, na.rm= T),
              ym_slope= mean(.data$ym_slope, na.rm= T))
}

#' @rdname mod_stem_diameter
#' @export
mod_petiole_width_C= function(data, control= nlme::lmeControl(maxIter=500000,
                                                              niterEM=25000)){
  data%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nlme::lme(Petiole_width_C_cm~CtoA,
                          data= .,
                          random= ~1 + CtoA |TreeNumber,method='ML',
                          control= control))
}

#' @rdname mod_stem_diameter
#' @export
mod_petiole_height= function(data,control= nlme::nlmeControl(maxIter=500000,niterEM=25000)){
  data%>%
    filter(!is.na(.data$Petiole_relative_height))%>%
    group_by(.data$Progeny)%>%
    dplyr::do(mod=
                nlme::nlme(Petiole_relative_height~
                             section_height(Position= Position_rachis_rel, a),
                           data= .,
                           start= list(fixed= c(a= 0.1)),
                           fixed= a~1,
                           random= a~1|TreeNumber,
                           control= control))%>%
    mutate(coef_mean= nlme::fixed.effects(.data$mod),
           sigma= summary(.data$mod)$sigma,
           SdG1= nlme::intervals(.data$mod)$reStruct$TreeNumber[1,'est.'],
           cov= list(stats::vcov(.data$mod)))
}


#' Fit all models
#'
#' @param x  A list of all data (generally from [import_data()])
#'
#' @return A list of all model fits
#' @export
#'
mod_all= function(x){

  # STEM SCALE --------------------------------------------------------------

  # Stem diameter at the level of the considered leaf
  StemDiam.nls= mod_stem_diameter(x$DataAll)

  # Stem height
  model.stemHeight= mod_stem_height(x$DataAll)

  # LEAF SCALE --------------------------------------------------------------

  # Rachis length
  rachisLength.lme= mod_rachis_length(x$DataAll)

  # Ratio petiol/rachis
  Pet= mod_petiole_ratio(x$DataAll)

  # B point position
  Bpos= mod_B_position(x$DataAll)

  # Number of leaflets
  nbLeaflets.nls= mod_nb_leaflet(x$DataAll)

  # Declination at C point
  decliC.lme= mod_C_declination(x$declination)

  # Leaf curvature
  df_optim= mod_leaf_curvature(x$Curve)

  x$Curve=
    merge(x$Curve,df_optim%>%select(-.data$value,-.data$conv),
          by = c('Progeny','TreeNumber'),all.x = T, sort = F)%>%
    dplyr::arrange(.data$Progeny, .data$TreeNumber, .data$Rank,
                   .data$RelativePositionRachisEstim)

  # Declination at A point
  decliA_nls= mod_A_declination(x$Curve)

  # Rachis deviation
  Dev= mod_A_deviation(x$Curve)

  # LEAFLET SCALE -----------------------------------------------------------

  # Model the leaflet position along the leaf rachis
  dispo_nls= mod_leaflet_position(x$DataAll,x$Area)

  # Leaflet length at Bpoint
  leaflet_length_B.lme= mod_leaflet_length_B(x$DataAll)

  # Leaflet width at Bpoint
  leaflet_width_B.lme= mod_leaflet_width_B(x$DataAll)

  # Leaflet relative length from relative position on rachis
  leafleftLength.nlme= mod_leaflet_length(x$Area)

  # Leaflet relative widths
  leafletWidth.nlme= mod_leaflet_width(x$Area)

  # Leaflet axial angle
  axialAngle.nlme= mod_leaflet_axial_angle(x$LftAngle)

  # Leaflet radial angle (mean radial angle vs relative position on rachis)
  radial.nls= mod_leaflet_radial_angle(x$LftAngle)

  # Frequency of leaflets type depending on their radial angle
  Rep= mod_leaflet_type_freq(x$LftAngle)

  # Leaflet shape
  Shape= mod_leaflet_shape(x$Area)

  # NERVE SHAPE -------------------------------------------------------------

  # Nerve width at C point
  petioleWidthC.lme= mod_petiole_width_C(x$PetioleSectionC)

  # Nerve height
  rachisRelativeHeight.nlme= mod_petiole_height(x$RachisHeight)

  out=
    list(StemDiam.nls,model.stemHeight,rachisLength.lme,Pet,Bpos,nbLeaflets.nls,
         decliC.lme,df_optim,decliA_nls,Dev,dispo_nls,leaflet_length_B.lme,
         leaflet_width_B.lme,leafleftLength.nlme,leafletWidth.nlme,axialAngle.nlme,
         radial.nls,Rep,Shape,petioleWidthC.lme,rachisRelativeHeight.nlme)
  names(out)=
    c("StemDiam.nls","model.stemHeight","rachisLength.lme","Pet","Bpos","nbLeaflets.nls",
      "decliC.lme","df_optim","decliA_nls","Dev","dispo_nls","leaflet_length_B.lme",
      "leaflet_width_B.lme","leafleftLength.nlme","leafletWidth.nlme","axialAngle.nlme",
      "radial.nls","Rep","Shape","petioleWidthC.lme","rachisRelativeHeight.nlme")
  out
}
