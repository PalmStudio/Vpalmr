#' Fit the models
#'
#' @description These functions fit a given model from which the user can
#'    extract parameters used as VPalm inputs
#'
#' @param data      Input data, generally one of the objects from [import_data()]
#' @param decMaxA   The maximum allowed declination angle at the A point
#' @param decSlopeA The declination slope at point A
#' @param control   The control parameters for each kind of model
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
#' @importFrom dplyr group_by do mutate filter rename
#' @importFrom nlme VarCorr
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
                  slope= coef(mod)['StemDiamSlope'],
                  Infl= coef(mod)['StemDiamInfl'],
                  sigma= summary(mod)$sigma)
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
    dplyr::mutate(coef= coef(mod),
                  sigma= summary(mod)$sigma)
}

#' @rdname mod_stem_diameter
#' @export
mod_rachis_length= function(data,
                            control= nlme::lmeControl(maxIter=500000,
                                                      niterEM=25000)){
  data%>%
    filter(TotalEmitted<= Physio_age+30 & TotalEmitted>= Physio_age-30 &
             !is.na(RachisLength))%>%
    group_by(Progeny)%>%
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
mod_B_position= function(data){

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
    dplyr::do(mod=
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
           sigma= summary(mod)$sigma,
           coef_mean= list(c(nbMax, nbSlope, nbInfl)),
           cov= list(vcov(mod)))
}


#' @rdname mod_stem_diameter
#' @export
mod_C_declination= function(data, control= nlme::lmeControl(maxIter=500000,
                                                            niterEM=25000)){
  data%>%
    filter(!is.na(Decli_C))%>%
    group_by(Progeny)%>%
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
    group_by(Progeny, TreeNumber)%>%
    dplyr::do(optim_Leaf_Curv_Infl(data= .))
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
    dplyr::do(decliA_nls=
                nls(data = .,
                    formula =
                      angA~sigmoid(X= angC, max= decMaxA, slope= decSlopeA, infl= decAInfl),
                    start=list(decAInfl= mean(.$decAInfl)),upper=c(decAInfl =120),
                    lower= list(infl=0),algorithm='port'))%>%
    mutate(decInflA= coef(decliA_nls),
           sigma= summary(decliA_nls)$sigma)
}


#' @rdname mod_stem_diameter
#' @export
mod_A_deviation= function(data){
  data%>%
    filter(Point=='A' & !is.na(Z_distance_cm))%>%
    group_by(Progeny,TreeNumber,Rank)%>%
    mutate(DevA_deg= deviation_A(devC = Z_distance_cm,
                                 Length= rachisLength))%>%
    group_by(Progeny)%>%
    summarise(DevA_deg_sd= sd(DevA_deg),
              DevA_deg= mean(DevA_deg))
}

#' @rdname mod_stem_diameter
#' @export
mod_leaflet_position= function(data, Area,
                               control= stats::nls.control(maxiter=500000)){
  data%>%
    group_by(TreeNumber,LeafIndex)%>%
    summarise(LeafNumber= min(LeafNumber,na.rm = T))%>%
    merge(Area,., by=c("TreeNumber","LeafIndex"), sort= F)%>%
    filter(PositionOnLeaflet==0)%>%
    group_by(Progeny)%>%
    dplyr::do(dispo_nls=
                nls(data = .,
                    formula =
                      RelativePositionRachis~
                      leaflet_position(Rank= RelativeLeafletRank, param= coefDispo),
                    start= list(coefDispo=2),
                    control=control))%>%
    mutate(coefDispo= coef(dispo_nls),
           coefDispo_SD= summary(dispo_nls)$parameters[,'Std. Error'],
           sigma= summary(dispo_nls)$sigma)
}


#' @rdname mod_stem_diameter
#' @export
mod_leaflet_length_B= function(data,control= nlme::lmeControl(maxIter=500000,
                                                              niterEM=25000)){
  data%>%
    filter(!is.na(LeafletBLength))%>%
    group_by(Progeny)%>%
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
    filter(!is.na(LeafletBLength))%>%
    group_by(Progeny)%>%
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
    filter(Width==0 & PositionOnLeaflet!=0)%>%
    mutate(Position_rachis_rel= PositionOnRachis/LeafLength)%>%
    group_by(Progeny,TreeNumber,LeafIndex)%>%
    mutate(Max_length= max(PositionOnLeaflet))%>%
    ungroup()%>%
    mutate(Relative_length= PositionOnLeaflet/Max_length)%>%
    group_by(Progeny)%>%
    mutate(L0_obs=
             mean(Relative_length[Position_rachis_rel==
                                    min(Position_rachis_rel, na.rm=T)]),
           Pos_Lmax_obs=
             mean(Position_rachis_rel[Relative_length==
                                        max(Relative_length,na.rm=T)]))%>%
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
                           control= control))%>%
    mutate(L0= summary(mod)$coefficients$fixed['L0'],
           Lfin= summary(mod)$coefficients$fixed['Lfin'],
           Pos_Lmax= summary(mod)$coefficients$fixed['Pos_Lmax'],
           sigma= summary(mod)$sigma,
           SdG1= as.numeric(VarCorr(mod)['L0','StdDev']),
           SdG2= as.numeric(VarCorr(mod)['Lfin','StdDev']),
           SdG3= as.numeric(VarCorr(mod)['Pos_Lmax','StdDev']),
           corG12= as.numeric(VarCorr(mod)['Lfin','Corr']),
           corG13= as.numeric(VarCorr(mod)['Pos_Lmax','Corr']),
           corG23= as.numeric(VarCorr(mod)['Pos_Lmax',
                                           length(summary(mod)$coefficients$fixed)+1]),
           cov= list(vcov(mod)),
           coef_mean= list(c(L0,Lfin,Pos_Lmax)),
           label= list(colnames(summary(mod)$coefficients$random$TreeNumber)),
           MatG= list(
             matrix(data=c(SdG1^2, SdG1*SdG2*corG12, SdG1*SdG3*corG13,
                           SdG1*SdG2*corG12, SdG2^2, SdG2*SdG3*corG23,
                           SdG1*SdG3*corG13, SdG2*SdG3*corG23, SdG3^2),
                    nrow= length(coef_mean), ncol= length(coef_mean),
                    dimnames=list(label,label))))
}



#' @rdname mod_stem_diameter
#' @export
mod_leaflet_width= function(data, control= nlme::nlmeControl(maxIter=500000,
                                                             niterEM=25000)){
  data%>%
    filter(PositionOnLeaflet!=0)%>%
    group_by(Progeny,TreeNumber,LeafIndex,Section)%>%
    mutate(Leaflet_length= max(PositionOnLeaflet),
           Leaflet_max_width= max(Width),
           Position_rachis_rel= PositionOnRachis/LeafLength,
           Width_rel= Width/Leaflet_max_width)%>%
    filter(Width==0)%>%ungroup()%>%
    group_by(TreeNumber,LeafIndex)%>%
    mutate(Max_max_width= max(Leaflet_max_width),
           Relative_max_width= Leaflet_max_width/Max_max_width,
           leafID= paste(TreeNumber,'Leaf',LeafIndex))%>%
    group_by(Progeny)%>%
    mutate(Wfin_obs=
             mean(Relative_max_width[Position_rachis_rel==
                                       max(Position_rachis_rel,na.rm=T)]),
           Pos_Wmax_obs=
             mean(Position_rachis_rel[Relative_max_width==
                                        max(Relative_max_width,na.rm=T)])
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
                           control= control))%>%
    mutate(W0=summary(mod)$coefficients$fixed['W0'],
           Wfin=summary(mod)$coefficients$fixed['Wfin'],
           Pos_Wmax=summary(mod)$coefficients$fixed['Pos_Wmax'],
           Sigma= summary(mod)$sigma,
           SdG1= as.numeric(VarCorr(mod)['W0','StdDev']),
           SdG2= as.numeric(VarCorr(mod)['Wfin','StdDev']),
           SdG3= as.numeric(VarCorr(mod)['Pos_Wmax','StdDev']),
           corG12= as.numeric(VarCorr(mod)['Wfin','Corr']),
           corG13= as.numeric(VarCorr(mod)['Pos_Wmax','Corr']),
           corG23=
             as.numeric(VarCorr(mod)['Pos_Wmax',
                                     length(summary(mod)$coefficients$fixed)+1]),
           cov= list(vcov(mod)),
           coef_mean= list(c(W0,Wfin,Pos_Wmax)),
           label= list(colnames(summary(mod)$coefficients$random$TreeNumber)),
           MatG= list(
             matrix( data= c(SdG1^2, SdG1*SdG2*corG12, SdG1*SdG3*corG13,
                             SdG1*SdG2*corG12, SdG2^2, SdG2*SdG3*corG23,
                             SdG1*SdG3*corG13, SdG2*SdG3*corG23, SdG3^2),
                     nrow= length(coef_mean), ncol= length(coef_mean),
                     dimnames= list(label,label))))
}


#' @rdname mod_stem_diameter
#' @export
mod_leaflet_axial_angle= function(data,control= nlme::nlmeControl(maxIter=500000,
                                                                  niterEM=25000)){
  data%>%
    filter(!is.na(Axial))%>%
    group_by(Progeny)%>%
    dplyr::do(mod=
                nlme::nlme(Axial~
                             leaflet_axial_angle(position= Position_rel, angle_C=angleC,
                                                 slope_C= slopeC, angle_A= angleA),
                           data= .,
                           start= list(fixed=c(angleC=101,slopeC=-2,angleA=10)),
                           fixed= angleC+slopeC+angleA~1,
                           random= angleC+slopeC+angleA~1|TreeNumber,
                           control= control))%>%
    mutate(angleC= summary(mod)$coefficients$fixed['angleC'],
           slopeC= summary(mod)$coefficients$fixed['slopeC'],
           angleA= summary(mod)$coefficients$fixed['angleA'],
           sigma= summary(mod)$sigma,
           SdG1 =as.numeric(VarCorr(mod)['angleC','StdDev']),
           SdG2=as.numeric(VarCorr(mod)['slopeC','StdDev']),
           SdG3 =as.numeric(VarCorr(mod)['angleA','StdDev']),
           corG12=as.numeric(VarCorr(mod)['slopeC','Corr']),
           corG13=as.numeric(VarCorr(mod)['angleA','Corr']),
           corG23=
             as.numeric(VarCorr(mod)['angleA',
                                     length(summary(mod)$coefficients$fixed)+1]),
           cov= list(vcov(mod)),
           coef_mean= list(c(angleC,slopeC, angleA)),
           label= list(colnames(summary(mod)$coefficients$random$TreeNumber)),
           MatG= list(
             matrix(
               data= c(SdG1^2, SdG1*SdG2*corG12, SdG1*SdG3*corG13,
                       SdG1*SdG2*corG12, SdG2^2, SdG2*SdG3*corG23,
                       SdG1*SdG3*corG13, SdG2*SdG3*corG23, SdG3^2),
               nrow=length(coef_mean), ncol=length(coef_mean),
               dimnames=list(label, label))))
}



#' @rdname mod_stem_diameter
#' @export
mod_leaflet_radial_angle= function(data){
  dataAngle=
    data%>%
    group_by(Progeny,Type,Section)%>%
    summarise(position= unique(Section/10),
              mean= mean(Radial_deg,na.rm=T),
              sd= sd(Radial_deg,na.rm=T))%>%
    mutate(up= mean+2*sd,
           dwn= mean-2*sd)

  radialHigh.nls=
    dataAngle%>%filter(Type!=-1)%>%na.omit()%>%group_by(Progeny,Type)%>%
    dplyr::do(mod=
                nls(formula =
                      up~
                      leaflet_radial_angle(position= position, A0= A0, Amax= Amax, Xm= 0.5),
                    data = .,start= list(A0= 10, Amax= max(.$up)), trace= F,
                    control= list(maxiter= 5000000, minFactor= 0.000000000001,
                                  warnOnly= T)))%>%
    mutate(Mode= ifelse(Type==1,"sup","inf"),Type= "High")

  radialLow.nls=
    dataAngle%>%filter(Type!=1)%>%na.omit()%>%group_by(Progeny,Type)%>%
    dplyr::do(mod=
                nls(formula =
                      dwn~
                      leaflet_radial_angle(position= position, A0= A0, Amax= Amax, Xm= 0.5),
                    data = .,start= list(A0= -10, Amax= max(.$up)), trace= F,
                    control= list(maxiter= 5000000, minFactor= 0.000000000001, warnOnly= T)))%>%
    mutate(Mode= ifelse(Type==0,"sup","inf"),Type= "Low")

  radial= rbind(radialHigh.nls,radialLow.nls)%>%arrange(Progeny,Type,Mode)
  radial$A0= sapply(radial%>%.$mod,function(x){coef(x)['A0']})
  radial$Amax= sapply(radial%>%.$mod,function(x){coef(x)['Amax']})
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
    group_by(TreeNumber,Section)%>%
    mutate(Leaflet_length= max(PositionOnLeaflet),
           Leaflet_max_width= max(Width),
           Position_rachis_rel= PositionOnRachis/LeafLength,
           Position_leaflet_rel= PositionOnLeaflet/Leaflet_length,
           Width_rel= Width/Leaflet_max_width,
           PositionRelative= Section/10)%>%
    summarise(
      Progeny= unique(Progeny),
      PositionRelative= unique(PositionRelative),
      shape= list(suppressWarnings(leaflet_shape_nls(Width_rel,Position_leaflet_rel))))

  Shape=
    data.frame(Shape%>%select(-shape),
               unlist(Shape$shape)%>%matrix(ncol = 2, byrow = T))%>%
    rename(xm= X1, ym= X2)%>%
    group_by(Progeny)%>%
    mutate(xm.lm= list(lm(data=.,xm~PositionRelative)),
           ym.lm= list(lm(data=.,ym~PositionRelative)))

  Shape$xm_intercept= sapply(Shape$xm.lm, function(x){coef(x)[1]})
  Shape$xm_slope= sapply(Shape$xm.lm, function(x){coef(x)[2]})
  Shape$ym_intercept= sapply(Shape$ym.lm, function(x){coef(x)[1]})
  Shape$ym_slope= sapply(Shape$ym.lm, function(x){coef(x)[2]})
  Shape%>%
    group_by(Progeny)%>%
    summarise(xm_intercept= mean(xm_intercept, na.rm= T),
              xm_slope= mean(xm_slope, na.rm= T),
              ym_intercept= mean(ym_intercept, na.rm= T),
              ym_slope= mean(ym_slope, na.rm= T))
}

#' @rdname mod_stem_diameter
#' @export
mod_petiole_width_C= function(data, control= nlme::lmeControl(maxIter=500000,
                                                              niterEM=25000)){
  data%>%
    group_by(Progeny)%>%
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
    filter(!is.na(Petiole_relative_height))%>%
    group_by(Progeny)%>%
    dplyr::do(mod=
                nlme::nlme(Petiole_relative_height~
                             section_height(Position= Position_rachis_rel, a),
                           data= .,
                           start= list(fixed= c(a= 0.1)),
                           fixed= a~1,
                           random= a~1|TreeNumber,
                           control= control))%>%
    mutate(coef_mean= summary(mod)$coefficients$fixed,
           sigma= summary(mod)$sigma,
           SdG1= nlme::intervals(mod)$reStruct$TreeNumber[1,'est.'],
           cov= list(vcov(mod)))
}


#' Fit all models
#'
#' @param data  A list of all data (generally from [import_data()])
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
  decliC.lme= mod_C_declination(x$Dec)

  # Leaf curvature
  df_optim= mod_leaf_curvature(x$Curve)

  x$Curve=
    merge(x$Curve,df_optim%>%select(-value,-conv),
          by = c('Progeny','TreeNumber'),all.x = T, sort = F)%>%
    dplyr::arrange(Progeny, TreeNumber, Rank, RelativePositionRachisEstim)

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
