#' Plot the architectural parameters estimations
#'
#' @description Plot the estimations of the palm architectural parameters from
#'  their parameter estimation and input data.
#'
#' @param data  Input data, generally one of the objects from [import_data()]
#' @param model Model fitted from one of [models][mod_stem_diameter()]
#' @param mode  Information level of the plot: "Progeny" for progeny scale, or "Tree"
#'              for tree scale.
#' @param decMaxA   The maximum allowed declination angle at the A point
#' @param decSlopeA The declination slope at point A
#' @param nb_leaves Number of leaves required for model prediction
#'
#' @section Plotting functions (plot_\*):
#'
#' * [plot_stem_diameter()]: plot the predicted stem diameter at the level of
#' the considered leaf using the leaf length
#'
#' * [plot_stem_height()]: plot the predicted stem height using the total cumulated
#'  emitted leaves
#'
#' * [plot_rachis_length()]: plot the predicted leaf length according to its index (leaf
#'  number)
#'
#' * [plot_petiole_ratio()]: plot the observed petiole to rachis ratio along the average
#' value used for simulation.
#'
#' * [plot_nb_leaflet()]: plot the predicted number of leaflets in function of rachis
#' length
#'
#' * [plot_C_declination()]: plot the declination angle at C point in function of leaf rank
#'
#' * [plot_leaf_curvature()]: plot the leaf curvature in function of the relative position
#'  on the leaf
#'
#' * [plot_A_declination()]: plot the leaf declination angle at A point in function of the C point
#'
#' * [plot_A_deviation()]: plot the leaf deviation computed using [deviation_A()], compared to the
#'  value used by Vpalm (the average value).
#'
#' * [plot_leaflet_position()]: plot the leaflet position along the rachis relative rank
#'  for each progeny with `mode= "Progeny"` or for each Progeny and colored by Tree with
#'  `mode= "Tree"`
#'
#' * [plot_leaflet_length_B()]: plot the predicted leaflet length at B point
#'
#' * [plot_leaflet_width_B()]: plot the predicted leaflet width at B point
#'
#' * [plot_leaflet_length()]: plot the predicted relative leaflet length
#'
#' * [plot_leaflet_width()]: plot the predicted leaflet width
#'
#' * [plot_leaflet_axial_angle()]: plot the predicted leaflet axial angle from the relative
#'  position of the leaflet on the rachis
#'
#' * [plot_leaflet_radial_angle()]: plot the leaflet radial angle from its relative position
#'   on the rachis
#'
#' * [plot_leaflet_type_freq()]: plot the leaflet radial angle frequency along the rachis.
#'  Leaflets can be positionned up, medium or down.
#'
#' * [plot_leaflet_shape()]: plot the leaflet shape
#'
#' * [plot_petiole_width_C()]: plot the petiole width at C point
#'
#' * [plot_petiole_height()]: plot the petiole heigth
#'
#' @importFrom dplyr group_by do mutate filter rename
#' @importFrom tidyr unnest
#' @importFrom ggplot2 ggplot facet_wrap geom_point geom_line ylab xlab aes guide_legend geom_ribbon facet_grid
#' @importFrom rlang .data
#'
#' @return A ggplot
#'
#'
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' archi_param= estimate_archi(progeny= "DY",map= 47)
#' plot_params(archi_param)
#'}
#' @export
#'
#'
plot_stem_diameter= function(data,model){

  pred=
    model%>%
    dplyr::ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(df=
                data.frame(RachisLength= seq(50,700,1),
                           StemDiameter=
                             sigmoid(X= seq(50,700,1), max= .$finalStemDiam,
                                     slope= .$StemDiamSlope, infl= .$StemDiamInfl)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(!is.na(.data$RachisLength)&!is.na(.data$StemDiameter))%>%
    ggplot2::ggplot(aes(x = .data$RachisLength, y= .data$StemDiameter, group= .data$Progeny))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point(aes(color= "Observed"))+
    ggplot2::geom_line(data=pred, aes(color= "Predicted"))+
    ggplot2::ylab('Stem Basis Diameter (cm)')+
    ggplot2::xlab('Rachis Length (cm)')
}


#' @rdname plot_stem_diameter
#' @export
plot_stem_height= function(data,model){
  pred=
    model%>%
    dplyr::ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(df=
                data.frame(TotalEmitted= seq(0,300,1),
                           StemHeight17=
                             stem_height(X= seq(0,300,1), y0= 5,
                                         coef= .$coefStemHeight)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(!is.na(.data$StemHeight17)&!is.na(.data$TotalEmitted))%>%
    ggplot(aes(x = .data$StemHeight17, y= .data$TotalEmitted, group= .data$Progeny))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point(aes(color= "Observed"))+
    ggplot2::geom_line(data=pred, aes(color= "Predicted"))+
    ggplot2::ylab('Stem height (cm)')+
    ggplot2::xlab('Number of emitted leaves from planting')
}


#' @rdname plot_stem_diameter
#' @export
plot_rachis_length= function(data,model,nb_leaves){
  physio=
    data%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::summarise(Physio_age= unique(.data$Physio_age))

  pred=
    model%>%
    merge(physio,by= "Progeny")%>%
    dplyr::ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(df=
                data.frame(LeafNumber= seq(.$Physio_age-nb_leaves,.$Physio_age,1),
                           RachisLength=
                             .$rachisLength_intercept +
                             seq(.$Physio_age-nb_leaves,.$Physio_age,1)*.$rachisLength_slope))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(!is.na(.data$RachisLength)&!is.na(.data$LeafNumber))%>%
    ggplot2::ggplot(aes(y= .data$RachisLength, x= .data$LeafNumber, color= "Observed"))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point()+
    ggplot2::geom_line(data= .data$pred, aes(color= "Predicted"))+
    ggplot2::ylab('Rachis Length (cm)')+ggplot2::xlab('Number of emmitted leaves from planting')+
    ggplot2::ggtitle(paste('Palm age=',paste(physio$Physio_age,collapse = ", "),
                           '\n(leaves emitted since planting)'))
}

#' @rdname plot_stem_diameter
#' @export
plot_petiole_ratio= function(data,model){
  data%>%
    dplyr::filter(.data$TotalEmitted<= .data$Physio_age+60 &
                    .data$TotalEmitted>= .data$Physio_age-60 &
                    .data$FrondRank>17)%>%
    dplyr::group_by(.data$Progeny)%>%
    merge(model, by= "Progeny")%>%
    ggplot2::ggplot(aes(x= .data$RachisLength))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point(aes(y= .data$RatioPetiole.x, color= "Observed"))+
    ggplot2::geom_line(aes(y= .data$RatioPetiole.y, color= "Simulated"))+
    ggplot2::labs(y= 'Petiol/Rachis length ratio', x= 'Rachis length (cm)')

}


#' @rdname plot_stem_diameter
#' @export
plot_B_position= function(data,model){
  data%>%
    dplyr::filter(.data$TotalEmitted<= .data$Physio_age+60 &
                    .data$TotalEmitted>= .data$Physio_age-60 &
                    .data$FrondRank>17 & !is.na(.data$PosB))%>%
    dplyr::group_by(.data$Progeny)%>%
    merge(model, by= "Progeny")%>%
    ggplot2::ggplot(aes(x= .data$RachisLength))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point(aes(y= .data$PosB.x, color= "Observed"))+
    ggplot2::geom_line(aes(y= .data$PosB.y, color= "Simulated"))+
    ggplot2::labs(y= 'Relative position of B point', x= 'Rachis length (cm)')
}



#' @rdname plot_stem_diameter
#' @export
plot_nb_leaflet= function(data,model){
  model=
    model$nbLeaflets.nls%>%
    dplyr::ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(df=
                data.frame(RachisLength= seq(50,700,1),
                           Nb_leaflets=
                             sigmoid(X= seq(50,700,1), max = .$nbMax,
                                     slope = .$nbSlope,
                                     infl = .$nbInfl)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(!is.na(.data$Nb_leaflets))%>%
    ggplot2::ggplot(aes(x = .data$RachisLength, y= .data$Nb_leaflets, group= .data$Progeny))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point(aes(color= "Observed"))+
    ggplot2::geom_line(data= model, aes(color= "Predicted"))+
    ggplot2::ylab('Number of leaflets')+
    ggplot2::xlab('Rachis Length (cm)')
}


#' @rdname plot_stem_diameter
#' @export
plot_C_declination= function(data,model){
  pred=
    model%>%
    ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(df=
                data.frame(Rank= seq(0,60,1),
                           Decli_C=
                             .$decliC_intercept+.$decliC_slope*seq(0,60,1)))%>%
    tidyr::unnest()

  data%>%
    ggplot2::ggplot(aes(y= .data$Decli_C, x= .data$Rank, color= "Observed"))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point()+
    ggplot2::geom_line(data= pred, aes(color= "Predicted"))+
    ggplot2::ylab('Declination at C point (degree)')+ggplot2::xlab('Rank')
}


#' @rdname plot_stem_diameter
#' @export
plot_leaf_curvature= function(data,model){

  df_plot=
    merge(data,model%>%select(-.data$value,-.data$conv),
          by = c('Progeny','TreeNumber'),all.x = T, sort = F)%>%
    dplyr::arrange(.data$Progeny, .data$TreeNumber, .data$Rank,
                   .data$RelativePositionRachisEstim)%>%
    dplyr::group_by(.data$Progeny, .data$TreeNumber, .data$Rank)%>%
    dplyr::mutate(
      angA_sim= sigmoid(X= .data$angC, max= 160, slope= 0.02,infl= .data$decAInfl),
      X_sim= leaf_curvature(position = .data$RelativePositionRachisEstim,
                            angC = .data$angC, angA = .data$angA_sim,
                            coefCurv = .data$coef_mean,
                            Length = .data$rachisLength)$X,
      Y_sim= leaf_curvature(position = .data$RelativePositionRachisEstim,
                            angC = .data$angC, angA = .data$angA_sim,
                            coefCurv = .data$coef_mean,
                            Length = .data$rachisLength)$Y)

  df_plot_sim=
    df_plot%>%
    dplyr::select(.data$Progeny, .data$TreeNumber, .data$Rank, .data$angC,
                  .data$angA_sim, .data$coef_mean, .data$rachisLength)%>%
    dplyr::summarise_all(mean)%>%
    dplyr::group_by(.data$Progeny, .data$TreeNumber, .data$Rank)%>%
    dplyr::do(leaf_curvature(position = seq(0,1,0.01),
                             angC = .$angC, angA = .$angA_sim,
                             coefCurv = .$coef_mean,
                             Length = .$rachisLength))

  ggplot2::ggplot(df_plot, aes(group= paste(.data$Progeny,.data$TreeNumber,.data$Rank)))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_line(aes(x = .data$X_distance_cm, y= .data$Y_distance_cm, color= "Observed"))+
    ggplot2::geom_point(aes(x = .data$X_distance_cm, y= .data$Y_distance_cm, color= "Observed"))+
    ggplot2::geom_line(aes(x = .data$X_sim, y= .data$Y_sim, color= "Predicted", lty= "Positions as observed"))+
    ggplot2::geom_line(data = df_plot_sim,
                       aes(x= .data$X, y= .data$Y,color= "Predicted", lty= "New simulated positions"))+
    ggplot2::ylab('Y distance (m)')+ggplot2::xlab('X distance (m)')+
    ggplot2::theme(legend.position = "bottom")+
    ggplot2::guides(color = guide_legend(title= "Colour: "),
                    linetype = guide_legend(title= "Line type: "))
}

#' @rdname plot_stem_diameter
#' @export
plot_A_declination= function(data, model, decMaxA= 180, decSlopeA= 0.01){
  pred=
    model%>%
    ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(angA_sim=
                data.frame(angC= seq(0,120,1),
                           angA= sigmoid(X=seq(0,120,1),max= decMaxA,
                                         slope = decSlopeA, infl= .$decInflA)))%>%
    tidyr::unnest()

  ggplot2::ggplot(data, aes(x = .data$angC, y= .data$angA, group= .data$Progeny))+
    ggplot2::facet_wrap(Progeny~.)+
    ggplot2::geom_point(aes(color= "Observed"))+
    ggplot2::geom_line(data= pred, aes(color= "Predicted"))+
    ggplot2::ylab('Declination at A point (degree)')+
    ggplot2::xlab('Declination at C point (degree)')
}

#' @rdname plot_stem_diameter
#' @export
plot_A_deviation= function(data,model){
  data%>%
    dplyr::filter(.data$Point=='A' & !is.na(.data$Z_distance_cm))%>%
    dplyr::group_by(.data$Progeny, .data$TreeNumber, .data$Rank)%>%
    dplyr::mutate(DevA_deg= deviation_A(devC = .data$Z_distance_cm,
                                        Length= .data$rachisLength))%>%
    merge(model%>%dplyr::rename(DevA_deg_sim= .data$DevA_deg),
          by = c('Progeny'),all.x = T, sort = F)%>%
    ggplot2::ggplot(aes(x= .data$rachisLength))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point(aes(y= .data$DevA_deg, color= "Observed"))+
    ggplot2::geom_line(aes(y= .data$DevA_deg_sim, color= "Predicted"))+
    ggplot2::ylab('Leaf deviation (deg)')+ggplot2::xlab('Rachis length (cm)')
}



#' @rdname plot_stem_diameter
#' @export
plot_leaflet_position= function(data, model, mode= c("Progeny","Tree")){
  mode= match.arg(mode,c("Progeny","Tree"))

  data=
    data%>%
    dplyr::group_by(.data$TreeNumber,.data$LeafIndex)%>%
    summarise(LeafNumber= min(.data$LeafNumber,na.rm = T))%>%
    merge(data,., by=c("TreeNumber","LeafIndex"), sort= F)%>%
    dplyr::filter(.data$PositionOnLeaflet==0)

  pred=
    model%>%
    ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(posrachis_sim=
                data.frame(RelativeLeafletRank= seq(0,1,0.01),
                           RelativePositionRachis=
                             leaflet_position(Rank= seq(0,1,0.01), param= .$coefDispo)))%>%
    tidyr::unnest()

  if(mode=="Progeny"){
    plot_pos=
      data%>%
      ggplot2::ggplot(aes(x= .data$RelativeLeafletRank, y= .data$RelativePositionRachis))+
      ggplot2::facet_wrap(.data$Progeny~.)+
      ggplot2::geom_point(aes(color= "Observed"))+
      ggplot2::geom_line(data=pred, aes(color= "Predicted"))+
      ggplot2::xlab('Leaflet relative rank')+ggplot2::ylab('Leafelt relative position')
  }
  if(mode=="Tree"){
    plot_pos=
      data%>%
      dplyr::group_by(.data$Progeny,.data$TreeNumber)%>%
      summarise(Tree_Index= unique(.data$TreeNumber))%>%
      dplyr::mutate(Tree_Index= seq_along(.data$TreeNumber))%>% # compute a Tree Index
      merge(data,., sort= F)%>%ungroup()%>%
      select(.data$Progeny, .data$Tree_Index, .data$LeafIndex,
             .data$TreeNumber, .data$LeafLength, .data$PositionOnRachis,
             .data$RelativeLeafletRank)%>%
      merge(pred%>%
              dplyr::rename(RelativeLeafletRank_sim= .data$RelativeLeafletRank,
                            RelativePositionRachis_sim= .data$RelativePositionRachis),
            by="Progeny")%>%
      dplyr::mutate(PositionRachis_sim= .data$RelativePositionRachis_sim*.data$LeafLength)%>%
      ggplot2::ggplot()+
      ggplot2::facet_wrap(.data$Progeny~.)+
      ggplot2::geom_point(aes(x= .data$RelativeLeafletRank, y= .data$PositionOnRachis, color= "Observed"))+
      ggplot2::geom_line(aes(x= .data$RelativeLeafletRank_sim, y= .data$PositionRachis_sim,
                             color= as.factor(.data$Tree_Index)))+
      # NB: if the true Tree number is needed, replace Tree_Index by TreeNumber above
      ggplot2::xlab('Leaflet relative rank')+ggplot2::ylab('Leafelt position on rachis (cm)')+
      ggplot2::labs(color= "Tree index \nin progenies")
  }
  plot_pos
}

#' @rdname plot_stem_diameter
#' @export
plot_leaflet_length_B= function(data, model){
  pred=
    model%>%
    ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(leaflet_length=
                data.frame(RachisLength= seq(50,700,1),
                           LeafletBLength=
                             .$Length_B_intercept+.$Length_B_slope*seq(50,700,1)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(!is.na(.data$LeafletBLength))%>%
    ggplot2::ggplot(aes(y=.data$LeafletBLength, x= .data$RachisLength, color= "Observed"))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point()+
    ggplot2::geom_line(data= pred, aes(color= "Predicted"))+
    ggplot2::ylab('Lealet length at B point (cm)')+ggplot2::xlab('Rachis length (cm)')

}


#' @rdname plot_stem_diameter
#' @export
plot_leaflet_width_B= function(data, model){
  pred=
    model%>%
    ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(leaflet_width=
                data.frame(RachisLength= seq(50,700,1),
                           LeafletBWidth=
                             .$width_B_intercept+.$width_B_slope*seq(50,700,1)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(!is.na(.data$LeafletBWidth))%>%
    ggplot2::ggplot(aes(y=.data$LeafletBWidth, x= .data$RachisLength, color= "Observed"))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point()+
    ggplot2::geom_line(data= pred, aes(color= "Predicted"))+
    ggplot2::ylab('Lealet width at B point (cm)')+ggplot2::xlab('Rachis length (cm)')
}



#' @rdname plot_stem_diameter
#' @export
plot_leaflet_length= function(data, model){
  pred=
    model$leafleftLength.nlme%>%
    ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(leaflet_length=
                data.frame(Position_rachis_rel= seq(0,1,0.01),
                           Relative_length=
                             leaflet_length(X= seq(0,1,0.01),Ymax=1,
                                            Y0= .$L0,Yfin= .$Lfin,
                                            X_Ymax= .$Pos_Lmax)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(.data$Width==0 & .data$PositionOnLeaflet!=0)%>%
    dplyr::mutate(Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength)%>%
    dplyr::group_by(.data$Progeny, .data$TreeNumber, .data$LeafIndex)%>%
    dplyr::mutate(Max_length= max(.data$PositionOnLeaflet))%>%
    ungroup()%>%
    dplyr::mutate(Relative_length= .data$PositionOnLeaflet/.data$Max_length)%>%
    dplyr::group_by(.data$Progeny)%>%
    ggplot2::ggplot(aes(y= .data$Relative_length, x= .data$Position_rachis_rel, color= "Observed"))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point()+
    ggplot2::geom_line(data= pred, aes(color= "Predicted"))+
    ggplot2::ylab('Relative leaflet length')+ggplot2::xlab('Relative position on rachis')
}


#' @rdname plot_stem_diameter
#' @export
plot_leaflet_width= function(data, model){
  pred=
    model%>%
    ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(leaflet_width=
                data.frame(Position_rachis_rel= seq(0,1,0.01),
                           Relative_max_width=
                             leaflet_max_width(X= seq(0,1,0.01), Ymax= 1,
                                               Y0= .$W0, Yfin= .$Wfin,
                                               X_Ymax= .$Pos_Wmax)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(.data$PositionOnLeaflet!=0)%>%
    dplyr::group_by(.data$Progeny, .data$TreeNumber, .data$LeafIndex, .data$Section)%>%
    dplyr::mutate(Leaflet_length= max(.data$PositionOnLeaflet),
                  Leaflet_max_width= max(.data$Width),
                  Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength,
                  Width_rel= .data$Width/.data$Leaflet_max_width)%>%
    dplyr::filter(.data$Width==0)%>%ungroup()%>%
    dplyr::group_by(.data$TreeNumber, .data$LeafIndex)%>%
    dplyr::mutate(Max_max_width= max(.data$Leaflet_max_width),
                  Relative_max_width= .data$Leaflet_max_width/.data$Max_max_width)%>%
    ggplot2::ggplot(aes(y=.data$Relative_max_width, x= .data$Position_rachis_rel, color= "Observed"))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point()+
    ggplot2::geom_line(data= pred, aes(color= "Predicted"))+
    ggplot2::ylab('Relative position on rachis')+
    ggplot2::xlab('Relative leaflet maximum width')
}

#' @rdname plot_stem_diameter
#' @export
plot_leaflet_axial_angle= function(data, model){
  pred=
    model%>%
    ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(leaflet_width=
                data.frame(Position_rel= seq(0,1,0.01),
                           Axial=
                             leaflet_axial_angle(position= seq(0,1,0.01), angle_C= .$angleC,
                                                 slope_C= .$slopeC, angle_A= .$angleA)))%>%
    tidyr::unnest()

  data%>%
    ggplot2::ggplot(aes(y=.data$Axial, x= .data$Position_rel, color= "Observed"))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point()+
    ggplot2::geom_line(data= pred, aes(color= "Predicted"))+
    ggplot2::ylab('Leaflet axial angle (deg)')+
    ggplot2::xlab('Relative position on rachis')
}

#' @rdname plot_stem_diameter
#' @export
plot_leaflet_radial_angle= function(data, model){
  data_polygons=
    model%>%
    ungroup()%>%
    dplyr::group_by(.data$Progeny, .data$Type, .data$Mode)%>%
    dplyr::do(df=
                data.frame(Position_rel= seq(0,1,0.01),
                           Radial_deg=
                             leaflet_radial_angle(position= seq(0,1,0.01),
                                                  A0= stats::coef(.$mod[[1]])['A0'],
                                                  Amax= stats::coef(.$mod[[1]])['Amax'],Xm=0.5)))%>%
    tidyr::unnest()%>%
    tidyr::spread(.data$Mode, .data$Radial_deg, sep="_")

  data_polygons=
    data_polygons%>%
    dplyr::filter(.data$Type=="High")%>%select(-.data$Mode_sup)%>%
    dplyr::rename(Mode_sup= .data$Mode_inf)%>%
    cbind(.,data_polygons%>%
            dplyr::filter(.data$Type=="Low")%>%select(.data$Mode_sup)%>%
            dplyr::rename(Mode_inf= .data$Mode_sup))%>%
    dplyr::mutate(Type= "Mid")%>%
    rbind(data_polygons,.)

  data%>%
    dplyr::mutate(Type= ifelse(.data$Type==0,"Mid",ifelse(.data$Type==1,"High","Low")))%>%
    ggplot2::ggplot()+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point(aes(y= .data$Radial_deg, x= .data$Position_rel, color= .data$Type))+
    geom_ribbon(data = data_polygons,
                aes(x= .data$Position_rel, ymin= .data$Mode_inf, ymax= .data$Mode_sup,
                    color= .data$Type, fill= .data$Type), alpha= 0.5)+
    ggplot2::labs(x= "Relative position on rachis", y= "Leaflet axial angle (deg)",
                  title= "Observed (points) and simulated range (polygons) of axial angles")
}


#' @rdname plot_stem_diameter
#' @export
plot_leaflet_type_freq= function(model){
  ggplot2::ggplot(model,aes(y= .data$Prop, x= .data$Position_rel, color= .data$Observation))+
    ggplot2::geom_line()+ggplot2::geom_point()+
    ggplot2::labs(y='Leaflets relative frequency',x='Relative position on rachis')
}

#' @rdname plot_stem_diameter
#' @export
plot_leaflet_shape= function(model){
  model%>%
    dplyr::group_by(.data$Progeny)%>%
    ggplot2::ggplot(aes(x= .data$PositionRelative))+
    facet_grid(.data$Progeny~.)+
    ggplot2::geom_point(aes(y= .data$xm, color= "xm"))+
    ggplot2::geom_point(aes(y= .data$ym, color= "ym"))+
    ggplot2::geom_smooth(aes(y= .data$xm, color= "xm"),method='lm',formula=y~x)+
    ggplot2::geom_smooth(aes(y= .data$ym, color= "ym"),method='lm',formula=y~x)+
    ggplot2::labs(y= "Leaflet shape")
}

#' @rdname plot_stem_diameter
#' @export
plot_petiole_width_C= function(data,model){
  pred_width.C=
    model%>%
    dplyr::ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(leaflet_length=
                data.frame(CtoA= seq(0,700,1),
                           Petiole_width_C_cm=
                             .$widthC_intercept+.$widthC_slope*seq(0,700,1)))%>%
    tidyr::unnest()

  data%>%
    ggplot2::ggplot(aes(y= .data$Petiole_width_C_cm, x= .data$CtoA, color= "Observed"))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point()+
    ggplot2::geom_line(data= pred_width.C, aes(color= "Predicted"))+
    ggplot2::ylab('Width section at C point (cm)')+ggplot2::xlab('Rachis length (cm)')
}


#' @rdname plot_stem_diameter
#' @export
plot_petiole_height= function(data,model){
  pred=
    model%>%
    dplyr::ungroup()%>%
    dplyr::group_by(.data$Progeny)%>%
    dplyr::do(rachis_height=
                data.frame(Position_rachis_rel= seq(0,1,0.1),
                           Petiole_relative_height=
                             section_height(Position = seq(0,1,0.1),
                                            a = .$coef.rachisRelativeHeight_mean)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(!is.na(.data$Petiole_relative_height))%>%
    ggplot2::ggplot(aes(y= .data$Petiole_relative_height, x= .data$Position_rachis_rel,
                        color= "Observed"))+
    ggplot2::facet_wrap(.data$Progeny~.)+
    ggplot2::geom_point()+
    ggplot2::geom_line(data= pred, aes(color= "Predicted"))+
    ggplot2::ylab('Petiole relative height')+
    ggplot2::xlab('Relative position on rachis')
}



#' Plot all model outputs
#'
#' @param data  A list of all data (generally from [import_data()])
#' @param model The models fitted on the data (generally from any [model][mod_stem_diameter()])
#' @param nb_leaves The number of leaves to plot
#'
#' @return A list of all plots
#' @export
#'
plot_all= function(data,model,nb_leaves= 45){
  list(
    stem_diameter= plot_stem_diameter(data = data$DataAll, model = model$StemDiam.nls),
    stem_height= plot_stem_height(data$DataAll,model$model.stemHeight),
    rachis_length=
      plot_rachis_length(data = data$DataAll,model = model$rachisLength.lme,
                         nb_leaves = nb_leaves),
    petiole_ratio= plot_petiole_ratio(data$DataAll,model$Pet),
    B_position= plot_B_position(data$DataAll,model$Bpos),
    nb_leaflet= plot_nb_leaflet(data$DataAll,model$nbLeaflets.nls),
    C_declination= plot_C_declination(data$declination, model$decliC.lme),
    leaf_curvature= plot_leaf_curvature(data= data$Curve, model = model$df_optim),
    A_declination= plot_A_declination(data$Curve,model$decliA_nls),
    A_deviation= plot_A_deviation(data$Curve,model$Dev),
    leaflet_position= plot_leaflet_position(data$DataAll,model$dispo_nls),
    leaflet_length_B= plot_leaflet_length_B(data$DataAll,model$leaflet_length_B.lme),
    leaflet_width_B= plot_leaflet_width_B(data$DataAll,model$leaflet_width_B.lme),
    leaflet_length= plot_leaflet_length(data$Area, model$leafleftLength.nlme),
    leaflet_width= plot_leaflet_width(data$Area, model$leafletWidth.nlme),
    leaflet_axial_angle= plot_leaflet_axial_angle(data$LftAngle, model$axialAngle.nlme),
    leaflet_radial_angle= plot_leaflet_radial_angle(data$LftAngle, model$radial.nls),
    leaflet_type_freq= plot_leaflet_type_freq(model$Rep),
    leaflet_shape= plot_leaflet_shape(model$Shape),
    petiole_width_C= plot_petiole_width_C(data$PetioleSectionC, model$petioleWidthC.lme),
    petiole_height= plot_petiole_height(data$RachisHeight, model$rachisRelativeHeight.nlme)
  )
}
