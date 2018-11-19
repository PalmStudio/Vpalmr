#' Plot the architectural parameters estimations
#'
#' @description Plot the estimations of the palm architectural parameters from
#'  their parameter estimation and input data.
#'
#' @param data  Input data, generally one of the objects from [import_data()]
#' @param model Model fitted from one of [models][mod_stem_diameter()]
#' @param Physio_age Physiological age, in number of leaves emmitted from birth
#' @param mode  Information level of the plot: "Progeny" for progeny scale, or "Tree"
#'              for tree scale.
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
#'  @importFrom dplyr group_by do mutate filter
#'  @importFrom tidyr unnest
#'  @importFrom ggplot2 ggplot facet_wrap geom_point geom_line ylab xlab aes
#'
#' @return A ggplot
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
    dplyr::group_by(Progeny)%>%
    dplyr::do(df=
                data.frame(RachisLength= seq(50,700,1),
                           StemDiameter=
                             f.sigmo(X= seq(50,700,1), max= .$finalStemDiam,
                                     slope= .$StemDiamSlope, infl= .$StemDiamInfl)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(!is.na(RachisLength)&!is.na(StemDiameter))%>%
    ggplot2::ggplot(aes(x = RachisLength, y= StemDiameter, group= Progeny))+
    ggplot2::facet_wrap(Progeny~.)+
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
    dplyr::group_by(Progeny)%>%
    dplyr::do(df=
                data.frame(TotalEmitted= seq(0,300,1),
                           StemHeight17=
                             Stem_height(X= seq(0,300,1), y0= 5,
                                         coef= .$coefStemHeight)))%>%
    tidyr::unnest()

  data%>%
    dplyr::filter(!is.na(StemHeight17)&!is.na(TotalEmitted))%>%
    ggplot(aes(x = StemHeight17, y= TotalEmitted, group= Progeny))+
    ggplot2::facet_wrap(Progeny~.)+
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
    group_by(Progeny)%>%
    summarise(Physio_age= unique(Physio_age))

  pred=
    model%>%
    merge(physio,by= "Progeny")%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(df=
         data.frame(LeafNumber= seq(.$Physio_age-nb_leaves,.$Physio_age,1),
                    RachisLength=
                      .$rachisLength_intercept +
                      seq(.$Physio_age-nb_leaves,.$Physio_age,1)*.$rachisLength_slope))%>%
    unnest()

  data%>%
    dplyr::filter(!is.na(RachisLength)&!is.na(LeafNumber))%>%
    ggplot(aes(y= RachisLength, x= LeafNumber, color= "Observed"))+
    facet_wrap(Progeny~.)+
    geom_point()+
    geom_line(data= pred, aes(color= "Predicted"))+
    ylab('Rachis Length (cm)')+xlab('Number of emmitted leaves from planting')+
    ggtitle(paste('Palm age=',paste(nbLeafEmitted$Physio_age,collapse = ", "),
                  '\n(leaves emitted since planting)'))
}

#' @rdname plot_stem_diameter
#' @export
plot_petiole_ratio= function(data,model,Physio_age){
  data%>%
    filter(TotalEmitted<= Physio_age+60 &
             TotalEmitted>= Physio_age-60 &
             FrondRank>17)%>%
    group_by(Progeny)%>%
    merge(model, by= "Progeny")%>%
    ggplot(aes(x= RachisLength))+
    facet_wrap(Progeny~.)+
    geom_point(aes(y= RatioPetiole.x, color= "Observed"))+
    geom_line(aes(y= RatioPetiole.y, color= "Simulated"))+
    labs(y= 'Petiol/Rachis length ratio', x= 'Rachis length (cm)')

}


#' @rdname plot_stem_diameter
#' @export
plot_B_position= function(data,model,Physio_age){
  data%>%
    filter(TotalEmitted<= Physio_age+60 &
             TotalEmitted>= Physio_age-60 &
             FrondRank>17 & !is.na(PosB))%>%
    group_by(Progeny)%>%
    merge(model, by= "Progeny")%>%
    ggplot(aes(x= RachisLength))+
    facet_wrap(Progeny~.)+
    geom_point(aes(y= PosB.x, color= "Observed"))+
    geom_line(aes(y= PosB.y, color= "Simulated"))+
    labs(y= 'Relative position of B point', x= 'Rachis length (cm)')
}



#' @rdname plot_stem_diameter
#' @export
plot_nb_leaflet= function(data,model){
  model=
    nbLeaflets.nls%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(df=
         data.frame(RachisLength= seq(50,700,1),
                    Nb_leaflets=
                      sigmoid(X= seq(50,700,1), max = .$nbMax,
                              slope = .$nbSlope,
                              infl = .$nbInfl)))%>%
    unnest()

  data%>%
    filter(!is.na(Nb_leaflets))%>%
    ggplot(aes(x = RachisLength, y= Nb_leaflets, group= Progeny))+
    facet_wrap(Progeny~.)+
    geom_point(aes(color= "Observed"))+
    geom_line(data= model, aes(color= "Predicted"))+
    ylab('Number of leaflets')+
    xlab('Rachis Length (cm)')
}


#' @rdname plot_stem_diameter
#' @export
plot_C_declination= function(data,model){
  pred=
    model%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(df=
         data.frame(Rank= seq(0,60,1),
                    Decli_C=
                      f.linear(X= seq(0,60,1), intercept = .$decliC_intercept,
                               slope = .$decliC_slope)))%>%
    unnest()

  data%>%
    ggplot(aes(y= Decli_C, x= Rank, color= "Observed"))+
    facet_wrap(Progeny~.)+
    geom_point()+
    geom_line(data= pred, aes(color= "Predicted"))+
    ylab('Declination at C point (degree)')+xlab('Rank')
}


#' @rdname plot_stem_diameter
#' @export
plot_leaf_curvature= function(data,model){

  df_plot=
    merge(data,model%>%select(-value,-conv),
          by = c('Progeny','TreeNumber'),all.x = T, sort = F)%>%
    arrange(Progeny, TreeNumber, Rank, RelativePositionRachisEstim)%>%
    group_by(Progeny, TreeNumber, Rank)%>%
    mutate(
      angA_sim= sigmoid(X= angC, max= 160, slope= 0.02,infl= decAInfl),
      X_sim= leaf_curvature(position = RelativePositionRachisEstim,
                            angC = angC, angA = angA_sim,
                            coefCurv = coefCurv,
                            Length = rachisLength)$X,
      Y_sim= leaf_curvature(position = RelativePositionRachisEstim,
                            angC = angC, angA = angA_sim,
                            coefCurv = coefCurv,
                            Length = rachisLength)$Y)

  df_plot_sim=
    df_plot%>%
    select(Progeny, TreeNumber,Rank,angC,angA_sim,
           coefCurv,rachisLength)%>%
    summarise_all(mean)%>%
    group_by(Progeny, TreeNumber,Rank)%>%
    do(leaf_curvature(position = seq(0,1,0.01),
                      angC = .$angC, angA = .$angA_sim,
                      coefCurv = .$coefCurv,
                      Length = .$rachisLength))

  ggplot(df_plot, aes(group= paste(Progeny,TreeNumber,Rank)))+
    facet_wrap(Progeny~.)+
    geom_line(aes(x = X_distance_cm, y= Y_distance_cm, color= "Observed"))+
    geom_point(aes(x = X_distance_cm, y= Y_distance_cm, color= "Observed"))+
    geom_line(aes(x = X_sim, y= Y_sim, color= "Predicted", lty= "Positions as observed"))+
    geom_line(data = df_plot_sim,
              aes(x= X, y= Y,color= "Predicted", lty= "New simulated positions"))+
    ylab('Y distance (m)')+xlab('X distance (m)')+
    theme(legend.position = "bottom")+
    guides(color = guide_legend(title= "Colour: "),
           linetype = guide_legend(title= "Line type: "))
}

#' @rdname plot_stem_diameter
#' @export
plot_A_declination= function(data,model){
  pred=
    model%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(angA_sim=
         data.frame(angC= seq(0,120,1),
                    angA= sigmoid(X=seq(0,120,1),max= decMaxA,
                                  slope = decSlopeA, infl= .$decInflA)))%>%
    unnest()

  ggplot(data, aes(x = angC, y= angA, group= Progeny))+
    facet_wrap(Progeny~.)+
    geom_point(aes(color= "Observed"))+
    geom_line(data= pred, aes(color= "Predicted"))+
    ylab('Declination at A point (degree)')+
    xlab('Declination at C point (degree)')
}

#' @rdname plot_stem_diameter
#' @export
plot_A_deviation= function(data,model){
  data%>%
    filter(Point=='A' & !is.na(Z_distance_cm))%>%
    group_by(Progeny,TreeNumber,Rank)%>%
    mutate(DevA_deg= deviation_A(devC = Z_distance_cm,
                                 Length= rachisLength))%>%
    merge(model%>%rename(DevA_deg_sim= DevA_deg),
          by = c('Progeny'),all.x = T, sort = F)%>%
    ggplot(aes(x= rachisLength))+
    facet_wrap(Progeny~.)+
    geom_point(aes(y= DevA_deg, color= "Observed"))+
    geom_line(aes(y= DevA_deg_sim, color= "Predicted"))+
    ylab('Leaf deviation (deg)')+xlab('Rachis length (cm)')
}



#' @rdname plot_stem_diameter
#' @export
plot_leaflet_position= function(data, model, mode= c("Progeny","Tree")){
  mode= match.arg(mode,c("Progeny","Tree"))

  data=
    data%>%
    group_by(TreeNumber,LeafIndex)%>%
    summarise(LeafNumber= min(LeafNumber,na.rm = T))%>%
    merge(Area,., by=c("TreeNumber","LeafIndex"), sort= F)%>%
    filter(PositionOnLeaflet==0)

  pred=
    model%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(posrachis_sim=
         data.frame(RelativeLeafletRank= seq(0,1,0.01),
                    RelativePositionRachis=
                      leaflet_position(Rank= seq(0,1,0.01), param= .$coefDispo)))%>%
    unnest()

  if(mode=="Progeny"){
    plot_pos=
    data%>%
      ggplot(aes(x= RelativeLeafletRank, y= RelativePositionRachis))+
      facet_wrap(Progeny~.)+
      geom_point(aes(color= "Observed"))+
      geom_line(data=pred, aes(color= "Predicted"))+
      xlab('Leaflet relative rank')+ylab('Leafelt relative position')
  }
  if(mode=="Tree"){
    plot_pos=
    data%>%
      group_by(Progeny,TreeNumber)%>%
      summarise(Tree_Index= unique(TreeNumber))%>%
      mutate(Tree_Index= seq_along(TreeNumber))%>% # compute a Tree Index
      merge(data,., sort= F)%>%ungroup()%>%
      select(Progeny,Tree_Index,LeafIndex,TreeNumber,
             LeafLength,PositionOnRachis,RelativeLeafletRank)%>%
      merge(pred%>%
              rename(RelativeLeafletRank_sim= RelativeLeafletRank,
                     RelativePositionRachis_sim= RelativePositionRachis), by="Progeny")%>%
      mutate(PositionRachis_sim= RelativePositionRachis_sim * LeafLength)%>%
      ggplot()+
      facet_wrap(Progeny~.)+
      geom_point(aes(x= RelativeLeafletRank, y= PositionOnRachis, color= "Observed"))+
      geom_line(aes(x= RelativeLeafletRank_sim, y= PositionRachis_sim,
                    color= as.factor(Tree_Index)))+
      # NB: if the true Tree number is needed, replace Tree_Index by TreeNumber above
      xlab('Leaflet relative rank')+ylab('Leafelt position on rachis (cm)')+
      labs(color= "Tree index \nin progenies",
           title= paste('Palm age=',paste(nbLeafEmitted$Physio_age, collapse = ", ")),
           subtitle= "Physiologial age in number of leaves emitted since planting")
  }
  plot_pos
}

#' @rdname plot_stem_diameter
#' @export
plot_leaflet_length_B= function(data, model){
  pred=
    model%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(leaflet_length=
         data.frame(RachisLength= seq(50,700,1),
                    LeafletBLength=
                      .$Length_B_intercept+.$Length_B_slope*seq(50,700,1)))%>%
    unnest()

  data%>%
    filter(!is.na(LeafletBLength))%>%
    ggplot(aes(y=LeafletBLength, x= RachisLength, color= "Observed"))+
    facet_wrap(Progeny~.)+
    geom_point()+
    geom_line(data= pred, aes(color= "Predicted"))+
    ylab('Lealet length at B point (cm)')+xlab('Rachis length (cm)')

}


#' @rdname plot_stem_diameter
#' @export
plot_leaflet_width_B= function(data, model){
  pred=
    model%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(leaflet_width=
         data.frame(RachisLength= seq(50,700,1),
                    LeafletBWidth=
                      .$width_B_intercept+.$width_B_slope*seq(50,700,1)))%>%
    unnest()

  data%>%
    filter(!is.na(LeafletBWidth))%>%
    ggplot(aes(y=LeafletBWidth, x= RachisLength, color= "Observed"))+
    facet_wrap(Progeny~.)+
    geom_point()+
    geom_line(data= pred, aes(color= "Predicted"))+
    ylab('Lealet width at B point (cm)')+xlab('Rachis length (cm)')
}



#' @rdname plot_stem_diameter
#' @export
plot_leaflet_length= function(data, model){
  pred=
    leafleftLength.nlme%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(leaflet_length=
         data.frame(Position_rachis_rel= seq(0,1,0.01),
                    Relative_length=
                      leaflet_length(X= seq(0,1,0.01),Ymax=1,
                                     Y0= .$L0,Yfin= .$Lfin,
                                     X_Ymax= .$Pos_Lmax)))%>%
    unnest()

  data%>%
    filter(Width==0 & PositionOnLeaflet!=0)%>%
    mutate(Position_rachis_rel= PositionOnRachis/LeafLength)%>%
    group_by(Progeny,TreeNumber,LeafIndex)%>%
    mutate(Max_length= max(PositionOnLeaflet))%>%
    ungroup()%>%
    mutate(Relative_length= PositionOnLeaflet/Max_length)%>%
    group_by(Progeny)%>%
    ggplot(aes(y= Relative_length, x= Position_rachis_rel, color= "Observed"))+
    facet_wrap(Progeny~.)+
    geom_point()+
    geom_line(data= pred, aes(color= "Predicted"))+
    ylab('Relative leaflet length')+xlab('Relative position on rachis')+
    ggtitle(paste('Palm age=',paste(nbLeafEmitted$Physio_age, collapse = ", "),
                  'leaves emitted since planting'))
}


#' @rdname plot_stem_diameter
#' @export
plot_leaflet_width= function(data, model){
  pred=
    model%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(leaflet_width=
         data.frame(Position_rachis_rel= seq(0,1,0.01),
                    Relative_max_width=
                      leaflet_max_width(X= seq(0,1,0.01), Ymax= 1,
                                        Y0= .$W0, Yfin= .$Wfin,
                                        X_Ymax= .$Pos_Wmax)))%>%
    unnest()

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
           Relative_max_width= Leaflet_max_width/Max_max_width)%>%
    ggplot(aes(y=Relative_max_width, x= Position_rachis_rel, color= "Observed"))+
    facet_wrap(Progeny~.)+
    geom_point()+
    geom_line(data= pred, aes(color= "Predicted"))+
    ylab('Relative position on rachis')+
    xlab('Relative leaflet maximum width')+
    ggtitle(c(paste('Palm age=',paste(nbLeafEmitted$Physio_age, collapse = ", "),
                    'leaves emitted since planting')))
}

#' @rdname plot_stem_diameter
#' @export
plot_leaflet_axial_angle= function(data, model){
  pred=
    model%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(leaflet_width=
         data.frame(Position_rel= seq(0,1,0.01),
                    Axial=
                      f.axialAngle(X= seq(0,1,0.01), angleC= .$angleC,
                                   slopeC= .$slopeC, angleA= .$angleA)))%>%
    unnest()

  data%>%
    ggplot(aes(y=Axial, x= Position_rel, color= "Observed"))+
    facet_wrap(Progeny~.)+
    geom_point()+
    geom_line(data= pred, aes(color= "Predicted"))+
    ylab('Leaflet axial angle (deg)')+
    xlab('Relative position on rachis')+
    ggtitle(c(paste('Palm age=',paste(nbLeafEmitted$Physio_age, collapse = ", "),
                    'leaves emitted since planting')))
}

#' @rdname plot_stem_diameter
#' @export
plot_leaflet_radial_angle= function(data, model){
  data_polygons=
    model%>%
    ungroup()%>%
    group_by(Progeny,Type,Mode)%>%
    do(df=
         data.frame(Position_rel= seq(0,1,0.01),
                    Radial_deg=
                      f.radialAngle(x= seq(0,1,0.01), A0= coef(.$mod[[1]])['A0'],
                                    Amax= coef(.$mod[[1]])['Amax'],Xm=0.5)))%>%
    unnest()%>%
    tidyr::spread(Mode,Radial_deg, sep="_")

  data_polygons=
    data_polygons%>%
    filter(Type=="High")%>%select(-Mode_sup)%>%rename(Mode_sup= Mode_inf)%>%
    cbind(.,data_polygons%>%
            filter(Type=="Low")%>%select(Mode_sup)%>%rename(Mode_inf= Mode_sup))%>%
    mutate(Type= "Mid")%>%
    rbind(data_polygons,.)

  data%>%
    mutate(Type= ifelse(Type==0,"Mid",ifelse(Type==1,"High","Low")))%>%
    ggplot()+
    facet_wrap(Progeny~.)+
    geom_point(aes(y= Radial_deg, x= Position_rel, color= Type))+
    geom_ribbon(data = data_polygons,
                aes(x= Position_rel, ymin= Mode_inf, ymax= Mode_sup,
                    color= Type, fill= Type), alpha= 0.5)+
    labs(x= "Relative position on rachis", y= "Leaflet axial angle (deg)",
         title= "Observed (points) and simulated range (polygons) of axial angles")
}


#' @rdname plot_stem_diameter
#' @export
plot_leaflet_type_freq= function(model){
  ggplot(model,aes(y= Prop, x= Position_rel, color= Observation))+
    geom_line()+geom_point()+
    labs(y='Leaflets relative frequency',x='Relative position on rachis')
}

#' @rdname plot_stem_diameter
#' @export
plot_leaflet_shape= function(model){
  model%>%
    group_by(Progeny)%>%
    ggplot(aes(x= PositionRelative))+
    facet_grid(Progeny~.)+
    geom_point(aes(y= xm, color= "xm"))+
    geom_point(aes(y= ym, color= "ym"))+
    geom_smooth(aes(y= xm, color= "xm"),method='lm',formula=y~x)+
    geom_smooth(aes(y= ym, color= "ym"),method='lm',formula=y~x)+
    labs(y= "Leaflet shape")
}

#' @rdname plot_stem_diameter
#' @export
plot_petiole_width_C= function(data,model){
  pred_width.C=
    model%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(leaflet_length=
         data.frame(CtoA= seq(0,700,1),
                    Petiole_width_C_cm=
                      .$widthC_intercept+.$widthC_slope*seq(0,700,1)))%>%
    unnest()

  data%>%
    ggplot(aes(y= Petiole_width_C_cm, x= CtoA, color= "Observed"))+
    facet_wrap(Progeny~.)+
    geom_point()+
    geom_line(data= pred_width.C, aes(color= "Predicted"))+
    ylab('Width section at C point (cm)')+xlab('Rachis length (cm)')
}


#' @rdname plot_stem_diameter
#' @export
plot_petiole_height= function(data,model){
  pred=
    model%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(rachis_height=
         data.frame(Position_rachis_rel= seq(0,1,0.1),
                    Petiole_relative_height=
                      section_height(Position = seq(0,1,0.1),
                                     a = .$coef.rachisRelativeHeight_mean)))%>%
    unnest()

  data%>%
    filter(!is.na(Petiole_relative_height))%>%
    ggplot(aes(y= Petiole_relative_height, x= Position_rachis_rel, color= "Observed"))+
    facet_wrap(Progeny~.)+
    geom_point()+
    geom_line(data= pred, aes(color= "Predicted"))+
    ylab('Petiole relative height')+
    xlab('Relative position on rachis')+
    ggtitle(c(paste('Palm age=',paste(nbLeafEmitted$Physio_age, collapse = ", "),
                    'leaves emitted since planting')))
}


