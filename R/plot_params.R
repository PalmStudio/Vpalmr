#' Plot the architectural parameters estimations
#'
#' @description Plot the estimations of the palm architectural parameters from
#'  their parameter estimation and input data.
#'
#' @param data  Input data, generally one of the objects from [import_data()]
#' @param model Model fitted from one of [models][mod_stem_diameter()]
#' @param Physio_age Physiological age, in number of leaves emmitted from birth
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
#' * [plot_A_declination()]; plot the leaf declination angle at A point in function of the C point
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
                data.frame(RachisLength= rachisLengthSimu,
                           StemDiameter=
                             f.sigmo(X= rachisLengthSimu, max= .$finalStemDiam,
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
                data.frame(TotalEmitted= LeafNumberSimu,
                           StemHeight17=
                             Stem_height(X= LeafNumberSimu, y0= 5,
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
plot_rachis_length= function(data,model){
  pred=
    rachisLength.lme%>%
    right_join(nbLeafEmitted,by = "Progeny")%>%
    ungroup()%>%
    group_by(Progeny)%>%
    do(df=
         data.frame(LeafNumber= seq(.$Physio_age-nbFronds_M,.$Physio_age,1),
                    RachisLength=
                      .$rachisLength_intercept +
                      seq(.$Physio_age-nbFronds_M,.$Physio_age,1)*.$rachisLength_slope))%>%
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
         data.frame(RachisLength= rachisLengthSimu,
                    Nb_leaflets=
                      sigmoid(X= rachisLengthSimu, max = .$nbMax,
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
         data.frame(Rank= rankSimu,
                    Decli_C=
                      f.linear(X= rankSimu, intercept = .$decliC_intercept,
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
    do(leaf_curvature(position = relSimu,
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

