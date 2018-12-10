#' Import all data
#'
#' @description Import and pre-compute the input data for VPalm
#'
#' @param parameter     File path to the parameter file (ParameterSimu.csv)
#' @param development   File path to the development file (Development_Rep4_SMSE.csv)
#' @param phylotaxy     File path to the phylotaxy file (Stem_SMSE14.csv)
#' @param declination   File path to the declination file (A_SMSE_Nov14.csv)
#' @param curvature     File path to the leaf curvature file (LeafCurvature_SMSE14.csv)
#' @param torsion       File path to the leaf torsion file (AnglesC&A_SMSE_Nov14.csv)
#' @param leaf_area     File path to the leaf area file (LeafArea_monitoring_SMSE.csv)
#' @param axial_angle   File path to the leaf axial angle file (LeafDispositionComp_SMSE14.csv)
#' @param petiole_width File path to the leaf petiole width file (Petiole_SMSE14.csv)
#' @param twist         File path to the leaf twist file (Torsion_SMSE14.csv), used to compute
#'                      leaf nerve height.
#' @param map           Physiological age requested in month after planting
#'
#' @return A list of the imported and pre-processed files
#'
#' @importFrom dplyr ungroup group_by summarise select "%>%" mutate filter
#' @importFrom data.table fread
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' Vpalmr::import_data(parameter= 'Archi/ParameterSimu.csv',
#'                     development= 'Archi/Development_Rep4_SMSE.csv',
#'                     phylotaxy= "Archi/Stem_SMSE14.csv",
#'                     declination= "Archi/AnglesC&A_SMSE_Nov14.csv",
#'                     curvature= "Archi/LeafCurvature_SMSE14.csv",
#'                     torsion= 'Archi/AnglesC&A_SMSE_Nov14.csv',
#'                     leaf_area= 'Archi/LeafArea_monitoring_SMSE.csv',
#'                     axial_angle= "Archi/LeafDispositionComp_SMSE14.csv",
#'                     petiole_width= "Archi/Petiole_SMSE14.csv",
#'                     twist= "Archi/Torsion_SMSE14.csv")
#' }
#'
#' @export
#'
import_data= function(parameter,development,phylotaxy,declination,curvature,torsion,
                      leaf_area,axial_angle,petiole_width,twist,map){

  # Parameter ---------------------------------------------------------------

  Parameter=
    data.table::fread(parameter, data.table = F)%>%
    mutate(Date= as.Date(.data$Date,format='%d/%m/%Y'))

  nbLeafEmitted=
    Parameter%>%
    filter(.data$MAP==map)%>%
    group_by(.data$Progeny)%>%
    summarise(Physio_age= unique(.data$nbLeaves))

  if(nrow(nbLeafEmitted)<1){
    stop("map requested by the user not found in the data provided")
  }

  # DataAll -----------------------------------------------------------------

  DataAll=
    data.table::fread(development, dec=',', data.table= F)%>%
    mutate(NurseryPlanting_Date=as.Date(.data$NurseryPlanting_Date,format='%d/%m/%Y'),
           Transplanting_Date=as.Date(.data$Transplanting_Date,format='%d/%m/%Y'),
           Observation_Date=as.Date(.data$Observation_Date,format='%d/%m/%Y'),
           #____conversion in cm
           Leaflet_length1=.data$Leaflet_length1/10,
           Leaflet_length2=.data$Leaflet_length2/10,
           Leaflet_length3= .data$Leaflet_length3/10,
           Leaflet_length4= .data$Leaflet_length4/10,
           Leaflet_width1= .data$Leaflet_width1/10,
           Leaflet_width2= .data$Leaflet_width2/10,
           Leaflet_width3= .data$Leaflet_width3/10,
           Leaflet_width4= .data$Leaflet_width4/10,
           # StemDiameter------stem basis diameter
           StemDiameter= .data$StemCircumference/pi)

  DataAll%<>%
    merge(nbLeafEmitted, by="Progeny")

  # Nb_frond: nb of leaves emitted between 2 dates
  # + TotalEmitted: cumulative nb of leaves emitted from planting date
  # Some observations for Nb_frond are missing, computing it from LeafIndexRank1

  DataAll%<>%
    group_by(.data$TreeNumber,.data$MonthAfterPlanting)%>%
    summarise(Nb_frond_new= mean(.data$Nb_frond), LeafIndexRank1= max(.data$LeafIndexRank1))%>%
    mutate(Nb_frond_new= ifelse(is.na(.data$Nb_frond_new),
                                .data$LeafIndexRank1-lag(.data$LeafIndexRank1),
                                .data$Nb_frond_new),
           TotalEmitted= cumsum(.data$Nb_frond_new))%>%
    select(-.data$LeafIndexRank1)%>%
    merge(DataAll,., by= c("TreeNumber","MonthAfterPlanting"),sort = F)%>%
    mutate(Nb_frond= .data$Nb_frond_new)%>%select(-.data$Nb_frond_new)%>%
    mutate(RatioPetiole= .data$PetioleLength/.data$RachisLength,
           LeafNumber= .data$TotalEmitted - .data$FrondRank)

  #LastEmitted----------------last leaf emitted (= TotalEmitted at rank 1)
  DataAll%<>%
    group_by(.data$TreeNumber)%>%
    mutate(LastEmitted= max(.data$TotalEmitted,na.rm=T))%>%
    ungroup()

  #PosB----relative position of B point on rachis
  DataAll$PosB= DataAll$Bposition/DataAll$RachisLength

  #LeafletBLength----leaflet length at B point

  DataAll$LeafletBLength[is.na(DataAll$LeafletBLength)]=
    DataAll%>%
    filter(is.na(.data$LeafletBLength))%>%
    select(.data$Leaflet_length1,.data$Leaflet_length2,
           .data$Leaflet_length3,.data$Leaflet_length4)%>%
    rowMeans(.,na.rm=T)

  #LeafletBWidth----leaflet max width at B point, estimated from the mean
  # leaflets width from the 4 leaflets sample

  DataAll%<>%
    mutate(max_width= pmax(.data$Leaflet_width1,.data$Leaflet_width2,
                           .data$Leaflet_width3,.data$Leaflet_width4,na.rm=T))%>%
    mutate(LeafletBWidth= ifelse(is.na(.data$LeafletBWidth),.data$max_width,
                                 .data$LeafletBWidth))



  test_dev= test_development(x = DataAll)
  if(!is.null(test_dev)){
    warning(paste("Potential error in Development file:\n",test_dev))
  }

  # Phylotaxy data ----------------------------------------------------------

  Phylo= data.table::fread(phylotaxy, dec=',', data.table= F)
  Phylo$Phylo=(3*360+abs(Phylo$Azimuth_9-Phylo$Azimuth_17))/(17-9)


  # Declination at C point --------------------------------------------------

  Dec=
    data.table::fread(declination, dec= ',', data.table= F)%>%
    mutate(Decli_C= .data$BendingC+90)



  # Leaf curvature ----------------------------------------------------------

  Curve=
    data.table::fread(curvature, dec= ",", data.table= F)%>%
    mutate(Obs_Date=as.Date(.data$Obs_Date,'%d/%m/%y'),
           Y_distance_soil_cm= .data$Y_distance_cm + .data$Height_O)

  # Changing the reference point to C point, and removing O point
  Curve%<>%
    group_by(.data$TreeNumber,.data$Rank)%>%
    mutate(Y_distance_cm= .data$Y_distance_cm-.data$Y_distance_cm[.data$Point=='C'],
           X_distance_cm= .data$X_distance_cm-.data$X_distance_cm[.data$Point=='C'],
           Z_distance_cm= .data$Z_distance_cm-.data$Z_distance_cm[.data$Point=='C'])%>%
    filter(.data$Point!='O' & .data$Complet!='Non' & .data$Manipe!='Yes')

  #rachis length estimation
  Curve%<>%
    ungroup()%>%
    group_by(.data$TreeNumber,.data$Rank)%>%
    mutate(PositionRachisEstim=
             object_lenght(.data$X_distance_cm,.data$Y_distance_cm,
                           method= "smooth.spline")$Distance%>%cumsum,
           rachisLength= max(.data$PositionRachisEstim),
           RelativePositionRachisEstim= .data$PositionRachisEstim/.data$rachisLength)
  # NB: The rachis lenght can also be found by finding the value of

  Curve%<>%
    ungroup()%>%
    group_by(.data$Progeny, .data$TreeNumber, .data$Rank)%>%
    summarise(angC= find_angle(.data$X_distance_cm,.data$Y_distance_cm)[1],
              angA= find_angle(.data$X_distance_cm,.data$Y_distance_cm)[2],
              rachisLength= mean(.data$rachisLength))%>%
    select(-.data$rachisLength)%>%
    merge(Curve,.,by = c('Progeny','TreeNumber','Rank'),all.x = T, sort = F)


  # Leaf torsion -----------------------------------------------------------

  #####     Rachis twist  #####
  Tor= data.table::fread(torsion, dec= ',', data.table= F)
  Tor$Observation_Date= as.Date(Tor$Observation_Date,format='%d/%m/%Y')

  ###absolute value of the twist
  Tor$TwistA=abs(Tor$TwistA)


  # Leaf Area ---------------------------------------------------------------

  Area=
    data.table::fread(leaf_area, dec=',', data.table= F)%>%
    mutate(NurseryPlantingDate= as.Date(.data$NurseryPlantingDate,format='%d/%m/%Y'),
           FieldPlantingDate =as.Date(.data$FieldPlantingDate,format='%d/%m/%Y'),
           Obs_Date =as.Date(.data$Obs_Date,format='%d/%m/%Y'))

  # Fixing typos on the dataset:
  Area$LeafletRankOnSection[Area$TreeNumber=="90_20"&Area$LeafIndex==6&Area$Section==1]= 18
  Area$NbLeaflets[Area$TreeNumber=="104_31"&Area$LeafIndex==-28&Area$Section==5]= 10
  Area$NbLeaflets[Area$TreeNumber=="106_29"&Area$LeafIndex==-9&Area$Section==10]= 13
  Area$NbLeaflets[Area$TreeNumber=="107_19"&Area$LeafIndex==-28&Area$Section==9]= 12

  #####     Leaflets position on rachis    #####
  #### calibration from the MAP corresponding to nbLeafEmitted

  # Min and max measured MAP in DataAll:
  MAP_meas=
    DataAll%>%
    group_by(.data$Progeny)%>%
    filter(.data$TotalEmitted<= .data$Physio_age+20 &
             .data$TotalEmitted>= .data$Physio_age-20 &
             .data$TreeNumber%in%unique(Area$TreeNumber))%>%
    summarise(Map_Min= min(.data$MonthAfterPlanting),
              Map_Max= max(.data$MonthAfterPlanting))%>%
    mutate(Map_Min= ifelse(.data$Map_Min<min(Area$MAP),min(Area$MAP),.data$Map_Min),
           Map_Max= ifelse(.data$Map_Max>min(Area$MAP),max(Area$MAP),.data$Map_Max))
  # NB: if min in DataAll is lower than the one in Area, we take the min from Area
  # as a proxy, and same for the maximum.

  # Filtering Area only for rows with the desired MAP:
  Area%<>%
    merge(MAP_meas, by = "Progeny")%>%
    filter(.data$MAP<=.data$Map_Max & .data$MAP>= .data$Map_Min)%>%
    select(-.data$Map_Max,-.data$Map_Min)

  # Computing the total number of leaflets per leaf from the sum of the number of leaflets
  # that are on position 0
  Area=
    Area%>%
    group_by(.data$TreeNumber,.data$LeafIndex)%>%
    mutate(NbLeaflets_0= ifelse(.data$PositionOnLeaflet==0,.data$NbLeaflets,0),
           TotalLeaflets= sum(.data$NbLeaflets_0))%>%
    select(-.data$NbLeaflets_0)

  # Re-computing the leaflet rank on rachis:
  Area%<>%
    ungroup()%>%
    group_by(.data$Progeny,.data$TreeNumber,.data$LeafIndex,.data$Section)%>%
    arrange(.by_group = T)%>%
    summarise(LeafletRankOnSection= unique(.data$LeafletRankOnSection),
              NbLeaflets= unique(.data$NbLeaflets))%>%
    ungroup()%>%group_by(.data$Progeny,.data$TreeNumber,.data$LeafIndex)%>%
    mutate(cum_LeafletRank= cumsum(.data$NbLeaflets))%>%
    mutate(LeafletRank= lag(.data$cum_LeafletRank,default = 0)+
             .data$LeafletRankOnSection)%>%
    select(-.data$cum_LeafletRank,-.data$LeafletRankOnSection,-.data$NbLeaflets)%>%
    merge(Area,., by=c("Progeny","TreeNumber","LeafIndex","Section"),sort = F)

  # Deriving some new varibales:
  Area=
    Area%>%
    mutate(RelativeLeafletRank= .data$LeafletRank/.data$TotalLeaflets,
           RelativePositionRachis= .data$PositionOnRachis/.data$LeafLength)


  test_area= test_Area(x = Area)
  if(!is.null(test_area)){
    warning(paste("Potential error in Area file:\n",test_area))
  }
  # Leaf axial angle --------------------------------------------------------

  LftAngle=
    data.table::fread(axial_angle, dec=',', data.table= F)%>%
    mutate(Axial= AlKashi(.data$Axial_cm,32.5),
           Alpha= 90-AlKashi(.data$Radial_cm,32.5),
           #estimation of the projected radia angle:
           Radial_deg= atan(sin(.data$Axial*pi/180)*
                              tan(.data$Alpha*pi/180))*180/pi)%>%
    group_by(.data$TreeNumber)%>%
    mutate(Rachis_length= max(.data$Position, na.rm= T),
           Position_rel= .data$Position/.data$Rachis_length,
           Nb_leaflets= max(.data$Leaflet_rank, na.rm= T))



  # Leaf nerve width at C point ---------------------------------------------

  PetioleSectionC=
    data.table::fread(petiole_width, dec=',', data.table= F)%>%
    mutate(Obs_Date= as.Date(.data$Obs_Date,'%d/%m/%y'),
           Petiole_width_C_cm= .data$Petiole_width_C_mm/10,
           Petiole_height_C_cm= .data$Petiole_height_C_mm/10)%>%
    filter(!is.na(.data$Petiole_width_C_cm) & !is.na(.data$CtoA))



  # Nerve Height ------------------------------------------------------------

  RachisHeight=
    data.table::fread(twist,dec=',', data.table= F)%>%
    mutate(Obs_Date= as.Date(.data$Obs_Date,'%d/%m/%y'),
           Position_rachis_rel= .data$Position_rachis/.data$Rachis_length,
           Petiole_width_cm= .data$Petiole_width/10,
           Petiole_height_cm= .data$Petiole_height/10)%>%
    filter(!is.na(.data$Petiole_width_cm))%>%
    group_by(.data$TreeNumber)%>%
    mutate(Petiole_max_width= max(.data$Petiole_width, na.rm= T),
           Petiole_max_height= max(.data$Petiole_height, na.rm=T))%>%
    ungroup()%>%
    mutate(Petiole_relative_width= .data$Petiole_width/.data$Petiole_max_width,
           Petiole_relative_height= .data$Petiole_height/.data$Petiole_max_height)

  out= list(Parameter,DataAll,Phylo,Dec,Curve,Tor,Area,LftAngle,PetioleSectionC,
            RachisHeight,map)
  names(out)=
    c("Parameter","DataAll","Phylo","Dec","Curve","Tor","Area","LftAngle",
      "PetioleSectionC","RachisHeight","MAP_requested")

  return(out)
}
