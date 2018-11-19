#' Import all data
#'
#' @description Import and pre-compute the input data for VPalm
#'
#' @param parameter     File path to the parameter file (ParameterSimu.csv)
#' @param development   File path to the development file (Development_Rep4_SMSE.csv)
#' @param phylotaxy     File path to the phylotaxy file (Stem_SMSE14.csv)
#' @param declination   File path to the declination file (A_SMSE_Nov14.csv)
#' @param curvature     File path to the leaf curvature file (LeafCurvature_SMSE14.csv)
#' @param toricity      File path to the leaf toricity file (AnglesC&A_SMSE_Nov14.csv)
#' @param leaf_area     File path to the leaf area file (LeafArea_monitoring_SMSE.csv)
#' @param axial_angle   File path to the leaf axial angle file (LeafDispositionComp_SMSE14.csv)
#' @param petiole_width File path to the leaf petiole width file (Petiole_SMSE14.csv)
#' @param twist         File path to the leaf twist file (Torsion_SMSE14.csv), used to compute
#'                      leaf nerve height.
#' @param map           Physiological age requested in month after planting
#'
#' @return A list of the imported and pre-processed files
#'
#' @importFrom dplyr ungroup group_by summarise "%>%"
#' @importFrom data.table fread
#'
#' @examples
#' \dontrun{
#' Vpalmr::import_data(parameter= 'Archi/ParameterSimu.csv',
#'                     development= 'Archi/Development_Rep4_SMSE.csv',
#'                     phylotaxy= "Archi/Stem_SMSE14.csv",
#'                     declination= "Archi/AnglesC&A_SMSE_Nov14.csv",
#'                     curvature= "Archi/LeafCurvature_SMSE14.csv",
#'                     toricity= 'Archi/AnglesC&A_SMSE_Nov14.csv',
#'                     leaf_area= 'Archi/LeafArea_monitoring_SMSE.csv',
#'                     axial_angle= "Archi/LeafDispositionComp_SMSE14.csv",
#'                     petiole_width= "Archi/Petiole_SMSE14.csv",
#'                     twist= "Archi/Torsion_SMSE14.csv")
#' }
#'
#' @export
#'
import_data= function(parameter,development,phylotaxy,declination,curvature,toricity,
                      leaf_area,axial_angle,petiole_width,twist,map){
  # Parameter ---------------------------------------------------------------

  Parameter=
    data.table::fread(parameter, data.table = F)%>%
    mutate(Date= as.Date(Date,format='%d/%m/%Y'))

  nbLeafEmitted=
    Parameter%>%
    filter(MAP==map)%>%
    group_by(Progeny)%>%
    summarise(Physio_age= unique(nbLeaves))


  # DataAll -----------------------------------------------------------------

  DataAll=
    data.table::fread(development, dec=',', data.table= F)%>%
    mutate(NurseryPlanting_Date=as.Date(NurseryPlanting_Date,format='%d/%m/%Y'),
           Transplanting_Date=as.Date(Transplanting_Date,format='%d/%m/%Y'),
           Observation_Date=as.Date(Observation_Date,format='%d/%m/%Y'),
           #____conversion in cm
           Leaflet_length1=Leaflet_length1/10,
           Leaflet_length2=Leaflet_length2/10,
           Leaflet_length3=Leaflet_length3/10,
           Leaflet_length4=Leaflet_length4/10,
           Leaflet_width1=Leaflet_width1/10,
           Leaflet_width2=Leaflet_width2/10,
           Leaflet_width3=Leaflet_width3/10,
           Leaflet_width4=Leaflet_width4/10,
           # StemDiameter------stem basis diameter
           StemDiameter= StemCircumference/pi)

  DataAll%<>%
    merge(nbLeafEmitted, by="Progeny")

  # Nb_frond: nb of leaves emitted between 2 dates
  # + TotalEmitted: cumulative nb of leaves emitted from planting date
  # Some observations for Nb_frond are missing, computing it from LeafIndexRank1

  DataAll%<>%
    group_by(TreeNumber,MonthAfterPlanting)%>%
    summarise(Nb_frond_new= mean(Nb_frond), LeafIndexRank1= max(LeafIndexRank1))%>%
    mutate(Nb_frond_new= ifelse(is.na(Nb_frond_new),
                                LeafIndexRank1-lag(LeafIndexRank1),
                                Nb_frond_new),
           TotalEmitted= cumsum(Nb_frond_new))%>%
    select(-LeafIndexRank1)%>%
    merge(DataAll,., by= c("TreeNumber","MonthAfterPlanting"),sort = F)%>%
    mutate(Nb_frond= Nb_frond_new)%>%select(-Nb_frond_new)%>%
    mutate(RatioPetiole= PetioleLength/RachisLength,
           LeafNumber= TotalEmitted - FrondRank)

  #LastEmitted----------------last leaf emitted (= TotalEmitted at rank 1)
  DataAll%<>%
    group_by(TreeNumber)%>%
    mutate(LastEmitted= max(TotalEmitted,na.rm=T))%>%
    ungroup()

  #PosB----relative position of B point on rachis
  DataAll$PosB= DataAll$Bposition/DataAll$RachisLength

  #LeafletBLength----leaflet length at B point

  DataAll$LeafletBLength[is.na(DataAll$LeafletBLength)]=
    DataAll%>%
    filter(is.na(LeafletBLength))%>%
    select(Leaflet_length1,Leaflet_length2,
           Leaflet_length3,Leaflet_length4)%>%
    rowMeans(.,na.rm=T)

  #LeafletBWidth----leaflet max width at B point, estimated from the mean
  # leaflets width from the 4 leaflets sample

  DataAll%<>%
    mutate(max_width= pmax(Leaflet_width1,Leaflet_width2,
                           Leaflet_width3,Leaflet_width4,na.rm=T))%>%
    mutate(LeafletBWidth= ifelse(is.na(LeafletBWidth),max_width,LeafletBWidth))



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
    mutate(Decli_C= BendingC+90)



  # Leaf curvature ----------------------------------------------------------

  Curve=
    data.table::fread(curvature, dec= ",", data.table= F)%>%
    mutate(Obs_Date=as.Date(Obs_Date,'%d/%m/%y'),
           Y_distance_soil_cm= Y_distance_cm + Height_O)

  # Changing the reference point to C point, and removing O point
  Curve%<>%
    group_by(TreeNumber,Rank)%>%
    mutate(Y_distance_cm= Y_distance_cm-Y_distance_cm[Point=='C'],
           X_distance_cm= X_distance_cm-X_distance_cm[Point=='C'],
           Z_distance_cm= Z_distance_cm-Z_distance_cm[Point=='C'])%>%
    filter(Point!='O' & Complet!='Non' & Manipe!='Yes')

  #rachis length estimation
  Curve%<>%
    group_by(TreeNumber,Rank)%>%
    mutate(PositionRachisEstim=
             object_lenght(X_distance_cm,Y_distance_cm,
                           method= "smooth.spline")$Distance%>%cumsum,
           rachisLength= max(PositionRachisEstim),
           RelativePositionRachisEstim= PositionRachisEstim/rachisLength)
  # NB: The rachis lenght can also be found by finding the value of

  Curve%<>%
    group_by(Progeny, TreeNumber, Rank)%>%
    summarise(angC= find_angle(X_distance_cm,Y_distance_cm)[1],
              angA= find_angle(X_distance_cm,Y_distance_cm)[2],
              rachisLength= mean(rachisLength))%>%
    select(-rachisLength)%>%
    merge(Curve,.,by = c('Progeny','TreeNumber','Rank'),all.x = T, sort = F)


  # Leaf toricity -----------------------------------------------------------

  #####     Rachis twist  #####
  Tor= data.table::fread(toricity, dec= ',', data.table= F)
  Tor$Observation_Date= as.Date(Tor$Observation_Date,format='%d/%m/%Y')

  ###absolute value of the twist
  Tor$TwistA=abs(Tor$TwistA)


  # Leaf Area ---------------------------------------------------------------

  Area=
    data.table::fread(leaf_area, dec=',', data.table= F)%>%
    mutate(NurseryPlantingDate= as.Date(NurseryPlantingDate,format='%d/%m/%Y'),
           FieldPlantingDate =as.Date(FieldPlantingDate,format='%d/%m/%Y'),
           Obs_Date =as.Date(Obs_Date,format='%d/%m/%Y'))

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
    group_by(Progeny)%>%
    filter(TotalEmitted<= Physio_age+20 & TotalEmitted>= Physio_age-20 &
             TreeNumber%in%unique(Area$TreeNumber))%>%
    summarise(Map_Min= min(MonthAfterPlanting), Map_Max= max(MonthAfterPlanting))%>%
    mutate(Map_Min= ifelse(Map_Min<min(Area$MAP),min(Area$MAP),Map_Min),
           Map_Max= ifelse(Map_Max>min(Area$MAP),max(Area$MAP),Map_Max))
  # NB: if min in DataAll is lower than the one in Area, we take the min from Area
  # as a proxy, and same for the maximum.

  # Filtering Area only for rows with the desired MAP:
  Area%<>%
    merge(MAP_meas, by = "Progeny")%>%
    filter(MAP<=Map_Max & MAP>= Map_Min)%>%
    select(-Map_Max,-Map_Min)

  # Computing the total number of leaflets per leaf from the sum of the number of leaflets
  # that are on position 0
  Area=
    Area%>%
    group_by(TreeNumber,LeafIndex)%>%
    mutate(NbLeaflets_0= ifelse(PositionOnLeaflet==0,NbLeaflets,0),
           TotalLeaflets= sum(NbLeaflets_0))%>%
    select(-NbLeaflets_0)

  # Re-computing the leaflet rank on rachis:
  Area%<>%
    group_by(Progeny,TreeNumber,LeafIndex,Section)%>%
    arrange(.by_group = T)%>%
    summarise(LeafletRankOnSection= unique(LeafletRankOnSection),
              NbLeaflets= unique(NbLeaflets))%>%
    ungroup()%>%group_by(Progeny,TreeNumber,LeafIndex)%>%
    mutate(cum_LeafletRank= cumsum(NbLeaflets))%>%
    mutate(LeafletRank= lag(cum_LeafletRank,default = 0)+LeafletRankOnSection)%>%
    select(-cum_LeafletRank,-LeafletRankOnSection,-NbLeaflets)%>%
    merge(Area,., by=c("Progeny","TreeNumber","LeafIndex","Section"),sort = F)

  # Deriving some new varibales:
  Area=
    Area%>%
    mutate(RelativeLeafletRank= LeafletRank/TotalLeaflets,
           RelativePositionRachis= PositionOnRachis/LeafLength)


  test_area= test_Area(x = Area)
  if(!is.null(test_area)){
    warning(paste("Potential error in Area file:\n",test_area))
  }
  # Leaf axial angle --------------------------------------------------------

  LftAngle=
    data.table::fread(axial_angle, dec=',', data.table= F)%>%
    mutate(Axial= AlKashi(Axial_cm,32.5),
           Alpha= 90-AlKashi(Radial_cm,32.5),
           #estimation of the projected radia angle:
           Radial_deg= atan(sin(Axial*pi/180)*
                              tan(Alpha*pi/180))*180/pi)%>%
    group_by(TreeNumber)%>%
    mutate(Rachis_length= max(Position, na.rm= T),
           Position_rel= Position/Rachis_length,
           Nb_leaflets= max(Leaflet_rank, na.rm= T))



  # Leaf nerve width at C point ---------------------------------------------

  PetioleSectionC=
    data.table::fread(petiole_width, dec=',', data.table= F)%>%
    mutate(Obs_Date= as.Date(Obs_Date,'%d/%m/%y'),
           Petiole_width_C_cm= Petiole_width_C_mm/10,
           Petiole_height_C_cm= Petiole_height_C_mm/10)%>%
    filter(!is.na(Petiole_width_C_cm) & !is.na(CtoA))



  # Nerve Height ------------------------------------------------------------

  RachisHeight=
    data.table::fread(twist,dec=',', data.table= F)%>%
    mutate(Obs_Date= as.Date(Obs_Date,'%d/%m/%y'),
           Position_rachis_rel= Position_rachis/Rachis_length,
           Petiole_width_cm= Petiole_width/10,
           Petiole_height_cm= Petiole_height/10)%>%
    filter(!is.na(Petiole_width_cm))%>%
    group_by(TreeNumber)%>%
    mutate(Petiole_max_width= max(Petiole_width, na.rm= T),
           Petiole_max_height= max(Petiole_height, na.rm=T))%>%
    ungroup()%>%
    mutate(Petiole_relative_width= Petiole_width/Petiole_max_width,
           Petiole_relative_height= Petiole_height/Petiole_max_height)

  out= list(Parameter,DataAll,Phylo,Dec,Curve,Tor,Area,LftAngle,PetioleSectionC,
            RachisHeight)
  names(out)=
    c("Parameter","DataAll","Phylo","Dec","Curve","Tor","Area","LftAngle",
      "PetioleSectionC","RachisHeight")

  return(out)
}
