


#' Sigmoidal function
#'
#' @description Estimate a sigmoidal function from parameters.
#'
#' @param file_archi Path to the architectural file
#' @param file_la    Path to the leaf area file
#' @param file_phylo Path to the phylotaxy file
#' @param file_angles Path to the leaf declination angles file
#'
#' @return A sigmoid
#'
#' @importFrom dplyr ungroup group_by summarise "%>%"
#' @importFrom data.table fread
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' f.sigmo(X= 1/10,max= 3,slope= 1,infl=5)
#'}
#' @export
#'
import_data= function(file_archi,file_la,file_phylo,file_angles){

  # Architectural data ------------------------------------------------------

  DataAll= data.table::fread(file_archi,header=T,sep=";",dec=',',data.table = F)
  DataAll$NurseryPlanting_Date=as.Date(DataAll$NurseryPlanting_Date,format='%d/%m/%Y')
  DataAll$Transplanting_Date=as.Date(DataAll$Transplanting_Date,format='%d/%m/%Y')
  DataAll$Observation_Date=as.Date(DataAll$Observation_Date,format='%d/%m/%Y')

  #____conversion in cm
  DataAll$Leaflet_length1=DataAll$Leaflet_length1/10
  DataAll$Leaflet_length2=DataAll$Leaflet_length2/10
  DataAll$Leaflet_length3=DataAll$Leaflet_length3/10
  DataAll$Leaflet_length4=DataAll$Leaflet_length4/10

  DataAll$Leaflet_width1=DataAll$Leaflet_width1/10
  DataAll$Leaflet_width2=DataAll$Leaflet_width2/10
  DataAll$Leaflet_width3=DataAll$Leaflet_width3/10
  DataAll$Leaflet_width4=DataAll$Leaflet_width4/10

  # StemDiameter------stem basis diameter
  DataAll$StemDiameter= DataAll$StemCircumference/pi

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
    mutate(Nb_frond= Nb_frond_new)%>%select(-Nb_frond_new)

  DataAll$RatioPetiole= DataAll$PetioleLength/DataAll$RachisLength

  DataAll$LeafNumber= DataAll$TotalEmitted-DataAll$FrondRank

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

  # Leaf area data ----------------------------------------------------------

  Area= read.csv(file_la,header=T,sep=";",dec=',')

  Area$NurseryPlantingDate =as.Date(Area$NurseryPlantingDate,format='%d/%m/%Y')
  Area$FieldPlantingDate =as.Date(Area$FieldPlantingDate,format='%d/%m/%Y')
  Area$Obs_Date =as.Date(Area$Obs_Date,format='%d/%m/%Y')

  ### to know on which tree the area was measured
  DataAll$areaData=rep('No',nrow(DataAll))
  for (t in unique(Area$TreeNumber)){
    DataAll[DataAll$TreeNumber==t,]$areaData='Yes'
  }

  # Phylotaxy data ----------------------------------------------------------

  # Read phylotaxy data and merge it to the main table
  DataAll=
    data.table::fread(file_phylo,header=T,sep=";",dec=',',data.table = F)%>%
    mutate(Phylo=(3*360+abs(Phylo$Azimuth_9-Phylo$Azimuth_17))/(17-9))%>%
    group_by(TreeNumber)%>%
    summarise(phylotaxy= mean(Phylo))%>%
    merge(DataAll,., by= c("TreeNumber"),sort = F)


  # Leaf declination angles at remarkable points ----------------------------
  Dec= data.table::fread(file_angles,header=T,sep=";",dec=',',data.table = F)
  Dec$Decli_C=Dec$BendingC+90
  Dec=Dec[Dec$Progeny==progeny,]

  return(list(DataAll,Dec,Curve))
}
