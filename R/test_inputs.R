#' Test the input data
#'
#' @description Performs several tests on the VPalm input files to check for potential
#'  errors such as wrong format, missing essential data or typos.
#'
#' @param x The input data.frame
#'
#' @return `NULL` if no major issue was found in the file, or a character string of all
#' cumulated potential errors in the data.frame.
#'
#' @details The outputs of the functions are used to make global warnings, messages or errors
#'
#' @importFrom dplyr group_by "%>%" arrange summarise filter mutate
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' Area= data.table::fread('LeafArea_monitoring_SMSE.csv',data.table = F,dec=',')
#' test_Area(Area= Area)
#'}
#' @export
test_Area= function(x){
  # Testing that the number of leaflets and their rank are the same in each
  # section of the leaf:
  warn= NULL

  Test=
    x%>%
    group_by(TreeNumber,LeafIndex,Section)%>%
    arrange(.by_group = T)%>%
    summarise(LeafletRankOnSection= length(unique(LeafletRankOnSection)),
              NbLeaflets= length(unique(NbLeaflets)))%>%
    filter(NbLeaflets>1|LeafletRankOnSection>1)%>%
    mutate(LeafletRankOnSection= ifelse(LeafletRankOnSection>1,"Problem","OK"),
           NbLeaflets= ifelse(NbLeaflets>1,"Problem","OK"))

  if(nrow(Test)>0){
    warn_inc(paste("* Problem on input Area data: wrong value. Hint:",
                   "Please check that the number of leaflets and their rank",
                   "are the same in each section of the leaves mentionned"))
  }

  # Test for potential error on maximum length:
  Lf=
    x%>%
    filter(Width==0 & PositionOnLeaflet!=0)%>%
    mutate(Position_rachis_rel= PositionOnRachis/LeafLength)%>%
    group_by(Progeny,TreeNumber,LeafIndex)%>%
    mutate(Max_length= max(PositionOnLeaflet))%>%
    ungroup()%>%
    mutate(Relative_length= PositionOnLeaflet/Max_length,
           leafID= paste(TreeNumber,'Leaf', LeafIndex))

  potential_error=
    Lf[Lf$Relative_length==1 &
         Lf$Position_rachis_rel<0.2 |
         Lf$Relative_length==1 &
         Lf$Position_rachis_rel>0.8,]

  if(nrow(potential_error)>0){
    warn= warn_inc(warn,paste('* Probable error in max length for Tree:',
                              paste(potential_error$leafID, 'section',
                                    potential_error$Section,"\n")))
  }


  # Test for potential error in leaflets width:
  WidthAll=
    x%>%
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
           leafID= paste(TreeNumber,'Leaf',LeafIndex))

  potential_error=
    WidthAll[WidthAll$Relative_max_width ==1 &
               WidthAll$Position_rachis_rel<0.3 |
               WidthAll$Relative_max_width ==1 &
               WidthAll$Position_rachis_rel>0.8,]

  if(nrow(potential_error)>0){
    warn= warn_inc(warn,paste('* Probable error in width for Tree:',
                              paste(potential_error$leafID, 'section',
                                    potential_error$Section,"\n")))
  }

  return(warn)
}

#' @rdname test_Area
#' @export
test_development= function(x){
  warn= NULL

  Pet=
    x%>%
    filter(TotalEmitted <= Physio_age+60 &
             TotalEmitted >= Physio_age-60 &
             FrondRank > 17 &
             RatioPetiole > 0.5)%>%
    group_by(Progeny,TreeNumber,LeafIndex)%>%
    summarise(Pot_error= n())

  if(nrow(Pet)>0){
    ID= paste0('Tree',Pet$TreeNumber,', leaf', Pet$LeafIndex,"\n")
    warn= warn_inc(warn,
                   paste('* Probable error in petiole length or rachis length for:\n',
                         paste(ID)))
  }
  return(warn)
}


