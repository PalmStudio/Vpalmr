#' Test the input data
#'
#' @description Performs several tests on the VPalm input files to check for potential
#'  errors such as wrong format, missing essential data or typos.
#'
#' @param x The input data.frame
#' @param path The path to a directory were to write the plots
#'
#' @return `NULL` if no major issue was found in the file, or a character string of all
#' cumulated potential errors in the data.frame.
#'
#' @details The outputs of the functions are used to make global warnings, messages or errors
#'
#' @importFrom dplyr group_by "%>%" arrange summarise filter mutate n
#' @importFrom rlang .data
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' Area= data.table::fread('LeafArea_monitoring_SMSE.csv',data.table = F,dec=',')
#' test_Area(Area= Area)
#'}
#' @export
test_Area= function(x, path= NULL){
  # Testing that the number of leaflets and their rank are the same in each
  # section of the leaf:
  warn= NULL

  Test=
    x%>%
    group_by(.data$TreeNumber, .data$LeafIndex, .data$Section)%>%
    arrange(.by_group = T)%>%
    summarise(LeafletRankOnSection= length(unique(.data$LeafletRankOnSection)),
              NbLeaflets= length(unique(.data$NbLeaflets)))%>%
    filter(.data$NbLeaflets>1|.data$LeafletRankOnSection>1)%>%
    mutate(LeafletRankOnSection= ifelse(.data$LeafletRankOnSection>1,"Problem","OK"),
           NbLeaflets= ifelse(.data$NbLeaflets>1,"Problem","OK"))

  if(nrow(Test)>0){
    warn_inc(warn,paste("* Problem on input Area data: wrong value. Hint:",
                   "Please check that the number of leaflets and their rank",
                   "are the same in each section of the leaves mentionned"))
  }

  # Test for potential error on maximum length:
  Lf=
    x%>%
    filter(.data$Width==0 & .data$PositionOnLeaflet!=0)%>%
    mutate(Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength)%>%
    group_by(.data$Progeny,.data$TreeNumber,.data$LeafIndex)%>%
    mutate(Max_length= max(.data$PositionOnLeaflet, na.rm= T))%>%
    ungroup()%>%
    mutate(Relative_length= .data$PositionOnLeaflet/.data$Max_length,
           leafID= paste(.data$TreeNumber,'Leaf', .data$LeafIndex))

  potential_error=
    Lf[Lf$Relative_length==1 &
         Lf$Position_rachis_rel<0.2 |
         Lf$Relative_length==1 &
         Lf$Position_rachis_rel>0.8,]

  if(nrow(potential_error)>0){
    warn= warn_inc(warn,paste('* Probable error in max length for Tree:',
                              paste(potential_error$leafID, 'section',
                                    potential_error$Section,
                                    ". The longest leaflet is at leaf extremity","\n")))
  }


  # Test for potential error in leaflets width:
  WidthAll=
    x%>%
    filter(.data$PositionOnLeaflet!=0)%>%
    group_by(.data$Progeny, .data$TreeNumber, .data$LeafIndex, .data$Section)%>%
    mutate(Leaflet_length= max(.data$PositionOnLeaflet, na.rm= T),
           Leaflet_max_width= max(.data$Width, na.rm= T),
           Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength,
           Width_rel= .data$Width/.data$Leaflet_max_width)%>%
    filter(.data$Width==0)%>%ungroup()%>%
    group_by(.data$TreeNumber,.data$LeafIndex)%>%
    mutate(Max_max_width= max(.data$Leaflet_max_width, na.rm= T),
           Relative_max_width= .data$Leaflet_max_width/.data$Max_max_width,
           leafID= paste(.data$TreeNumber,'Leaf',.data$LeafIndex))

  potential_error=
    WidthAll[WidthAll$Relative_max_width ==1 &
               WidthAll$Position_rachis_rel<0.3 |
               WidthAll$Relative_max_width ==1 &
               WidthAll$Position_rachis_rel>0.8,]

  if(nrow(potential_error)>0){
    warn= warn_inc(warn,paste('* Probable error in width for Tree:',
                              paste(potential_error$leafID, 'section',
                                    potential_error$Section,
                                    "The larger leaflet is at leaf extremity.","\n")))
  }


  tryCatch({
    x%>%
      ungroup()%>%
      group_by(.data$Progeny,.data$TreeNumber,.data$MAP,.data$LeafIndex,.data$Section)%>%
      arrange(.by_group = T)%>%
      summarise(LeafletRankOnSection= unique(.data$LeafletRankOnSection))},
    error= function(cond){
      err=
        x%>%
        ungroup()%>%
        group_by(.data$Progeny,.data$TreeNumber,.data$LeafIndex,.data$Section)%>%
        arrange(.by_group = T)%>%
        summarise(N_LeafletRankOnSection= length(unique(.data$LeafletRankOnSection)))%>%
        filter(N_LeafletRankOnSection>1)
      message(paste("Please check the LeafletRankOnSection for:"))
      message(paste(paste("Progeny",err$Progeny, "Tree", err$TreeNumber, "Leaf index",err$LeafIndex,
                          "Section", err$Section), collapse ="\n"))
      message("Here is the original message:")
      message(cond)
    })


  if(!is.null(path)){
    testplot=
      x%>%
      group_by(.data$TreeNumber,.data$Section)%>%
      mutate(Leaflet_length= max(.data$PositionOnLeaflet),
             Leaflet_max_width= max(.data$Width),
             Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength,
             Position_leaflet_rel= .data$PositionOnLeaflet/.data$Leaflet_length,
             Width_rel= .data$Width/.data$Leaflet_max_width,
             PositionRelative= .data$Section/10)%>%
      mutate(info= paste("MAP", MAP,"section", Section))%>%
      ggplot(aes(x= PositionOnLeaflet, y= Width))+
      facet_wrap(.~TreeNumber+MAP)+
      geom_line(aes(color= as.factor(Section)))+
      labs(color = 'Section')+
      ggtitle("Each line should be smooth and follow sort of a bell curve")

    ggsave(filename = file.path(path,"Leaflet_position_vs_width.png"), plot= testplot,
           width = 60, height = 80, units = "cm",dpi = 120)

  }

  return(warn)
}

#' @rdname test_Area
#' @export
test_development= function(x, path= NULL){
  warn= NULL

  Pet=
    x%>%
    filter(.data$TotalEmitted <= .data$Physio_age+60 &
             .data$TotalEmitted >= .data$Physio_age-60 &
             .data$FrondRank > 17 &
             .data$RatioPetiole > 0.5)%>%
    group_by(.data$Progeny, .data$TreeNumber, .data$LeafIndex)%>%
    summarise(Pot_error= n())

  if(nrow(Pet)>0){
    ID= paste0('Tree',Pet$TreeNumber,', leaf', Pet$LeafIndex,"\n")
    warn= warn_inc(warn,
                   paste('* Probable error in petiole length or rachis length for:\n',
                         paste(ID)))
  }


  if(!is.null(path)){
    testplot=
      ggplot(data = development, aes(x = MonthAfterPlanting, y= TotalEmitted))+
      facet_wrap(.~Progeny)+
      geom_line(aes(color= TreeNumber))+
      labs(color = 'Tree')+
      ggtitle("Each line should be smooth and only be increasing")

    ggsave(filename = file.path(path,"MAP_vs_TotalEmitted.png"), plot= testplot,
           width = 60, height = 80, units = "cm",dpi = 120)

  }

  return(warn)
}


