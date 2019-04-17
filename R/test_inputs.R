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
    dplyr::group_by(.data$TreeNumber, .data$LeafIndex, .data$Section)%>%
    dplyr::arrange(.by_group = T)%>%
    dplyr::summarise(LeafletRankOnSection= length(unique(.data$LeafletRankOnSection)),
                     NbLeaflets= length(unique(.data$NbLeaflets)))%>%
    dplyr::filter(.data$NbLeaflets>1|.data$LeafletRankOnSection>1)%>%
    dplyr::mutate(LeafletRankOnSection= ifelse(.data$LeafletRankOnSection>1,"Problem","OK"),
                  NbLeaflets= ifelse(.data$NbLeaflets>1,"Problem","OK"))

  if(nrow(Test)>0){
    warn_inc(warn,paste("* Problem on input Area data: wrong value. Hint:",
                        "Please check that the number of leaflets and their rank",
                        "are the same in each section of the leaves mentionned"))
  }

  # Test for potential error on maximum length:
  Lf=
    x%>%
    dplyr::filter(.data$Width==0 & .data$PositionOnLeaflet!=0)%>%
    dplyr::mutate(Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength)%>%
    dplyr::group_by(.data$Progeny,.data$TreeNumber,.data$LeafIndex)%>%
    dplyr::mutate(Max_length= max(.data$PositionOnLeaflet, na.rm= T))%>%
    dplyr::ungroup()%>%
    dplyr::mutate(Relative_length= .data$PositionOnLeaflet/.data$Max_length,
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
    dplyr::filter(.data$PositionOnLeaflet!=0)%>%
    dplyr::group_by(.data$Progeny, .data$TreeNumber, .data$LeafIndex, .data$Section)%>%
    dplyr::mutate(Leaflet_length= max(.data$PositionOnLeaflet, na.rm= T),
                  Leaflet_max_width= max(.data$Width, na.rm= T),
                  Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength,
                  Width_rel= .data$Width/.data$Leaflet_max_width)%>%
    dplyr::filter(.data$Width==0)%>%ungroup()%>%
    dplyr::group_by(.data$TreeNumber,.data$LeafIndex)%>%
    dplyr::mutate(Max_max_width= max(.data$Leaflet_max_width, na.rm= T),
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
      dplyr::ungroup()%>%
      dplyr::group_by(.data$Progeny,.data$TreeNumber,.data$MAP,.data$LeafIndex,.data$Section)%>%
      dplyr::arrange(.by_group = T)%>%
      dplyr::summarise(LeafletRankOnSection= unique(.data$LeafletRankOnSection))},
    error= function(cond){
      err=
        x%>%
        dplyr::ungroup()%>%
        dplyr::group_by(.data$Progeny,.data$TreeNumber,.data$LeafIndex,.data$Section)%>%
        dplyr::arrange(.by_group = T)%>%
        dplyr::summarise(N_LeafletRankOnSection= length(unique(.data$LeafletRankOnSection)))%>%
        dplyr::filter(.data$N_LeafletRankOnSection>1)
      message(paste("Please check the LeafletRankOnSection for:"))
      message(paste(paste("Progeny",err$Progeny, "Tree", err$TreeNumber, "Leaf index",err$LeafIndex,
                          "Section", err$Section), collapse ="\n"))
      message("Here is the original message:")
      message(cond)
    })


  if(!is.null(path)){
    testplot=
      x%>%
      dplyr::group_by(.data$TreeNumber,.data$Section)%>%
      dplyr::mutate(Leaflet_length= max(.data$PositionOnLeaflet),
                    Leaflet_max_width= max(.data$Width),
                    Position_rachis_rel= .data$PositionOnRachis/.data$LeafLength,
                    Position_leaflet_rel= .data$PositionOnLeaflet/.data$Leaflet_length,
                    Width_rel= .data$Width/.data$Leaflet_max_width,
                    PositionRelative= .data$Section/10)%>%
      dplyr::mutate(info= paste("MAP", .data$MAP,"section", .data$Section))%>%
      ggplot2::ggplot(aes(x= .data$PositionOnLeaflet, y= .data$Width))+
      ggplot2::facet_wrap(.~.data$TreeNumber+.data$MAP)+
      ggplot2::geom_line(aes(color= as.factor(.data$Section)))+
      ggplot2::labs(color = 'Section')+
      ggplot2::ggtitle("Each line should be smooth and follow sort of a bell curve")

    ggplot2::ggsave(filename = file.path(path,"Leaflet_position_vs_width.png"), plot= testplot,
                    width = 60, height = 80, units = "cm",dpi = 120)

  }


  # Testing if the number of leaflets per section is unique:
  tryCatch(expr = {
    x%>%
      group_by(.data$TreeNumber,.data$MAP,.data$LeafIndex,.data$Section)%>%
      summarise(NbLeaflets_section= unique(.data$NbLeaflets))
  },
  error= function(cond){
    message("NbLeaflets should be unique for each section")
  })


  return(warn)
}

#' @rdname test_Area
#' @export
test_development= function(x, path= NULL){
  warn= NULL

  Pet=
    x%>%
    dplyr::filter(.data$TotalEmitted <= .data$Physio_age+60 &
                    .data$TotalEmitted >= .data$Physio_age-60 &
                    .data$FrondRank > 17 &
                    .data$RatioPetiole > 0.5)%>%
    dplyr::group_by(.data$Progeny, .data$TreeNumber, .data$LeafIndex)%>%
    dplyr::summarise(Pot_error= n())

  if(nrow(Pet)>0){
    ID= paste0('Tree',Pet$TreeNumber,', leaf', Pet$LeafIndex,"\n")
    warn= warn_inc(warn,
                   paste('* Probable error in petiole length or rachis length for:\n',
                         paste(ID)))
  }


  if(!is.null(path)){
    testplot=
      ggplot2::ggplot(data = x, aes(x = .data$MonthAfterPlanting, y= .data$TotalEmitted))+
      ggplot2::facet_wrap(.~.data$Progeny)+
      ggplot2::geom_line(aes(color= .data$TreeNumber))+
      ggplot2::labs(color = 'Tree')+
      ggplot2::ggtitle("Each line should be smooth and only be increasing")

    ggplot2::ggsave(filename = file.path(path,"MAP_vs_TotalEmitted.png"), plot= testplot,
                    width = 60, height = 80, units = "cm",dpi = 120)

  }

  # Testing that Petiole length is less than Rachis length
  df_err= x%>%filter(.data$PetioleLength>.data$RachisLength)

  if(nrow(df_err)>0){
    df_err%>%glue::glue_data("Error in RachisLength or PetioleLength for Tree {TreeNumber}, ",
                             "MAP {MonthAfterPlanting}, leaf index {LeafIndex}")%>%
      warn_inc(warn,.)
  }

  # Testing that PosB is not > 1 (it is a relative position)
  df_err= x%>%filter(.data$PosB>1.0&.data$PosB<0.2)

  if(nrow(df_err)>0){
    df_err%>%glue::glue_data("Error in RachisLength or Bposition for Tree {TreeNumber}, ",
                             "MAP {MonthAfterPlanting}, leaf index {LeafIndex}")%>%
      warn_inc(warn,.)
  }

  return(warn)
}


#' Test the position on leaflet
#'
#' @description Test the position on leaflet (`PositionOnLeaflet`) format and
#' values from one section of a leaf in a given palm tree.
#'
#' @param data A sample of the Area data file
#'
#' @return The input data.frame with one more column for the test: `err`. If the section
#' pass the test, `err` have the value "no error", else it will take a character value of
#' the given error.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Data frame with the right positions:
#' df_sub= data.frame(TreeNumber= "90_15", LeafIndex= -1, Section= 1,
#'                    PositionOnLeaflet= Vpalmr::position_on_leaflet(56))
#' # Testing (all should be ok):
#' test_pos_on_leaflet(df_sub)
#'
#' # Adding a wrong value:
#' df_sub[5,'PositionOnLeaflet']= 20
#' test_pos_on_leaflet(df_sub)
#' }
test_pos_on_leaflet= function(data){
  err=
    try({
      data%>%
        dplyr::select(.data$TreeNumber, .data$LeafIndex, .data$Section, .data$PositionOnLeaflet)%>%
        dplyr::mutate(PositionOnLeaflet_TRUE=
                        Vpalmr::position_on_leaflet(x = dplyr::last(.data$PositionOnLeaflet))
        )
      # NB: last(PositionOnLeaflet) gives the leaflet length
    }, silent = TRUE)

  # Tests:
  if(inherits(err, "try-error")){
    data$err= "Number of PositionOnLeaflet, should have 8 values"
  }else{
    if(!all(err$PositionOnLeaflet==err$PositionOnLeaflet_TRUE)){
      data$err= ifelse(err$PositionOnLeaflet==err$PositionOnLeaflet_TRUE,
                       "no error", "PositionOnLeaflet sequence is wrong")
      # data$err= "PositionOnLeaflet sequence is wrong"
    }else{
      data$err= "no error"
    }
  }
  data
}


#' Test all positions on leaflet
#'
#' @description Test the format and values of the positions on leaflets of
#' all sections of all leaves on all trees
#'
#' @param data The area data.frame
#'
#' @return A data.frame showing each line with potential error
#' @export
#'
#' @examples
#' \dontrun{
#' test_pos_on_leaflet_all(Area)
#' }
test_pos_on_leaflet_all= function(data){
  data%>%
    dplyr::filter(.data$MAP>47)%>% # The measurements were different before MAP 47
    dplyr::group_by(.data$MAP, .data$TreeNumber, .data$LeafIndex, .data$Section)%>%
    dplyr::arrange(.data$PositionOnLeaflet,.by_group = TRUE)%>%
    dplyr::do(test_pos_on_leaflet(.))%>%
    dplyr::summarise(err= unique(.data$err))%>%
    dplyr::filter(.data$err!="no error")
}
