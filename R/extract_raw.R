#' Extract section
#'
#' @description Extract a section from the form
#'
#' @param Row The row index of the section
#' @param Col The column index of the section
#'
#' @details The row and column of the section is the one from the upper left corner
#'
#' @return A data.frame with the section form data
#' @export
#'
#' @examples
#' \dontrun{
#' extract_section_area(5,1)
#' }
extract_section_area= function(df_archi,Row,Col){
  data.frame(
    NbLeaflets= as.numeric(df_archi[Row+1,Col+1]),
    LeafletRankOnSection= as.numeric(df_archi[Row+2,Col+1]),
    PositionOnLeaflet= as.numeric(df_archi[(Row+2):(Row+8),Col+3]),
    Width= as.numeric(df_archi[(Row+2):(Row+8),Col+4]),
    Rachis_section_FW= as.numeric(df_archi[Row+3,Col+1]),
    Rachis_sample_FW= as.numeric(df_archi[Row+4,Col+1]),
    Rachis_sample_DW= as.numeric(df_archi[Row+5,Col+1]),
    Leaflet_sample_FW= as.numeric(df_archi[Row+6,Col+1]),
    Leaflet_sample_DW= as.numeric(df_archi[Row+7,Col+1]),
    Leaflet_section_FW= as.numeric(df_archi[Row+8,Col+1]),
    Leaflet_section_DW= as.numeric(df_archi[Row+9,Col+1]),
    Petiole_section_Width= as.numeric(df_archi[Row+9,Col+3]),
    Petiole_section_Thickness= as.numeric(df_archi[Row+9,Col+4])
  )
}




#' Extract Area form
#'
#' @description Extract the data from a form of one leaf from one Tree. This function is called
#' by [extract_sheets()]. Users should use [extract_sheets()] instead of this function.
#'
#' @param form The form (usually a sheet from an excel file)
#' @return A data.frame with the different measurements
#' @export
#'
#' @examples
#' \dontrun{
#' extract_form_area(form)
#' }
#'
extract_form_area= function(form){
  form= apply(form, 2, as.character)

  # Position of the sections in the form:
  sections_positions= data.frame(Rows= c(rep(seq(5,45,10),2)), Cols= c(rep(1,5),rep(7,5)))

  # Extract each section:
  sections_list=
    lapply(1:nrow(sections_positions), function(x){
      extract_section_area(df_archi = form, Row = sections_positions[x,1],Col = sections_positions[x,2])%>%
        dplyr::mutate(Section= x)
    })
  sections_df= data.table::rbindlist(sections_list)%>%tibble::as_tibble()

  Obs_once=
    data.frame(obs_date= dmy(form[1,2]),
               Progeny= form[2,2],
               TreeNumber= form[3,2],
               LeafIndex= as.numeric(form[1,5]),
               LeafIndexRank1= as.numeric(form[2,5]),
               RachisLength= as.numeric(form[1,8]))

  Obs_once[rep(1,nrow(sections_df)),]%>%
    cbind(sections_df)%>%
    dplyr::mutate(PositionOnRachis= .data$RachisLength*(2.0*.data$Section-1)/20)
}


#' Extract development form
#'
#' @description Extract the data from a leaf development form. This function is called
#' by [extract_sheets()]. Users should use [extract_sheets()] instead of this function.
#'
#' @param form The data form on which extraction is made
#'
#' @return A data.frame with the data from all sheets (all leaves)
#' @export
#'
#' @examples
#' \dontrun{
#' extract_form_dev(form)
#' }
extract_form_dev= function(form){
  nbRow=nrow(form)
  extr=data.frame(
    Observation_Date=as.Date(as.numeric(as.character(form[1,2])),origin='1899-12-30'),
    Progeny=as.character(form[1,4]),
    TreeNumber=as.character(form[3:nbRow,1]),
    LeafIndexRank1=as.numeric(as.character((form[3:nbRow,2]))),
    StemHeight17=as.numeric(as.character((form[3:nbRow,3]))),
    LeafIndex=as.numeric(as.character((form[3:nbRow,4]))),
    PetioleLength=as.numeric(as.character((form[3:nbRow,5]))),
    Bposition=as.numeric(as.character((form[3:nbRow,6]))),
    RachisLength=as.numeric(as.character((form[3:nbRow,7]))),
    LeafletBLength=as.numeric(as.character((form[3:nbRow,8]))),
    LeafletBWidth=as.numeric(as.character((form[3:nbRow,9])))
  )
  return(extr)
}


#' Extract sheets
#'
#' @description Extract the data from all forms and all sheets from the leaf area/biomass or
#' the development excel files
#'
#' @param path The path to the excel file
#' @param form The form type (see form section)
#' @param sheet The sheets index on the files (see details)
#'
#' @details If the sheet is not given, the function will import data from all sheets. Otherwise,
#' the user can provide a given sheet, or a vector of sheets (see examples).
#'
#' @section form
#' The form argument is either 'area' or 'development', and correspond to the type of form that is
#' needed: either the ones from leaf area and biomass, or the ones from leaf development.
#'
#' @note If you have an error as the following one: "JAVA_HOME cannot be determined from the Registry",
#' please check that you have a JAVA version that match you R one, *i.e.* that your JAVA is 64-bits if
#' your R is 64-bits also. If JAVA is 64-bits, it should be installed in "C:/Program Files/Java",
#' otherwise it is located in "C:/Program Files (x86)/Java"
#'
#' @return A data.frame with the data from all sheets (all leaves)
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract all sheets:
#' extract_sheets(path= "1-Data/raw/Form Archi.xlsx",form='area')
#'
#' # Extract only the first sheet:
#' form(path= "1-Data/raw/Form Archi.xlsx", form='area', sheet= 1)
#'
#' # Extract the second and third sheets:
#' extract_sheets(path= "1-Data/raw/Form Archi.xlsx",form='area', sheet= c(2,3))
#' }
extract_sheets= function(path,form=c('area','development'),sheet=NULL){
  form=match.arg(form,c('area','development'))
  if(is.null(sheet)){
    N_sheets= 1:(xlsx::getSheets(xlsx::loadWorkbook(path))%>%length())
  }
  lapply(N_sheets, function(x){
    tryCatch(expr = {
      df= xlsx::read.xlsx(path,sheetIndex = x, header = F, colClasses = "character")
      if(form=='area'){
        extract_form_area(form = df)
      }else if(form=='development') {
        extract_form_dev(form = df)
      }
    },
    error=function(cond) {
      message(paste("Error during sheet extraction for sheet:",x))
      message("Here's the original error message:")
      message(cond)
      return(NULL)
    })
  })%>%data.table::rbindlist()
}

