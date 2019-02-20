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
#' extract_section(5,1)
#' }
extract_section= function(df_archi,Row,Col){
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




#' Extract form
#'
#' @description Extract the data from a form of one leaf from one Tree
#' @param form The form (usually a sheet from an excel file)
#' @return A data.frame with the different measurements
#' @export
#'
#' @examples
#' \dontrun{
#' extract_form(form)
#' }
#'
extract_form= function(form){
  form= apply(form, 2, as.character)

  # Position of the sections in the form:
  sections_positions= data.frame(Rows= c(rep(seq(5,45,10),2)), Cols= c(rep(1,5),rep(7,5)))

  # Extract each section:
  sections_list=
    lapply(1:nrow(sections_positions), function(x){
      extract_section(df_archi = form, Row = sections_positions[x,1],Col = sections_positions[x,2])%>%
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



#' Extract sheets
#'
#' @description Extract the data from all forms and all sheets from the leaf area and
#'  biomass excel files
#'
#' @param path The path to the excel file
#' @param sheet The sheets index on the files (see details)
#'
#' @details If the sheet is not given, the function will import data from all sheets. Otherwise,
#' the user can provide a given sheet, or a vector of sheets (see examples).
#'
#' @note If you have an error as the following one: "JAVA_HOME cannot be determined from the Registry",
#' please chack that you have a JAVA version that match you R one, *i.e.* that your JAVA is 64-bits if
#' your R is 64-bits also. If JAVA is 64-bits, it should be installed in "C:/Program Files/Java",
#' otherwise it is located in "C:/Program Files (x86)/Java"
#'
#' @return A data.frame with the data from all sheets (all leaves)
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract all sheets:
#' extract_sheets(path= "1-Data/raw/Form Archi.xlsx")
#'
#' # Extract only the first sheet:
#' extract_sheets(path= "1-Data/raw/Form Archi.xlsx", sheet= 1)
#'
#' # Extract the second and third sheets:
#' extract_sheets(path= "1-Data/raw/Form Archi.xlsx", sheet= c(2,3))
#' }
extract_sheets= function(path, sheet=NULL){
  if(is.null(sheet)){
    N_sheets= 1:(xlsx::getSheets(xlsx::loadWorkbook(path))%>%length())
  }
  lapply(N_sheets, function(x){
    tryCatch(expr = {
      df_archi= xlsx::read.xlsx(path,sheetIndex = x, header = F, colClasses = "character")
      extract_form(form = df_archi)
    },
    error=function(cond) {
      message(paste("Error during sheet extraction for sheet:",x))
      message("Here's the original error message:")
      message(cond)
      return(NULL)
    })
  })%>%data.table::rbindlist()
}
