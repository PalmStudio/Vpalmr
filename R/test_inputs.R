#' Test the Area data.frame
#'
#' @description Performs several tests on the Area file to check for potential errors or typos
#'
#' @param Area      The Area data.frame
#'
#' @return `TRUE` if the file had no major problem, or a special output that points to the potential
#' error
#' @importFrom dplyr group_by "%>%" arrange summarise filter mutate
#'
#' @examples
#'\dontrun{
#' library(Vpalmr)
#' Area= data.table::fread('LeafArea_monitoring_SMSE.csv',data.table = F,dec=',')
#' test_Area(Area= Area)
#'}
#' @export
test_Area= function(Area){
  # Testing that the number of leaflets and their rank are the same in each
  # section of the leaf:
  Test=
    Area%>%
    group_by(TreeNumber,LeafIndex,Section)%>%
    arrange(.by_group = T)%>%
    summarise(LeafletRankOnSection= length(unique(LeafletRankOnSection)),
              NbLeaflets= length(unique(NbLeaflets)))%>%
    filter(NbLeaflets>1|LeafletRankOnSection>1)%>%
    mutate(LeafletRankOnSection= ifelse(LeafletRankOnSection>1,"Problem","OK"),
           NbLeaflets= ifelse(NbLeaflets>1,"Problem","OK"))

  if(nrow(Test)>0){
    message("Problem on input Area data: wrong value")
    message(paste("Please check that the number of leaflets and their rank",
                  "are the same in each section of the leaves mentionned"))
    return(Test)
  }else{
    return(TRUE)
  }
}

#' @rdname test_Area
#' @export
test_petiole_width= function(x){
  if(nrow(na.omit(x[x$RatioPetiole>0.5,]))>0){
    ID=paste('tree ',x$TreeNumber,' leaf ', x$LeafIndex,sep='')
    message('Probable error in petiole length or rachis length for:\n',
            paste(ID[x$RatioPetiole>0.5]))
  }
}
