#' Compute MAP for development data
#'
#' @description Compute the months after planting for the development data from
#' the transplanting date and the observation date. The function also corrects the
#' MAP if the field campaign is made over 2 consecutive months (should have the same MAP).
#'
#'
#' @param x The development data.frame
#'
#' @return The development data.frame with new or updated MAP.
#' @export
#'
#' @examples
#' \dontrun{
#' compute_MAP(development)
#' }
compute_MAP= function(x){

  # Lookup table for transplanting date for each progeny:
  Planting_date_df=
    x%>%
    select(.data$TreeNumber,.data$Transplanting_Date)%>%
    stats::na.omit()%>%
    group_by(.data$TreeNumber)%>%
    summarise(Transplanting_Date= unique(.data$Transplanting_Date))

  # (Re-)computing the MAP:
  dplyr::left_join(x%>%dplyr::select(-.data$Transplanting_Date),Planting_date_df)%>%
    dplyr::mutate(MAP_comp= lubridate::interval(.data$Transplanting_Date,.data$Observation_Date)%/%
                    months(1))%>%
    dplyr::group_by(.data$TreeNumber)%>%
    dplyr::arrange(.data$Observation_Date)%>%
    dplyr::mutate(MAP_comp= ifelse(!is.na(.data$Nb_frond)&!is.na(lag(.data$Nb_frond))&
                                     .data$Nb_frond==lag(.data$Nb_frond)&
                                     !is.na(.data$LeafIndexRank1)&!is.na(lag(.data$LeafIndexRank1))&
                                     .data$LeafIndexRank1==lag(.data$LeafIndexRank1),
                                   lag(.data$MAP_comp),.data$MAP_comp))%>%
    dplyr::mutate(MAP_comp= ifelse(!is.na(.data$Nb_frond)&!is.na(lead(.data$Nb_frond))&
                                     .data$Nb_frond==lead(.data$Nb_frond)&
                                     !is.na(.data$LeafIndexRank1)&!is.na(lead(.data$LeafIndexRank1))&
                                     .data$LeafIndexRank1==lead(.data$LeafIndexRank1),
                                   lead(.data$MAP_comp),.data$MAP_comp))%>%ungroup()%>%
    dplyr::mutate(MonthAfterPlanting= .data$MAP_comp)%>%
    dplyr::select(-.data$MAP_comp)
  # NB: The last two mutates are used for the case when one session is made on several different days
}
