#' Parameter initialization
#'
#' @description Return the default initialisation for the parameters
#'
#' @param seed The seed to use for random generation
#'
#' @return The default parameter list
#'
#' @examples
#'\dontrun{
#' init_list(10)
#'}
#' @export
init_list= function(){
  list(
    trunkBending_M=0,
    trunkBending_SD=0,
    # internodeFinalLength=2,
    laminaAngle = 140,
    leafletStiffness= 1500,    # leafletStifness= 0     --> soft (pending) leaflets
                               # leafletStifness= 10000 --> erect leaflets
    leafletStiffness_SD= 7000, # nerve shape section at leaf extremities
    frondBaseWidth = 30,
    frondtipWidth = 0.3,
    frondBaseHeight = 10,
    frondTipHeight = 0.5,
    nbInflorescences = 0,      # not yet implemented
    rachisTwistCoef=3,
    rachisDevCoef=2
  )
}
