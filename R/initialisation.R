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
init_list= function(seed= 1){
  list(
    seed= seed,
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
    rachisDevCoef=2,
    #number of iterations for parameters estimation:
    param_controle= nlme::lmeControl(maxIter=500000, niterEM=25000),
    epsilon=10^-6,             # epsilon value for matrix inversion
    nbFronds_M= 45,            # number of fronds in 3D mock ups
    nbTreeSimu=25,             # number of trees simulated
    side=1                     # leaf side
  )
}
