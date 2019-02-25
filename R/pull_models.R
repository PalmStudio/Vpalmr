#' Simulation coefficient
#'
#' @description Computes the simulation coefficients from a sampling in
#' the variance-covariance matrix for VPalm input.
#'
#' @param data The output from one of the [lme models][mod_stem_diameter()]
#' @param epsilon epsilon value for matrix inversion
#' @param type The type of sampling performed on the parameter distribution (see details)
#'
#' @details If `type = sample`, `data` must have a `MatG` column that represents the variance-covariance
#' matrix, and a `coef_mean` for the mean coefficient values. It `type = mean`, `data` must have
#' a `coef_mean` column only.
#' Epsilon is used to avoid non semi-positive matrices, and is added to the matrix
#' diagonal to make it closer to the identity matrix.
#' The sampling is made using a normal distribution: `rnorm(mean= 0,sd= 1)`. If the type
#' is equal to 'sample', the function makes a random sample in the distribution of the
#' parameter using the mean coefficient and a sampled standard deviation. If it is set to
#' 'mean', only the mean coefficient is used.
#'
#' @return The coefficients as a [tibble::tibble()]
#'
#' @importFrom rlang .data
#'
#' @export
coef_sample= function(data, epsilon, type= c('sample','mean')){

  type= match.arg(type,c('sample','mean'))

  if(type=="sample"){
    data$MatG=
      lapply(data$MatG, function(x){
        diag(x)= diag(x)+epsilon
        x
      })

    data$coef_sd=
      lapply(data$MatG, function(x){
        matrix(data= stats::rnorm(n=ncol(x)),ncol=ncol(x))%*%chol(x)
      })

    data$coef_simu=
      data%>%
      dplyr::do(coef_simu= .$coef_mean + .$coef_sd)%>%
      dplyr::pull(.data$coef_simu)
  }else{
    data$coef_simu=
      data%>%
      dplyr::do(coef_simu= .$coef_mean)%>%
      dplyr::pull(.data$coef_simu)
  }
  data
}

#' Pull (n)lme output
#'
#' @description Format lme or nlme outputs in a standard [tibble::tibble()] and computes
#' some usefull informations such as covariance matrices. Carefull, this function is
#' poorly designed for the moment and is only applicable on the (n)lme used in the context
#' of the package because only two or three factors (lme and nlme resp.) can be used.
#'
#' @param data The output from one of the [lme models][mod_stem_diameter()] or
#' [nlme models][mod_stem_diameter()]
#' @param epsilon epsilon value for matrix inversion
#' @param type The type of sampling performed on the parameter distribution (see
#'  [coef_sample()])
#'
#' @details Epsilon is used to avoid non semi-positive matrices, and is added to the matrix
#' diagonal to make it closer to the identity matrix
#'
#' @return A [tibble::tibble()] with model outputs
#'
#' @export
pull_lme= function(data, epsilon= 10^-6, type= c('sample','mean')){
  data%>%
    dplyr::mutate(intercept= nlme::fixed.effects(.data$mod)[1],
                  slope= nlme::fixed.effects(.data$mod)[2],
                  cov= list(stats::vcov(.data$mod)),
                  coef_mean= list(nlme::fixed.effects(.data$mod)),
                  sigma= summary(.data$mod)$sigma,
                  label= list(colnames(nlme::random.effects(.data$mod))),
                  SdG1= as.numeric(nlme::VarCorr(.data$mod)['(Intercept)','StdDev']),
                  SdG2= as.numeric(nlme::VarCorr(.data$mod)[2,'StdDev']),
                  corG= as.numeric(nlme::VarCorr(.data$mod)[2,'Corr']),
                  MatG= list(matrix(
                    data=c(.data$SdG1^2,
                           .data$SdG1*.data$SdG2*.data$corG,
                           .data$SdG1*.data$SdG2*.data$corG,
                           .data$SdG2^2),
                    nrow= length(.data$coef_mean), ncol= length(.data$coef_mean),
                    dimnames=list(.data$label, .data$label))))%>%
    coef_sample(epsilon, type= type)
}

#' @rdname pull_lme
#' @export
pull_nlme= function(data, epsilon= 10^-6, type= c('sample','mean')){
  data%>%
    mutate(sigma= summary(.data$mod)$sigma,
           SdG1= as.numeric(nlme::VarCorr(.data$mod)[1,'StdDev']),
           SdG2= as.numeric(nlme::VarCorr(.data$mod)[2,'StdDev']),
           SdG3= as.numeric(nlme::VarCorr(.data$mod)[3,'StdDev']),
           corG12= as.numeric(nlme::VarCorr(.data$mod)[2,'Corr']),
           corG13= as.numeric(nlme::VarCorr(.data$mod)[3,'Corr']),
           corG23= as.numeric(utils::tail(nlme::VarCorr(.data$mod)[3,],1)),
           cov= list(stats::vcov(.data$mod)),
           coef_mean= list(nlme::fixed.effects(.data$mod)),
           label= list(colnames(nlme::random.effects(.data$mod))),
           MatG= list(
             matrix(
               data= c(.data$SdG1^2, .data$SdG1*.data$SdG2*.data$corG12,
                       .data$SdG1*.data$SdG3*.data$corG13,
                       .data$SdG1*.data$SdG2*.data$corG12,
                       .data$SdG2^2, .data$SdG2*.data$SdG3*.data$corG23,
                       .data$SdG1*.data$SdG3*.data$corG13,
                       .data$SdG2*.data$SdG3*.data$corG23, .data$SdG3^2),
               nrow=length(.data$coef_mean), ncol=length(.data$coef_mean),
               dimnames=list(.data$label, .data$label))))%>%
    coef_sample(epsilon, type= type)
}

#' @rdname pull_lme
#' @export
pull_others= function(data, type= c('sample','mean'),sd){
  type= match.arg(type,c('sample','mean'))
  if(nrow(data)>0){
    if(type=="sample"){
      data$coef_simu=
        stats::rnorm(n=1,mean= mean(data$coef_mean),
                     sd= sd)
    }else{
      data$coef_simu= mean(data$coef_mean)
    }
  }else{
    warning("No data to pull")
  }
  data
}
