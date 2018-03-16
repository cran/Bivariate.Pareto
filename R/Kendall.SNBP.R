#' Kendall's tau under the SNBP distribution
#'
#' @param Alpha0 Copula parameter \eqn{\alpha_{0}} with restricted range.
#' @param Alpha1 Positive scale parameter \eqn{\alpha_{1}} for the Pareto margin.
#' @param Alpha2 Positive scale parameter \eqn{\alpha_{2}} for the Pareto margin.
#' @param Gamma Common positive shape parameter \eqn{\gamma} for the Pareto margins.
#' @description Compute Kendall's tau under the Sankaran and Nair bivairate Pareto (SNBP) distribution (Sankaran and Nair, 1993) by numerical integration.
#' @details The admissible range of \code{Alpha0} (\eqn{\alpha_{0}}) is \eqn{0 \leq \alpha_{0} \leq (\gamma+1) \alpha_{1} \alpha_{2}.}
#' @return \item{tau}{Kendall's tau.}
#'
#' @references Sankaran and Nair (1993), A bivariate Pareto model and its applications to reliability, Naval Research Logistics, 40(7): 1013-1020.
#' @references Shih et al. (2018), Fitting competing risks data to bivariate Pareto models, Communications in Statistics - Theory and Methods, doi: 10.1080/03610926.2018.1425450.
#' @importFrom stats integrate
#' @export
#'
#' @examples
#' library(Bivariate.Pareto)
#' Kendall.SNBP(7e-5,0.0036,0.0075,1.8277)

Kendall.SNBP = function(Alpha0,Alpha1,Alpha2,Gamma) {

  ### checking inputs ###
  if (Alpha1 <= 0) {stop("Alpha1 must be positive")}
  if (Alpha2 <= 0) {stop("Alpha2 must be positive")}
  if (Gamma <= 0) {stop("Gamma must be positive")}
  if (Alpha0 > Alpha1*Alpha2*(Gamma+1) | Alpha0 < 0) {stop("Alpha0 is invalid")} else {Delta = Alpha0/(Alpha1*Alpha2)}

  douint = function(f,u.lower,u.upper,v.lower,v.upper) {

    f1 = function(u) integrate(f,lower = v.lower,upper = v.upper,u = u,rel.tol = 1e-9)$value
    f2 = function(u) sapply(u,f1)
    return(integrate(f2,u.lower,u.upper,rel.tol = 1e-9)$value)

  }

  delta = Alpha0/(Alpha1*Alpha2)

  Kendall.int = function(u,v) {

    A = (u^(-1/Gamma)+v^(-1/Gamma)-1+delta*(u^(-1/Gamma)-1)*(v^(-1/Gamma)-1))^(-Gamma-2)
    B = (1+delta*(u^(-1/Gamma)-1))*(1+delta*(v^(-1/Gamma)-1))
    C = (u^(-1/Gamma)+v^(-1/Gamma)-1+delta*(u^(-1/Gamma)-1)*(v^(-1/Gamma)-1))^(-Gamma-1)
    D = u^(-1/Gamma-1)*v^(-1/Gamma-1)

    den = (Gamma+1)/Gamma*D*B*A-delta/Gamma*D*C
    dis = (u^(-1/Gamma)+v^(-1/Gamma)-1+delta*(u^(-1/Gamma)-1)*(v^(-1/Gamma)-1))^(-Gamma)

    return(den*dis)

  }

  if (delta == 1) {

    return(0)

  } else if (Alpha0 == 0) {

    return((1/Gamma)/((1/Gamma)+2))

  } else {

    return(4*douint(Kendall.int,0,1,0,1)-1)

  }

}
