#' Generate samples from the SNBP distribution
#'
#' @param n Sample size.
#' @param Alpha0 Copula parameter \eqn{\alpha_{0}} with restricted range.
#' @param Alpha1 Positive scale parameter \eqn{\alpha_{1}} for the Pareto margin.
#' @param Alpha2 Positive scale parameter \eqn{\alpha_{2}} for the Pareto margin.
#' @param Gamma Common positive shape parameter \eqn{\gamma} for the Pareto margins.
#' @description Generate samples from the Sankaran and Nair bivairate Pareto (SNBP) distribution (Sankaran and Nair, 1993).
#' @details The admissible range of \code{Alpha0} (\eqn{\alpha_{0}}) is \eqn{0 \leq \alpha_{0} \leq (\gamma+1) \alpha_{1} \alpha_{2}.}
#'
#' @return \item{X}{\code{X} is asscoiated with the parameters \code{Alpha1} and \code{Gamma}.}
#' \item{Y}{\code{Y} is asscoiated with the parameters \code{Alpha2} and \code{Gamma}.}
#'
#' @references Sankaran and Nair (1993), A bivariate Pareto model and its applications to reliability, Naval Research Logistics, 40(7): 1013-1020.
#' @references Shih et al. (2018), Fitting competing risks data to bivariate Pareto models, Communications in Statistics - Theory and Methods, to appear.
#' @importFrom stats runif uniroot
#' @export
#'
#' @examples
#' library(Bivariate.Pareto)
#' SN.Pareto(5,2,1,1,1)

SN.Pareto = function(n,Alpha0,Alpha1,Alpha2,Gamma) {

  ### checking inputs ###
  if (n < 1 | n != round(n)) {stop("sample size n must be greater than or equal to 1 (integer)")}
  if (Alpha1 <= 0) {stop("Alpha1 must be positive")}
  if (Alpha2 <= 0) {stop("Alpha2 must be positive")}
  if (Gamma <= 0) {stop("Gamma must be positive")}
  if (Alpha0 > Alpha1*Alpha2*(Gamma+1) | Alpha0 < 0) {stop("Alpha0 is invalid")} else {Delta = Alpha0/(Alpha1*Alpha2)}

  U = rep(0,n)
  V = rep(0,n)
  for (i in 1:n) {

    u = runif(1)
    a = runif(1)
    f = function(v) {

      f1 = (u^(-1/Gamma)+v^(-1/Gamma)-1+Delta*(u^(-1/Gamma)-1)*(v^(-1/Gamma)-1))^(-Gamma-1)
      u^(-1/Gamma-1)*(1+Delta*(v^(-1/Gamma)-1))*f1-a

    }
    v = uniroot(f,c(1e-10,1),tol = 1e-9)

    U[i] = u
    V[i] = v$root

  }
  X = (U^(-1/Gamma)-1)/Alpha1
  Y = (V^(-1/Gamma)-1)/Alpha2

  return(cbind(X,Y))

}
