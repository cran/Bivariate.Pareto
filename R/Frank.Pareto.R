#' Generate samples from the Frank copula with the Pareto margins
#'
#' @param n Sample size.
#' @param Theta Copula parameter \eqn{\theta}.
#' @param Alpha1 Positive scale parameter \eqn{\alpha_{1}} for the Pareto margin.
#' @param Alpha2 Positive scale parameter \eqn{\alpha_{2}} for the Pareto margin.
#' @param Gamma1 Positive shape parameter \eqn{\gamma_{1}} for the Pareto margin.
#' @param Gamma2 Positive shape parameter \eqn{\gamma_{2}} for the Pareto margin.
#' @description Generate samples from the Frank copula with the Pareto margins.
#' @return \item{X}{\code{X} is asscoiated with the parameters \code{Alpha1} and \code{Gamma1}.}
#' \item{Y}{\code{Y} is asscoiated with the parameters \code{Alpha2} and \code{Gamma2}.}
#'
#' @references Shih et al. (2018), Fitting competing risks data to bivariate Pareto models, Communications in Statistics - Theory and Methods, to appear.
#' @importFrom stats runif
#' @export
#'
#' @examples
#' library(Bivariate.Pareto)
#' Frank.Pareto(5,5,1,1,1,1)

Frank.Pareto = function(n,Theta,Alpha1,Alpha2,Gamma1,Gamma2) {

  ### checking inputs ###
  if (n < 1 | n != round(n)) {stop("sample size n must be greater than or equal to 1 (integer)")}
  if (Theta == 0) {stop("Theta cannot be zero")}
  if (Alpha1 <= 0) {stop("Alpha1 must be positive")}
  if (Alpha2 <= 0) {stop("Alpha2 must be positive")}
  if (Gamma1 <= 0) {stop("Gamma1 must be positive")}
  if (Gamma2 <= 0) {stop("Gamma2 must be positive")}

  U = runif(n)
  a = runif(n)
  V = (-1/Theta)*log(1+a*(exp(-Theta)-1)/(exp(-Theta*U)-a*(exp(-Theta*U)-1)))
  X = (1/Alpha1)*(U^(-1/Gamma1)-1)
  Y = (1/Alpha2)*(V^(-1/Gamma2)-1)

  return(cbind(X,Y))

}
