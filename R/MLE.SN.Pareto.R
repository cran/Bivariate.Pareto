#' Maximum likelihood estimation for bivariate dependent competing risks data under the SNBP distribution
#'
#' @param t.event Vector of the observed failure times.
#' @param event1 Vector of the indicators for the failure cause 1.
#' @param event2 Vector of the indicators for the failure cause 2.
#' @param Alpha0 Copula parameter \eqn{\alpha_{0}} with restricted range.
#' @param Alpha1.0 Initial guess for the scale parameter \eqn{\alpha_{1}} with default value 1.
#' @param Alpha2.0 Initial guess for the scale parameter \eqn{\alpha_{2}} with default value 1.
#' @param Gamma.0 Initial guess for the common shape parameter \eqn{\gamma} with default value 1.
#' @param epsilon Positive tunning parameter in the NR algorithm with default value \eqn{10^{-5}}.
#' @param d Positive tunning parameter in the NR algorithm with default value \eqn{e^{10}}.
#' @param r.1 Positive tunning parameter in the NR algorithm with default value 1.
#' @param r.2 Positive tunning parameter in the NR algorithm with default value 1.
#' @param r.3 Positive tunning parameter in the NR algorithm with default value 1.
#' @description Maximum likelihood estimation for bivariate dependent competing risks data under the SNBP distribution (Sankaran and Nair, 1993).
#' @details The admissible range of \code{Alpha0} (\eqn{\alpha_{0}}) is \eqn{0 \leq \alpha_{0} \leq (\gamma+1) \alpha_{1} \alpha_{2}.}
#'
#' To adapt our functions to dependent censoring models in Emura and Chen (2018), one can simply set \code{event2} = \code{1-event1} (See examples).
#'
#' @return \item{n}{Sample size.}
#' \item{count}{Iteration number.}
#' \item{random}{Randomization number.}
#' \item{Alpha1}{Positive scale parameter for the Pareto margin (failure cause 1).}
#' \item{Alpha2}{Positive scale parameter for the Pareto margin (failure cause 2).}
#' \item{Gamma}{Common positive shape parameter for the Pareto margins.}
#' \item{MedX}{Median lifetime due to failure cause 1.}
#' \item{MedY}{Median lifetime due to failure cause 2.}
#' \item{MeanX}{Mean lifetime due to failure cause 1.}
#' \item{MeanY}{Mean lifetime due to failure cause 2.}
#' \item{logL}{Log-likelihood value under the fitted model.}
#' \item{AIC}{AIC value under the fitted model.}
#' \item{BIC}{BIC value under the fitted model.}
#'
#' @references Sankaran PG, Nair NU (1993), A bivariate Pareto model and its applications to reliability, Naval Research Logistics, 40(7): 1013-1020.
#' @references Shih J-H, Lee W, Sun L-H, Emura T (2018), Fitting competing risks data to bivariate Pareto models, Communications in Statistics - Theory and Methods, doi: 10.1080/03610926.2018.1425450.
#' @references Emura T, Chen Y-H (2018) Analysis of Survival Data with Dependent Censoring, Copula-Based Approaches, JSS Research Series in Statistics, Springer, in press.
#' @importFrom stats qnorm runif
#' @importFrom utils globalVariables
#' @import compound.Cox
#' @export
#'
#' @examples
#' t.event = c(72,40,20,65,24,46,62,61,60,60,59,59,49,20, 3,58,29,26,52,20,
#'             51,51,31,42,38,69,39,33, 8,13,33, 9,21,66, 5,27, 2,20,19,60,
#'             32,53,53,43,21,74,72,14,33, 8,10,51, 7,33, 3,43,37, 5, 6, 2,
#'             5,64, 1,21,16,21,12,75,74,54,73,36,59, 6,58,16,19,39,26,60,
#'             43, 7, 9,67,62,17,25, 0, 5,34,59,31,58,30,57, 5,55,55,52, 0,
#'             51,17,70,74,74,20, 2, 8,27,23, 1,52,51, 6, 0,26,65,26, 6, 6,
#'             68,33,67,23, 6,11, 6,57,57,29, 9,53,51, 8, 0,21,27,22,12,68,
#'             21,68, 0, 2,14,18, 5,60,40,51,50,46,65, 9,21,27,54,52,75,30,
#'             70,14, 0,42,12,40, 2,12,53,11,18,13,45, 8,28,67,67,24,64,26,
#'             57,32,42,20,71,54,64,51, 1, 2, 0,54,69,68,67,66,64,63,35,62,
#'             7,35,24,57, 1, 4,74, 0,51,36,16,32,68,17,66,65,19,41,28, 0,
#'             46,63,60,59,46,63, 8,74,18,33,12, 1,66,28,30,57,50,39,40,24,
#'             6,30,58,68,24,33,65, 2,64,19,15,10,12,53,51, 1,40,40,66, 2,
#'             21,35,29,54,37,10,29,71,12,13,27,66,28,31,12, 9,21,19,51,71,
#'             76,46,47,75,75,49,75,75,31,69,74,25,72,28,36, 8,71,60,14,22,
#'             67,62,68,68,27,68,68,67,67, 3,49,12,30,67, 5,65,24,66,36,66,
#'             40,13,40, 0,14,45,64,13,24,15,26, 5,63,35,61,61,50,57,21,26,
#'             11,59,42,27,50,57,57, 0, 1,54,53,23, 8,51,27,52,52,52,45,48,
#'             18, 2, 2,35,75,75, 9,39, 0,26,17,43,53,47,11,65,16,21,64, 7,
#'             38,55, 5,28,38,20,24,27,31, 9, 9,11,56,36,56,15,51,33,70,32,
#'             5,23,63,30,53,12,58,54,36,20,74,34,70,25,65, 4,10,58,37,56,
#'             6, 0,70,70,28,40,67,36,23,23,62,62,62, 2,34, 4,12,56, 1, 7,
#'             4,70,65, 7,30,40,13,22, 0,18,64,13,26, 1,16,33,22,30,53,53,
#'             7,61,40, 9,59, 7,12,46,50, 0,52,19,52,51,51,14,27,51, 5, 0,
#'             41,53,19)
#'
#' event1 = c(0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
#'            0,0,1,0,0,0,1,0,1,1,0,1,1,1,1,0,0,1,1,0,
#'            1,0,0,1,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,1,
#'            1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
#'            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'            0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,1,0,0,
#'            0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,
#'            0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'            0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,
#'            1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
#'            0,0,0,0,0,0,0,1,0,0,1,1,0,1,0,0,1,1,0,0,
#'            1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
#'            0,0,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,0,0,
#'            1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,0,1,1,0,1,
#'            0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,
#'            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
#'            0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,
#'            0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
#'            1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,0,1,
#'            1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,
#'            0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,1,1,0,1,0,
#'            1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,
#'            1,0,0,1,0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,1,
#'            1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,
#'            0,0,1)
#'
#' event2 = c(0,1,1,0,0,1,0,0,0,0,0,0,0,1,1,0,1,1,0,1,
#'            0,0,0,1,1,0,0,1,0,0,1,0,0,0,0,1,1,0,0,0,
#'            0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,
#'            0,0,1,0,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,
#'            1,1,1,0,1,1,1,1,1,1,0,1,0,1,0,1,0,0,0,1,
#'            0,1,1,0,0,1,0,0,1,1,1,0,0,0,0,1,1,0,1,1,
#'            0,1,0,0,1,1,0,0,0,1,1,0,0,1,1,1,0,1,0,0,
#'            1,0,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,0,0,1,
#'            0,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,0,1,
#'            0,0,1,1,0,1,0,1,1,1,0,1,0,0,0,0,0,0,1,0,
#'            1,1,1,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,
#'            0,0,0,0,1,0,1,0,1,1,1,1,0,1,1,1,0,1,1,1,
#'            1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,
#'            0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
#'            0,0,1,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,1,1,
#'            1,1,0,0,1,0,0,0,0,1,1,1,1,0,1,1,1,0,1,0,
#'            1,1,1,1,1,1,0,1,1,1,1,0,0,1,0,0,1,1,1,0,
#'            1,0,0,1,1,0,0,1,1,0,0,1,1,1,1,0,0,0,1,1,
#'            0,1,1,1,0,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,
#'            0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,
#'            1,1,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
#'            0,1,0,0,1,1,0,1,1,1,0,0,0,1,0,1,0,0,1,1,
#'            0,0,0,0,1,1,1,0,1,0,1,1,0,1,1,1,0,0,1,0,
#'            0,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,1,1,
#'            1,0,0)
#'
#' library(Bivariate.Pareto)
#' set.seed(10)
#' MLE.SN.Pareto(t.event,event1,event2,Alpha0 = 7e-5)
#'
#' ### dependent censoring ###
#'
#' library(compound.Cox)
#' library(Bivariate.Pareto)
#' data(Lung)
#' t.vec = Lung$t.vec
#' d.vec = Lung$d.vec
#'
#' MLE.SN.Pareto(t.vec,d.vec,1-d.vec,Alpha0 = 1e-5,d = exp(15))

MLE.SN.Pareto = function(t.event,event1,event2,Alpha0,Alpha1.0 = 1,Alpha2.0 = 1,Gamma.0 = 1,
                         epsilon = 1e-5,d = exp(10),r.1 = 6,r.2 = 6,r.3 = 6) {

  ### checking inputs ###
  n = length(t.event)
  if (length(t.event[t.event < 0]) != 0) {stop("t.event must be non-negative")}
  if (length(event1) != n) {stop("the length of event1 is different from t.event")}
  if (length(event2) != n) {stop("the length of event2 is different from t.event")}
  if (length(event1[event1 == 0 | event1 == 1]) != n) {stop("elements in event1 must be either 0 or 1")}
  if (length(event2[event2 == 0 | event2 == 1]) != n) {stop("elements in event2 must be either 0 or 1")}
  temp.event = event1+event2
  if (length(temp.event[temp.event == 2]) != 0) {stop("event1 and event2 cannot be 1 simultaneously")}
  if (Alpha0 < 0) {stop("Alpha0 cannot be negative")}
  if (Alpha1.0 <= 0) {stop("Alpha1.0 must be positive")}
  if (Alpha2.0 <= 0) {stop("Alpha2.0 must be positive")}
  if (Gamma.0 <= 0) {stop("Gamma.0 must be positive")}
  if (epsilon <= 0) {stop("epsilon must be positive")}
  if (d <= 0) {stop("d must be positive")}
  if (r.1 <= 0) {stop("r.1 must be positive")}
  if (r.2 <= 0) {stop("r.2 must be positive")}
  if (r.3 <= 0) {stop("r.3 must be positive")}

  ### functions ###
  log_L = function(par) {

    Alpha1 = exp(par[1])
    Alpha2 = exp(par[2])
    Gamma  = exp(par[3])

    p0 = 1+Alpha1*t.event+Alpha2*t.event+Alpha0*t.event^2
    p1 = Alpha1+Alpha0*t.event
    p2 = Alpha2+Alpha0*t.event

    L1 = sum( event1+event2 )*log( Gamma )
    L2 = sum( event1*log(p1) )
    L3 = sum( event2*log(p2) )
    L4 = sum((Gamma+event1+event2)*log(p0))

    L1+L2+L3-L4

  }
  SL_function = function(par) {

    Alpha1 = exp(par[1])
    Alpha2 = exp(par[2])
    Gamma  = exp(par[3])

    p0 = 1+Alpha1*t.event+Alpha2*t.event+Alpha0*t.event^2
    p1 = Alpha1+Alpha0*t.event
    p2 = Alpha2+Alpha0*t.event

    d1 = sum(event1/p1)-sum((Gamma+event1+event2)*t.event/p0)
    d2 = sum(event2/p2)-sum((Gamma+event1+event2)*t.event/p0)
    d3 = sum(event1+event2)/Gamma-sum(log(p0))

    c(exp(par[1])*d1,exp(par[2])*d2,exp(par[3])*d3)

  }
  HL_function = function(par) {

    Alpha1 = exp(par[1])
    Alpha2 = exp(par[2])
    Gamma  = exp(par[3])

    p0 = 1+Alpha1*t.event+Alpha2*t.event+Alpha0*t.event^2
    p1 = Alpha1+Alpha0*t.event
    p2 = Alpha2+Alpha0*t.event

    d1 = sum(event1/p1)-sum((Gamma+event1+event2)*t.event/p0)
    d2 = sum(event2/p2)-sum((Gamma+event1+event2)*t.event/p0)
    d3 = sum(event1+event2)/Gamma-sum(log(p0))

    d11 = sum((event1+event2+Gamma)*t.event^2/p0^2-event1/p1^2)
    d22 = sum((event1+event2+Gamma)*t.event^2/p0^2-event2/p2^2)
    d33 = sum(-(event1+event2))/Gamma^2
    d23 = sum(-t.event/p0)
    d12  =sum((event1+event2+Gamma)*t.event^2/p0^2)
    d13 = sum(-t.event/p0)

    D11 = exp(2*par[1])*d11+exp(par[1])*d1
    D12 = exp(par[1])*exp(par[2])*d12
    D13 = exp(par[1])*exp(par[3])*d13
    D22 = exp(2*par[2])*d22+exp(par[2])*d2
    D23 = exp(par[2])*exp(par[3])*d23
    D33 = exp(2*par[3])*d33+exp(par[3])*d3

    matrix(c(D11,D12,D13,D12,D22,D23,D13,D23,D33),3,3)

  }
  H_function = function(par) {

    Alpha1 = par[1]
    Alpha2 = par[2]
    Gamma  = par[3]

    p0 = 1+Alpha1*t.event+Alpha2*t.event+Alpha0*t.event^2
    p1 = Alpha1+Alpha0*t.event
    p2 = Alpha2+Alpha0*t.event

    d1 = sum(event1/p1)-sum((Gamma+event1+event2)*t.event/p0)
    d2 = sum(event2/p2)-sum((Gamma+event1+event2)*t.event/p0)
    d3 = sum(event1+event2)/Gamma-sum(log(p0))

    d11 = sum((event1+event2+Gamma)*t.event^2/p0^2-event1/p1^2)
    d22 = sum((event1+event2+Gamma)*t.event^2/p0^2-event2/p2^2)
    d33 = sum(-(event1+event2))/Gamma^2
    d23 = sum(-t.event/p0)
    d12 = sum((event1+event2+Gamma)*t.event^2/p0^2)
    d13 = sum(-t.event/p0)

    matrix(c(d11,d12,d13,d12,d22,d23,d13,d23,d33),3,3)

  }

  par_old = c(log(Alpha1.0),log(Alpha2.0),log(Gamma.0))
  count  = 0
  random = 0

  repeat{

    temp = try(solve(HL_function(par_old),silent = TRUE))
    if (class(temp) == "try-error"){

      random = random+1
      count = 0
      par_old = c(log(Alpha1.0*exp(runif(1,-r.1,r.1))),
                  log(Alpha2.0*exp(runif(1,-r.2,r.2))),
                  log(Gamma.0*exp(runif(1,-r.3,r.3))))
      next

    }

    par_new = par_old-solve(HL_function(par_old))%*%SL_function(par_old)
    count = count+1

    if (is.na(sum(par_new)) |
        max(abs(par_new)) > log(d)) {

      random = random+1
      count  = 0
      par_old = c(log(Alpha1.0*exp(runif(1,-r.1,r.1))),
                  log(Alpha2.0*exp(runif(1,-r.2,r.2))),
                  log(Gamma.0*exp(runif(1,-r.3,r.3))))
      next

    }

    if (max(abs(exp(par_old)-exp(par_new))) < epsilon) {break}
    par_old = par_new

  }

  Alpha1_hat = exp(par_new[1])
  Alpha2_hat = exp(par_new[2])
  Gamma_hat  = exp(par_new[3])

  Info = solve(-H_function(exp(par_new)))
  Alpha1_se = sqrt(Info[1,1])
  Alpha2_se = sqrt(Info[2,2])
  Gamma_se  = sqrt(Info[3,3])

  InfoL = solve(-HL_function(par_new))
  CI_Alpha1 = c(Alpha1_hat*exp(-qnorm(0.975)*sqrt(InfoL[1,1])),
                Alpha1_hat*exp(+qnorm(0.975)*sqrt(InfoL[1,1])))
  CI_Alpha2 = c(Alpha2_hat*exp(-qnorm(0.975)*sqrt(InfoL[2,2])),
                Alpha2_hat*exp(+qnorm(0.975)*sqrt(InfoL[2,2])))
  CI_Gamma  = c(Gamma_hat*exp(-qnorm(0.975)*sqrt(InfoL[3,3])),
                Gamma_hat*exp(+qnorm(0.975)*sqrt(InfoL[3,3])))

  MedX_hat = (2^(1/Gamma_hat)-1)/Alpha1_hat
  MedY_hat = (2^(1/Gamma_hat)-1)/Alpha2_hat
  transX = c((1-2^(1/Gamma_hat))/Alpha1_hat^2,0,-2^(1/Gamma_hat)*log(2)/(Alpha1_hat*Gamma_hat^2))
  transY = c(0,(1-2^(1/Gamma_hat))/Alpha2_hat^2,-2^(1/Gamma_hat)*log(2)/(Alpha2_hat*Gamma_hat^2))
  MedX_se = sqrt(t(transX)%*%Info%*%transX)
  MedY_se = sqrt(t(transY)%*%Info%*%transY)

  temp_transX = c(-1,0,-2^(1/Gamma_hat)*log(2)/((2^(1/Gamma_hat)-1)*Gamma_hat))
  temp_transY = c(0,-1,-2^(1/Gamma_hat)*log(2)/((2^(1/Gamma_hat)-1)*Gamma_hat))
  temp_MedX_se = sqrt(t(temp_transX)%*%InfoL%*%temp_transX)
  temp_MedY_se = sqrt(t(temp_transY)%*%InfoL%*%temp_transY)

  CI_MedX = c(MedX_hat*exp(-qnorm(0.975)*temp_MedX_se),
              MedX_hat*exp(+qnorm(0.975)*temp_MedX_se))
  CI_MedY = c(MedY_hat*exp(-qnorm(0.975)*temp_MedY_se),
              MedY_hat*exp(+qnorm(0.975)*temp_MedY_se))

  Alpha1.res = c(Estimate = Alpha1_hat,SE = Alpha1_se,CI.lower = CI_Alpha1[1],CI.upper = CI_Alpha1[2])
  Alpha2.res = c(Estimate = Alpha2_hat,SE = Alpha2_se,CI.lower = CI_Alpha2[1],CI.upper = CI_Alpha2[2])
  Gamma.res  = c(Estimate = Gamma_hat, SE = Gamma_se, CI.lower = CI_Gamma[1], CI.upper = CI_Gamma[2] )

  MedX.res = c(Estimate = MedX_hat,SE = MedX_se,CI.lower = CI_MedX[1],CI.upper = CI_MedX[2])
  MedY.res = c(Estimate = MedY_hat,SE = MedY_se,CI.lower = CI_MedY[1],CI.upper = CI_MedY[2])

  if (Gamma_hat < 1) {

    return(list(n = n,Iteration = count,Randomization = random,
                Alpha1 = Alpha1.res,Alpha2 = Alpha2.res,Gamma = Gamma.res,
                MedX = MedX.res,MedY = MedY.res,MeanX = "Unavaliable",MeanY = "Unavaliable",
                logL = log_L(par_new),AIC = 2*length(par_new)-2*log_L(par_new),
                BIC = length(par_new)*log(length(t.event))-2*log_L(par_new)))

  } else {

    MeanX_hat = 1/(Alpha1_hat*(Gamma_hat-1))
    MeanY_hat = 1/(Alpha2_hat*(Gamma_hat-1))
    trans2X = c(-1/(Alpha1_hat^2*(Gamma_hat-1)),0,-1/(Alpha1_hat*(Gamma_hat-1)^2))
    trans2Y = c(0,-1/(Alpha2_hat^2*(Gamma_hat-1)),-1/(Alpha2_hat*(Gamma_hat-1)^2))
    MeanX_se = sqrt(t(trans2X)%*%Info%*%trans2X)
    MeanY_se = sqrt(t(trans2Y)%*%Info%*%trans2Y)

    temp_trans2X = c(-1,0,-Gamma_hat/(Gamma_hat-1))
    temp_trans2Y = c(0,-1,-Gamma_hat/(Gamma_hat-1))
    temp_MeanX_se = sqrt(t(temp_trans2X)%*%InfoL%*%temp_trans2X)
    temp_MeanY_se = sqrt(t(temp_trans2Y)%*%InfoL%*%temp_trans2Y)

    CI_MeanX = c(MeanX_hat*exp(-qnorm(0.975)*temp_MeanX_se),
                 MeanX_hat*exp(+qnorm(0.975)*temp_MeanX_se))
    CI_MeanY = c(MeanY_hat*exp(-qnorm(0.975)*temp_MeanY_se),
                 MeanY_hat*exp(+qnorm(0.975)*temp_MeanY_se))

    MeanX.res = c(Estimate = MeanX_hat,SE = MeanX_se,CI.lower = CI_MeanX[1],CI.upper = CI_MeanX[2])
    MeanY.res = c(Estimate = MeanY_hat,SE = MeanY_se,CI.lower = CI_MeanY[1],CI.upper = CI_MeanY[2])

    return(list(n = n,Iteration = count,Randomization = random,
                Alpha1 = Alpha1.res,Alpha2 = Alpha2.res,Gamma = Gamma.res,
                MedX = MedX.res,MedY = MedY.res,
                MeanX = MeanX.res,MeanY = MeanY.res,logL = log_L(par_new),
                AIC = 2*length(par_new)-2*log_L(par_new),
                BIC = length(par_new)*log(length(t.event))-2*log_L(par_new)))

  }

}

