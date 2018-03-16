#' Maximum likelihood estimation for bivariate dependent competing risks data under the Frank copula with the common Pareto margins
#'
#' @param t.event Vector of the observed failure times.
#' @param event1 Vector of the indicators for the failure cause 1.
#' @param event2 Vector of the indicators for the failure cause 2.
#' @param Theta.0 Initial guess for the copula parameter \eqn{\theta}.
#' @param Alpha.0 Initial guess for the common scale parameter \eqn{\alpha} with default value 1.
#' @param Gamma.0 Initial guess for the common shape parameter \eqn{\gamma} with default value 1.
#' @param epsilon Positive tunning parameter in the NR algorithm with default value \eqn{10^{-5}}.
#' @param r.1 Positive tunning parameter in the NR algorithm with default value 1.
#' @param r.2 Positive tunning parameter in the NR algorithm with default value 1.
#' @param r.3 Positive tunning parameter in the NR algorithm with default value 1.
#' @param bootstrap Perform parametric bootstrap if \code{TRUE}.
#' @param B Number of bootstrap replications.
#' @description Maximum likelihood estimation for bivariate dependent competing risks data under the Frank copula with the common Pareto margins.
#' @details The parametric bootstrap method requires the assumption of the uniform censoring distribution. One must notice that such assumption is not always true in real data analysis.
#'
#' @return \item{n}{Sample size.}
#' \item{count}{Iteration number.}
#' \item{random}{Randomization number.}
#' \item{Theta}{Copula parameter.}
#' \item{Theta.B}{Copula parameter (SE and CI are calculated by parametric bootstrap method).}
#' \item{Alpha}{Common positive scale parameter for the Pareto margin.}
#' \item{Alpha.B}{Common positive scale parameter for the Pareto margin (SE and CI are calculated by parametric bootstrap method).}
#' \item{Gamma}{Common positive shape parameter for the Pareto margin.}
#' \item{Gamma.B}{Common positive shape parameter for the Pareto margin (SE and CI are calculated by parametric bootstrap method).}
#' \item{logL}{Log-likelihood value under the fitted model.}
#' \item{AIC}{AIC value under the fitted model.}
#' \item{BIC}{BIC value under the fitted model.}
#'
#' @references Shih et al. (2018), Fitting competing risks data to bivariate Pareto models, Communications in Statistics - Theory and Methods, doi: 10.1080/03610926.2018.1425450.
#' @importFrom stats qnorm runif sd
#' @importFrom utils globalVariables
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
#' MLE.Frank.Pareto.com(t.event,event1,event2,bootstrap = FALSE)

MLE.Frank.Pareto.com = function(t.event,event1,event2,Theta.0 = 1,Alpha.0 = 1,Gamma.0 = 1,
                                epsilon = 1e-5,r.1 = 13,r.2 = 3,r.3 = 3,
                                bootstrap = FALSE,B = 200) {

  ### checking inputs ###
  n = length(t.event)
  if (length(t.event[t.event < 0]) != 0) {stop("t.event must be non-negative")}
  if (length(event1) != n) {stop("the length of event1 is different from t.event")}
  if (length(event2) != n) {stop("the length of event2 is different from t.event")}
  if (length(event1[event1 == 0 | event1 == 1]) != n) {stop("elements in event1 must be either 0 or 1")}
  if (length(event2[event2 == 0 | event2 == 1]) != n) {stop("elements in event2 must be either 0 or 1")}
  temp.event = event1+event2
  if (length(temp.event[temp.event == 2]) != 0) {stop("event1 and event2 cannot be 1 simultaneously")}
  if (Theta.0 == 0) {stop("Theta.0 cannot be zero")}
  if (Alpha.0 <= 0) {stop("Alpha.0 must be positive")}
  if (Gamma.0 <= 0) {stop("Gamma.0 must be positive")}
  if (epsilon <= 0) {stop("epsilon must be positive")}
  if (r.1 <= 0) {stop("r.1 must be positive")}
  if (r.2 <= 0) {stop("r.2 must be positive")}
  if (r.3 <= 0) {stop("r.3 must be positive")}
  if (B != round(B)) {stop("B must be integer")}

  ### functions ###
  log_L = function(par) {

    Theta = par[1]
    Alpha = exp(par[2])
    Gamma = exp(par[3])

    h  = Alpha*Gamma/(1+Alpha*t.event)
    S  = (1+Alpha*t.event)^(-Gamma)
    p0 = exp(-Theta)
    p1 = exp(-Theta*S)
    p2 = 1+(p1-1)^2/(p0-1)
    ST = -(1/Theta)*log(p2)

    M1 = (event1+event2)*(log(h)+log(S)-Theta*S+log((p1-1)/(p0-1))+Theta*ST)
    M2 = (1-event1-event2)*log(-log(p2)/Theta)

    sum(M1)+sum(M2)

  }
  SL_function = function(par) {

    Theta = par[1]
    Alpha = exp(par[2])
    Gamma = exp(par[3])

    h  = Alpha*Gamma/(1+Alpha*t.event)
    S  = (1+Alpha*t.event)^(-Gamma)
    p0 = exp(-Theta)
    p1 = exp(-Theta*S)
    p2 = 1+(p1-1)^2/(p0-1)
    ST = -(1/Theta)*log(p2)

    der_h_Alpha = Gamma/(1+Alpha*t.event)^2
    der_S_Alpha = -Gamma*t.event*(1+Alpha*t.event)^(-Gamma-1)
    der_h_Gamma = Alpha/(1+Alpha*t.event)
    der_S_Gamma = -(1+Alpha*t.event)^(-Gamma)*log(1+Alpha*t.event)

    der_ST_Theta = Theta^(-2)*log(p2)-(p1-1)*(-2*S*p1*(p0-1)+p0*(p1-1))/(Theta*(p0-1+(p1-1)^2)*(p0-1))
    der_ST_Alpha = (p1-1)*2*der_S_Alpha*p1/(p0-1+(p1-1)^2)
    der_ST_Gamma = (p1-1)*2*der_S_Gamma*p1/(p0-1+(p1-1)^2)

    d0_1 = event1*(-S-S*p1/(p1-1)+p0/(p0-1)+ST+Theta*der_ST_Theta)
    d0_2 = event2*(-S-S*p1/(p1-1)+p0/(p0-1)+ST+Theta*der_ST_Theta)
    d0_3 = (1-event1-event2)*((-1/Theta)+((p0-1)*(-2)*S*p1*(p1-1)+p0*(p1-1)^2)/((p0-1)*(p0-1+(p1-1)^2)*log(p2)))
    d0   = sum(d0_1)+sum(d0_2)+sum(d0_3)

    d1_1 = event1*(der_h_Alpha/h+der_S_Alpha/S-Theta*der_S_Alpha-(der_S_Alpha*Theta*p1)/(p1-1)+Theta*der_ST_Alpha)
    d1_2 = event2*(der_h_Alpha/h+der_S_Alpha/S-Theta*der_S_Alpha-(der_S_Alpha*Theta*p1)/(p1-1)+Theta*der_ST_Alpha)
    d1_3 = (1-event1-event2)*(-2*Theta*der_S_Alpha*p1*(p1-1)/(log(p2)*(p0-1+(p1-1)^2)))
    d1  = sum(d1_1)+sum(d1_2)+sum(d1_3)

    d2_1 = event1*(der_h_Gamma/h+der_S_Gamma/S-Theta*der_S_Gamma-(der_S_Gamma*Theta*p1)/(p1-1)+Theta*der_ST_Gamma)
    d2_2 = event2*(der_h_Gamma/h+der_S_Gamma/S-Theta*der_S_Gamma-(der_S_Gamma*Theta*p1)/(p1-1)+Theta*der_ST_Gamma)
    d2_3 = (1-event1-event2)*(-2*Theta*der_S_Gamma*p1*(p1-1)/(log(p2)*(p0-1+(p1-1)^2)))
    d2   = sum(d2_1)+sum(d2_2)+sum(d2_3)

    s0 = d0
    s1 = d1*exp(par[2])
    s2 = d2*exp(par[3])

    c(s0,s1,s2)

  }
  HL_function = function (par) {

    Theta = par[1]
    Alpha = exp(par[2])
    Gamma = exp(par[3])

    h = Alpha*Gamma/(1+Alpha*t.event)
    S = (1+Alpha*t.event)^(-Gamma)
    p0 = exp(-Theta)
    p1 = exp(-Theta*S)
    p2 = (1 + (p1 - 1)^2/(p0 - 1))
    ST = -(1/Theta)*log(1+(exp(-Theta*S)-1)*(exp(-Theta*S)-1)/(exp(-Theta)-1))

    der_h_Alpha = Gamma/(1+Alpha*t.event)^2
    der_S_Alpha = -Gamma*t.event*(1+Alpha*t.event)^(-Gamma-1)
    der_h_Gamma = Alpha/(1+Alpha*t.event)
    der_S_Gamma = -(1+Alpha*t.event)^(-Gamma)*log(1+Alpha*t.event)

    der_ST_Theta = (Theta^(-2)*log(1+(p1-1)^2/(p0-1))-(p1-1)*(-2*S*p1*(p0-1)+p0*(p1-1))/(Theta*(p0-1)*((p1-1)^2+p0-1)))
    der_ST_Alpha = (p1-1)*2*der_S_Alpha*p1/(p0-1+(p1-1)^2)
    der_ST_Gamma = (p1-1)*2*der_S_Gamma*p1/(p0-1+(p1-1)^2)

    der_h_Alpha_Alpha = -2*Gamma*t.event/(1+Alpha*t.event)^3
    der_S_Alpha_Alpha = Gamma*(Gamma+1)*t.event^2*(1+Alpha*t.event)^(-Gamma-2)
    der_h_Gamma_Gamma = 0
    der_S_Gamma_Gamma = (1+Alpha*t.event)^(-Gamma)*(log(1+Alpha*t.event))^2
    der_h_Alpha_Gamma = (1+Alpha*t.event)^(-2)
    der_S_Alpha_Gamma = t.event*(Gamma*log(1+Alpha*t.event)-1)/(1+Alpha*t.event)^(Gamma+1)

    der_ST_Theta_Theta = -2*Theta^(-3)*log(p2)+(p0*(p1-1)^2-2*S*p1*(p1-1)*(p0-1))/(Theta^2*(p0-1)*(p0-1+(p1-1)^2))-((Theta*(p0-1)*(p0-1+(p1-1)^2))*(2*S^2*p1*(p0-1)*(2*p1-1)-p0*(p1-1)^2)-(-2*S*p1*(p1-1)*(p0-1)+p0*(p1-1)^2)*((p0-1)*(p0-1+(p1-1)^2)+Theta*(-p0*(p0-1+(p1-1)^2)+(p0-1)*(-p0-2*S*p1*(p1-1)))))/(Theta*(p0-1)*(p0-1+(p1-1)^2))^2
    der_ST_Alpha_Alpha = ((p0-1+(p1-1)^2)*(2*der_S_Alpha_Alpha*p1*(p1-1)+2*der_S_Alpha*(-2*Theta*der_S_Alpha*p1^2+Theta*der_S_Alpha*p1))-2*der_S_Alpha*(p1^2-p1)*(-2*Theta*der_S_Alpha*p1*(p1-1)))/(p0-1+(p1-1)^2)^2
    der_ST_Gamma_Gamma = ((p0-1+(p1-1)^2)*(2*der_S_Gamma_Gamma*p1*(p1-1)+2*der_S_Gamma*(-2*Theta*der_S_Gamma*p1^2+Theta*der_S_Gamma*p1))-2*der_S_Gamma*(p1^2-p1)*(-2*Theta*der_S_Gamma*p1*(p1-1)))/(p0-1+(p1-1)^2)^2
    der_ST_Theta_Alpha = ((p0-1+(p1-1)^2)*(2*der_S_Alpha*S*p1*(1-2*p1))-(2*der_S_Alpha*p1*(p1-1))*(-p0-2*S*p1*(p1-1)))/(p0-1+(p1-1)^2)^2
    der_ST_Theta_Gamma = ((p0-1+(p1-1)^2)*(2*der_S_Gamma*S*p1*(1-2*p1))-(2*der_S_Gamma*p1*(p1-1))*(-p0-2*S*p1*(p1-1)))/(p0-1+(p1-1)^2)^2
    der_ST_Alpha_Gamma = ((p0-1+(p1-1)^2)*(2*der_S_Alpha_Gamma*p1*(p1-1)+2*der_S_Alpha*(-2*Theta*der_S_Gamma*p1^2+Theta*der_S_Gamma*p1))-2*der_S_Alpha*(p1^2-p1)*(-2*Theta*der_S_Gamma*p1*(p1-1)))/(p0-1+(p1-1)^2)^2

    D0_1 = event1*(-S-S*p1/(p1-1)+p0/(p0-1)+ST+Theta*der_ST_Theta)
    D0_2 = event2*(-S-S*p1/(p1-1)+p0/(p0-1)+ST+Theta*der_ST_Theta)
    D0_3 = (1-event1-event2)*((-1/Theta)+((p0-1)*(-2)*S*p1*(p1-1)+p0*(p1-1)^2)/((p0-1)*(p0-1+(p1-1)^2)*log(p2)))
    D0   = sum(D0_1)+sum(D0_2)+sum(D0_3)

    D1_1 = event1*(der_h_Alpha/h+der_S_Alpha/S-Theta*der_S_Alpha-(der_S_Alpha*Theta*p1)/(p1-1)+Theta*der_ST_Alpha)
    D1_2 = event2*(der_h_Alpha/h+der_S_Alpha/S-Theta*der_S_Alpha-(der_S_Alpha*Theta*p1)/(p1-1)+Theta*der_ST_Alpha)
    D1_3 = (1-event1-event2)*(-2*Theta*der_S_Alpha*p1*(p1-1)/(log(p2)*(p0-1+(p1-1)^2)))
    D1  = sum(D1_1)+sum(D1_2)+sum(D1_3)

    D2_1 = event1*(der_h_Gamma/h+der_S_Gamma/S-Theta*der_S_Gamma-(der_S_Gamma*Theta*p1)/(p1-1)+Theta*der_ST_Gamma)
    D2_2 = event2*(der_h_Gamma/h+der_S_Gamma/S-Theta*der_S_Gamma-(der_S_Gamma*Theta*p1)/(p1-1)+Theta*der_ST_Gamma)
    D2_3 = (1-event1-event2)*(-2*Theta*der_S_Gamma*p1*(p1-1)/(log(p2)*(p0-1+(p1-1)^2)))
    D2   = sum(D2_1)+sum(D2_2)+sum(D2_3)

    d0_1 = event1*(-S-S*p1/(p1-1)+exp(-Theta)/(exp(-Theta)-1)+ST+Theta*der_ST_Theta)
    d0_2 = event2*(-S-S*p1/(p1-1)+exp(-Theta)/(exp(-Theta)-1)+ST+Theta*der_ST_Theta)
    d0_3 = (1-event1-event2)*((-1/Theta)+((p0-1)*(-2*S*p1*(p1-1))+p0*(p1-1)^2)/((p0-1)*(p0-1+(p1-1)^2)*log(p2)))
    d0   = sum(d0_1)+sum(d0_2)+sum(d0_3)

    d1_1 = event1*(der_h_Alpha/h+der_S_Alpha/S-Theta*der_S_Alpha-(der_S_Alpha*Theta*p1)/(p1-1)+Theta*der_ST_Alpha)
    d1_2 = event2*(der_h_Alpha/h+der_S_Alpha/S-Theta*der_S_Alpha-(der_S_Alpha*Theta*p1)/(p1-1)+Theta*der_ST_Alpha)
    d1_3 = (1-event1-event2)*(-2*Theta*der_S_Alpha*p1*(p1-1)/(log(p2)*(p0-1+(p1-1)^2)))
    d1   = sum(d1_1)+sum(d1_2)+sum(d1_3)

    d2_1 = event1*(der_h_Gamma/h+der_S_Gamma/S-Theta*der_S_Gamma-(der_S_Gamma*Theta*p1)/(p1-1)+Theta*der_ST_Gamma)
    d2_2 = event2*(der_h_Gamma/h+der_S_Gamma/S-Theta*der_S_Gamma-(der_S_Gamma*Theta*p1)/(p1-1)+Theta*der_ST_Gamma)
    d2_3 = (1-event1-event2)*(1-event1-event2)*(-2*Theta*der_S_Gamma*p1*(p1-1)/(log(p2)*(p0-1+(p1-1)^2)))
    d2   = sum(d2_1)+sum(d2_2)+sum(d2_3)

    d00_1 = event1*((S^2*p1*(p1-1)-S^2*p1^2)/(p1-1)^2+(-p0*(p0-1)+p0^2)/(p0-1)^2+2*der_ST_Theta+Theta*der_ST_Theta_Theta)
    d00_2 = event2*((S^2*p1*(p1-1)-S^2*p1^2)/(p1-1)^2+(-p0*(p0-1)+p0^2)/(p0-1)^2+2*der_ST_Theta+Theta*der_ST_Theta_Theta)
    d00_3 = (1-event1-event2)*(Theta^(-2)+((log(p2)*(p0-1+(p1-1)^2)*(p0-1))*(2*S^2*p1*(p0-1)*(2*p1-1)-p0*(p1-1)^2)-((p0-1)*(-2*S*p1*(p1-1))+p0*(p1-1)^2)*(p0*(p1-1)^2-2*S*(p1-1)*(p0-1)*p1+log(p2)*(-p0*(p0-1+(p1-1)^2)+(p0-1)*(-p0-2*S*p1*(p1-1)))))/(log(p2)*(p0-1+(p1-1)^2)*(p0-1))^2)
    d00   = sum(d00_1)+sum(d00_2)+sum(d00_3)

    d11_1 = event1*((h*der_h_Alpha_Alpha-(der_h_Alpha)^2)/(h)^2+(S*der_S_Alpha_Alpha-(der_S_Alpha)^2)/(S)^2-Theta*der_S_Alpha_Alpha-((Theta*der_S_Alpha_Alpha*p1-Theta^2*der_S_Alpha^2*p1)*(p1-1)+Theta^2*der_S_Alpha^2*p1^2)/(p1-1)^2+Theta*der_ST_Alpha_Alpha)
    d11_2 = event2*((h*der_h_Alpha_Alpha-(der_h_Alpha)^2)/(h)^2+(S*der_S_Alpha_Alpha-(der_S_Alpha)^2)/(S)^2-Theta*der_S_Alpha_Alpha-((Theta*der_S_Alpha_Alpha*p1-Theta^2*der_S_Alpha^2*p1)*(p1-1)+Theta^2*der_S_Alpha^2*p1^2)/(p1-1)^2+Theta*der_ST_Alpha_Alpha)
    d11_3 = (1-event1-event2)*((-2*Theta*p1*(der_S_Alpha_Alpha*(p1-1)+Theta*der_S_Alpha^2*(-2*p1+1))*(p0-1+(p1-1)^2)*log(p2)-(-2*Theta*der_S_Alpha*p1*(p1-1))*(log(p2)*(-2*(p1-1)*Theta*der_S_Alpha*p1)-2*Theta*der_S_Alpha*p1*(p1-1)))/((p0-1+(p1-1)^2)*log(p2))^2)
    d11   = sum(d11_1)+sum(d11_2)+sum(d11_3)

    d22_1 = event1*((h*der_h_Gamma_Gamma-(der_h_Gamma)^2)/(h)^2+(S*der_S_Gamma_Gamma-(der_S_Gamma)^2)/(S)^2-Theta*der_S_Gamma_Gamma-((Theta*der_S_Gamma_Gamma*p1-Theta^2*der_S_Gamma^2*p1)*(p1-1)+Theta^2*der_S_Gamma^2*p1^2)/(p1-1)^2+Theta*der_ST_Gamma_Gamma)
    d22_2 = event2*((h*der_h_Gamma_Gamma-(der_h_Gamma)^2)/(h)^2+(S*der_S_Gamma_Gamma-(der_S_Gamma)^2)/(S)^2-Theta*der_S_Gamma_Gamma-((Theta*der_S_Gamma_Gamma*p1-Theta^2*der_S_Gamma^2*p1)*(p1-1)+Theta^2*der_S_Gamma^2*p1^2)/(p1-1)^2+Theta*der_ST_Gamma_Gamma)
    d22_3 = (1-event1-event2)*((-2*Theta*p1*(der_S_Gamma_Gamma*(p1-1)+Theta*der_S_Gamma^2*(-2*p1+1))*(p0-1+(p1-1)^2)*log(p2)-(-2*Theta*der_S_Gamma*p1*(p1-1))*(log(p2)*(-2*(p1-1)*Theta*der_S_Gamma*p1)-2*Theta*der_S_Gamma*p1*(p1-1)))/((p0-1+(p1-1)^2)*log(p2))^2)
    d22   = sum(d22_1)+sum(d22_2)+sum(d22_3)

    d12_1 = event1*((der_h_Alpha_Gamma*h-der_h_Alpha*der_h_Gamma)/h^2+(der_S_Alpha_Gamma*S-der_S_Alpha*der_S_Gamma)/S^2-Theta*der_S_Alpha_Gamma-(Theta*(p1-1)*p1*(der_S_Alpha_Gamma-Theta*der_S_Alpha*der_S_Gamma)+Theta^2*der_S_Alpha*der_S_Gamma*p1^2)/(p1-1)^2+Theta*der_ST_Alpha_Gamma)
    d12_2 = event2*((der_h_Alpha_Gamma*h-der_h_Alpha*der_h_Gamma)/h^2+(der_S_Alpha_Gamma*S-der_S_Alpha*der_S_Gamma)/S^2-Theta*der_S_Alpha_Gamma-(Theta*(p1-1)*p1*(der_S_Alpha_Gamma-Theta*der_S_Alpha*der_S_Gamma)+Theta^2*der_S_Alpha*der_S_Gamma*p1^2)/(p1-1)^2+Theta*der_ST_Alpha_Gamma)
    d12_3 = (1-event1-event2)*((log(p2)*(p0-1+(p1-1)^2)*(-2*Theta*der_S_Alpha_Gamma*p1*(p1-1)-2*Theta^2*der_S_Alpha*der_S_Gamma*p1*(1-2*p1))-(-2*Theta*der_S_Alpha*p1*(p1-1))*(-2*Theta*der_S_Gamma*p1*(p1-1))*(1+log(p2)))/(log(p2)*(p0-1+(p1-1)^2))^2)
    d12   = sum(d12_1)+sum(d12_2)+sum(d12_3)

    d01_1 = event1*(-der_S_Alpha-((der_S_Alpha*p1-Theta*S*der_S_Alpha*p1)*(p1-1)+Theta*der_S_Alpha*S*p1^2)/(p1-1)^2+der_ST_Alpha+Theta*der_ST_Theta_Alpha)
    d01_2 = event2*(-der_S_Alpha-((der_S_Alpha*p1-Theta*S*der_S_Alpha*p1)*(p1-1)+Theta*der_S_Alpha*S*p1^2)/(p1-1)^2+der_ST_Alpha+Theta*der_ST_Theta_Alpha)
    d01_3 = (1-event1-event2)*((log(p2)*(p0-1+(p1-1)^2)*(-2*der_S_Alpha*p1*(p1-1)-2*Theta*der_S_Alpha*S*p1*(1-2*p1))-(-2*Theta*der_S_Alpha*p1*(p1-1))*((-2*S*p1*(p1-1)*(p0-1)+p0*(p1-1)^2)/(p0-1)+log(p2)*(-p0-2*S*p1*(p1-1))))/(log(p2)*(p0-1+(p1-1)^2))^2)
    d01   = sum(d01_1)+sum(d01_2)+sum(d01_3)

    d02_1 = event1*(-der_S_Gamma-((der_S_Gamma*p1-Theta*S*der_S_Gamma*p1)*(p1-1)+Theta*der_S_Gamma*S*p1^2)/(p1-1)^2+der_ST_Gamma+Theta*der_ST_Theta_Gamma)
    d02_2 = event2*(-der_S_Gamma-((der_S_Gamma*p1-Theta*S*der_S_Gamma*p1)*(p1-1)+Theta*der_S_Gamma*S*p1^2)/(p1-1)^2+der_ST_Gamma+Theta*der_ST_Theta_Gamma)
    d02_3 = (1-event1-event2)*((log(p2)*(p0-1+(p1-1)^2)*(-2*der_S_Gamma*p1*(p1-1)-2*Theta*der_S_Gamma*S*p1*(1-2*p1))-(-2*Theta*der_S_Gamma*p1*(p1-1))*((-2*S*p1*(p1-1)*(p0-1)+p0*(p1-1)^2)/(p0-1)+log(p2)*(-p0-2*S*p1*(p1-1))))/(log(p2)*(p0-1+(p1-1)^2))^2)
    d02   = sum(d02_1)+sum(d02_2)+sum(d02_3)

    M00 = d00
    M01 = d01*exp(par[2])
    M02 = d02*exp(par[3])
    M11 = d11*exp(2*par[2])+D1*exp(par[2])
    M12 = d12*exp(par[2]+par[3])
    M22 = d22*exp(2*par[3])+D2*exp(par[3])

    matrix(c(M00,M01,M02,M01,M11,M12,M02,M12,M22),3,3)

  }
  H_function = function (par) {

    Theta = par[1]
    Alpha = par[2]
    Gamma = par[3]

    h = Alpha*Gamma/(1+Alpha*t.event)
    S = (1+Alpha*t.event)^(-Gamma)
    p0 = exp(-Theta)
    p1 = exp(-Theta*S)
    p2 = (1 + (p1 - 1)^2/(p0 - 1))
    ST = -(1/Theta)*log(1+(exp(-Theta*S)-1)*(exp(-Theta*S)-1)/(exp(-Theta)-1))

    der_h_Alpha = Gamma/(1+Alpha*t.event)^2
    der_S_Alpha = -Gamma*t.event*(1+Alpha*t.event)^(-Gamma-1)
    der_h_Gamma = Alpha/(1+Alpha*t.event)
    der_S_Gamma = -(1+Alpha*t.event)^(-Gamma)*log(1+Alpha*t.event)

    der_ST_Theta = (Theta^(-2)*log(1+(p1-1)^2/(p0-1))-(p1-1)*(-2*S*p1*(p0-1)+p0*(p1-1))/(Theta*(p0-1)*((p1-1)^2+p0-1)))
    der_ST_Alpha = (p1-1)*2*der_S_Alpha*p1/(p0-1+(p1-1)^2)
    der_ST_Gamma = (p1-1)*2*der_S_Gamma*p1/(p0-1+(p1-1)^2)

    der_h_Alpha_Alpha = -2*Gamma*t.event/(1+Alpha*t.event)^3
    der_S_Alpha_Alpha = Gamma*(Gamma+1)*t.event^2*(1+Alpha*t.event)^(-Gamma-2)
    der_h_Gamma_Gamma = 0
    der_S_Gamma_Gamma = (1+Alpha*t.event)^(-Gamma)*(log(1+Alpha*t.event))^2
    der_h_Alpha_Gamma = (1+Alpha*t.event)^(-2)
    der_S_Alpha_Gamma = t.event*(Gamma*log(1+Alpha*t.event)-1)/(1+Alpha*t.event)^(Gamma+1)

    der_ST_Theta_Theta = -2*Theta^(-3)*log(p2)+(p0*(p1-1)^2-2*S*p1*(p1-1)*(p0-1))/(Theta^2*(p0-1)*(p0-1+(p1-1)^2))-((Theta*(p0-1)*(p0-1+(p1-1)^2))*(2*S^2*p1*(p0-1)*(2*p1-1)-p0*(p1-1)^2)-(-2*S*p1*(p1-1)*(p0-1)+p0*(p1-1)^2)*((p0-1)*(p0-1+(p1-1)^2)+Theta*(-p0*(p0-1+(p1-1)^2)+(p0-1)*(-p0-2*S*p1*(p1-1)))))/(Theta*(p0-1)*(p0-1+(p1-1)^2))^2
    der_ST_Alpha_Alpha = ((p0-1+(p1-1)^2)*(2*der_S_Alpha_Alpha*p1*(p1-1)+2*der_S_Alpha*(-2*Theta*der_S_Alpha*p1^2+Theta*der_S_Alpha*p1))-2*der_S_Alpha*(p1^2-p1)*(-2*Theta*der_S_Alpha*p1*(p1-1)))/(p0-1+(p1-1)^2)^2
    der_ST_Gamma_Gamma = ((p0-1+(p1-1)^2)*(2*der_S_Gamma_Gamma*p1*(p1-1)+2*der_S_Gamma*(-2*Theta*der_S_Gamma*p1^2+Theta*der_S_Gamma*p1))-2*der_S_Gamma*(p1^2-p1)*(-2*Theta*der_S_Gamma*p1*(p1-1)))/(p0-1+(p1-1)^2)^2
    der_ST_Theta_Alpha = ((p0-1+(p1-1)^2)*(2*der_S_Alpha*S*p1*(1-2*p1))-(2*der_S_Alpha*p1*(p1-1))*(-p0-2*S*p1*(p1-1)))/(p0-1+(p1-1)^2)^2
    der_ST_Theta_Gamma = ((p0-1+(p1-1)^2)*(2*der_S_Gamma*S*p1*(1-2*p1))-(2*der_S_Gamma*p1*(p1-1))*(-p0-2*S*p1*(p1-1)))/(p0-1+(p1-1)^2)^2
    der_ST_Alpha_Gamma = ((p0-1+(p1-1)^2)*(2*der_S_Alpha_Gamma*p1*(p1-1)+2*der_S_Alpha*(-2*Theta*der_S_Gamma*p1^2+Theta*der_S_Gamma*p1))-2*der_S_Alpha*(p1^2-p1)*(-2*Theta*der_S_Gamma*p1*(p1-1)))/(p0-1+(p1-1)^2)^2

    d0_1 = event1*(-S-S*p1/(p1-1)+exp(-Theta)/(exp(-Theta)-1)+ST+Theta*der_ST_Theta)
    d0_2 = event2*(-S-S*p1/(p1-1)+exp(-Theta)/(exp(-Theta)-1)+ST+Theta*der_ST_Theta)
    d0_3 = (1-event1-event2)*((-1/Theta)+((p0-1)*(-2*S*p1*(p1-1))+p0*(p1-1)^2)/((p0-1)*(p0-1+(p1-1)^2)*log(p2)))
    d0   = sum(d0_1)+sum(d0_2)+sum(d0_3)

    d1_1 = event1*(der_h_Alpha/h+der_S_Alpha/S-Theta*der_S_Alpha-(der_S_Alpha*Theta*p1)/(p1-1)+Theta*der_ST_Alpha)
    d1_2 = event2*(der_h_Alpha/h+der_S_Alpha/S-Theta*der_S_Alpha-(der_S_Alpha*Theta*p1)/(p1-1)+Theta*der_ST_Alpha)
    d1_3 = (1-event1-event2)*(-2*Theta*der_S_Alpha*p1*(p1-1)/(log(p2)*(p0-1+(p1-1)^2)))
    d1   = sum(d1_1)+sum(d1_2)+sum(d1_3)

    d2_1 = event1*(der_h_Gamma/h+der_S_Gamma/S-Theta*der_S_Gamma-(der_S_Gamma*Theta*p1)/(p1-1)+Theta*der_ST_Gamma)
    d2_2 = event2*(der_h_Gamma/h+der_S_Gamma/S-Theta*der_S_Gamma-(der_S_Gamma*Theta*p1)/(p1-1)+Theta*der_ST_Gamma)
    d2_3 = (1-event1-event2)*(1-event1-event2)*(-2*Theta*der_S_Gamma*p1*(p1-1)/(log(p2)*(p0-1+(p1-1)^2)))
    d2   = sum(d2_1)+sum(d2_2)+sum(d2_3)

    d00_1 = event1*((S^2*p1*(p1-1)-S^2*p1^2)/(p1-1)^2+(-p0*(p0-1)+p0^2)/(p0-1)^2+2*der_ST_Theta+Theta*der_ST_Theta_Theta)
    d00_2 = event2*((S^2*p1*(p1-1)-S^2*p1^2)/(p1-1)^2+(-p0*(p0-1)+p0^2)/(p0-1)^2+2*der_ST_Theta+Theta*der_ST_Theta_Theta)
    d00_3 = (1-event1-event2)*(Theta^(-2)+((log(p2)*(p0-1+(p1-1)^2)*(p0-1))*(2*S^2*p1*(p0-1)*(2*p1-1)-p0*(p1-1)^2)-((p0-1)*(-2*S*p1*(p1-1))+p0*(p1-1)^2)*(p0*(p1-1)^2-2*S*(p1-1)*(p0-1)*p1+log(p2)*(-p0*(p0-1+(p1-1)^2)+(p0-1)*(-p0-2*S*p1*(p1-1)))))/(log(p2)*(p0-1+(p1-1)^2)*(p0-1))^2)
    d00   = sum(d00_1)+sum(d00_2)+sum(d00_3)

    d11_1 = event1*((h*der_h_Alpha_Alpha-(der_h_Alpha)^2)/(h)^2+(S*der_S_Alpha_Alpha-(der_S_Alpha)^2)/(S)^2-Theta*der_S_Alpha_Alpha-((Theta*der_S_Alpha_Alpha*p1-Theta^2*der_S_Alpha^2*p1)*(p1-1)+Theta^2*der_S_Alpha^2*p1^2)/(p1-1)^2+Theta*der_ST_Alpha_Alpha)
    d11_2 = event2*((h*der_h_Alpha_Alpha-(der_h_Alpha)^2)/(h)^2+(S*der_S_Alpha_Alpha-(der_S_Alpha)^2)/(S)^2-Theta*der_S_Alpha_Alpha-((Theta*der_S_Alpha_Alpha*p1-Theta^2*der_S_Alpha^2*p1)*(p1-1)+Theta^2*der_S_Alpha^2*p1^2)/(p1-1)^2+Theta*der_ST_Alpha_Alpha)
    d11_3 = (1-event1-event2)*((-2*Theta*p1*(der_S_Alpha_Alpha*(p1-1)+Theta*der_S_Alpha^2*(-2*p1+1))*(p0-1+(p1-1)^2)*log(p2)-(-2*Theta*der_S_Alpha*p1*(p1-1))*(log(p2)*(-2*(p1-1)*Theta*der_S_Alpha*p1)-2*Theta*der_S_Alpha*p1*(p1-1)))/((p0-1+(p1-1)^2)*log(p2))^2)
    d11   = sum(d11_1)+sum(d11_2)+sum(d11_3)

    d22_1 = event1*((h*der_h_Gamma_Gamma-(der_h_Gamma)^2)/(h)^2+(S*der_S_Gamma_Gamma-(der_S_Gamma)^2)/(S)^2-Theta*der_S_Gamma_Gamma-((Theta*der_S_Gamma_Gamma*p1-Theta^2*der_S_Gamma^2*p1)*(p1-1)+Theta^2*der_S_Gamma^2*p1^2)/(p1-1)^2+Theta*der_ST_Gamma_Gamma)
    d22_2 = event2*((h*der_h_Gamma_Gamma-(der_h_Gamma)^2)/(h)^2+(S*der_S_Gamma_Gamma-(der_S_Gamma)^2)/(S)^2-Theta*der_S_Gamma_Gamma-((Theta*der_S_Gamma_Gamma*p1-Theta^2*der_S_Gamma^2*p1)*(p1-1)+Theta^2*der_S_Gamma^2*p1^2)/(p1-1)^2+Theta*der_ST_Gamma_Gamma)
    d22_3 = (1-event1-event2)*((-2*Theta*p1*(der_S_Gamma_Gamma*(p1-1)+Theta*der_S_Gamma^2*(-2*p1+1))*(p0-1+(p1-1)^2)*log(p2)-(-2*Theta*der_S_Gamma*p1*(p1-1))*(log(p2)*(-2*(p1-1)*Theta*der_S_Gamma*p1)-2*Theta*der_S_Gamma*p1*(p1-1)))/((p0-1+(p1-1)^2)*log(p2))^2)
    d22   = sum(d22_1)+sum(d22_2)+sum(d22_3)

    d12_1 = event1*((der_h_Alpha_Gamma*h-der_h_Alpha*der_h_Gamma)/h^2+(der_S_Alpha_Gamma*S-der_S_Alpha*der_S_Gamma)/S^2-Theta*der_S_Alpha_Gamma-(Theta*(p1-1)*p1*(der_S_Alpha_Gamma-Theta*der_S_Alpha*der_S_Gamma)+Theta^2*der_S_Alpha*der_S_Gamma*p1^2)/(p1-1)^2+Theta*der_ST_Alpha_Gamma)
    d12_2 = event2*((der_h_Alpha_Gamma*h-der_h_Alpha*der_h_Gamma)/h^2+(der_S_Alpha_Gamma*S-der_S_Alpha*der_S_Gamma)/S^2-Theta*der_S_Alpha_Gamma-(Theta*(p1-1)*p1*(der_S_Alpha_Gamma-Theta*der_S_Alpha*der_S_Gamma)+Theta^2*der_S_Alpha*der_S_Gamma*p1^2)/(p1-1)^2+Theta*der_ST_Alpha_Gamma)
    d12_3 = (1-event1-event2)*((log(p2)*(p0-1+(p1-1)^2)*(-2*Theta*der_S_Alpha_Gamma*p1*(p1-1)-2*Theta^2*der_S_Alpha*der_S_Gamma*p1*(1-2*p1))-(-2*Theta*der_S_Alpha*p1*(p1-1))*(-2*Theta*der_S_Gamma*p1*(p1-1))*(1+log(p2)))/(log(p2)*(p0-1+(p1-1)^2))^2)
    d12   = sum(d12_1)+sum(d12_2)+sum(d12_3)

    d01_1 = event1*(-der_S_Alpha-((der_S_Alpha*p1-Theta*S*der_S_Alpha*p1)*(p1-1)+Theta*der_S_Alpha*S*p1^2)/(p1-1)^2+der_ST_Alpha+Theta*der_ST_Theta_Alpha)
    d01_2 = event2*(-der_S_Alpha-((der_S_Alpha*p1-Theta*S*der_S_Alpha*p1)*(p1-1)+Theta*der_S_Alpha*S*p1^2)/(p1-1)^2+der_ST_Alpha+Theta*der_ST_Theta_Alpha)
    d01_3 = (1-event1-event2)*((log(p2)*(p0-1+(p1-1)^2)*(-2*der_S_Alpha*p1*(p1-1)-2*Theta*der_S_Alpha*S*p1*(1-2*p1))-(-2*Theta*der_S_Alpha*p1*(p1-1))*((-2*S*p1*(p1-1)*(p0-1)+p0*(p1-1)^2)/(p0-1)+log(p2)*(-p0-2*S*p1*(p1-1))))/(log(p2)*(p0-1+(p1-1)^2))^2)
    d01   = sum(d01_1)+sum(d01_2)+sum(d01_3)

    d02_1 = event1*(-der_S_Gamma-((der_S_Gamma*p1-Theta*S*der_S_Gamma*p1)*(p1-1)+Theta*der_S_Gamma*S*p1^2)/(p1-1)^2+der_ST_Gamma+Theta*der_ST_Theta_Gamma)
    d02_2 = event2*(-der_S_Gamma-((der_S_Gamma*p1-Theta*S*der_S_Gamma*p1)*(p1-1)+Theta*der_S_Gamma*S*p1^2)/(p1-1)^2+der_ST_Gamma+Theta*der_ST_Theta_Gamma)
    d02_3 = (1-event1-event2)*((log(p2)*(p0-1+(p1-1)^2)*(-2*der_S_Gamma*p1*(p1-1)-2*Theta*der_S_Gamma*S*p1*(1-2*p1))-(-2*Theta*der_S_Gamma*p1*(p1-1))*((-2*S*p1*(p1-1)*(p0-1)+p0*(p1-1)^2)/(p0-1)+log(p2)*(-p0-2*S*p1*(p1-1))))/(log(p2)*(p0-1+(p1-1)^2))^2)
    d02   = sum(d02_1)+sum(d02_2)+sum(d02_3)

    matrix(c(d00,d01,d02,d01,d11,d12,d02,d12,d22),3,3)

  }

  par_old = c(Theta.0,log(Alpha.0),log(Gamma.0))
  count  = 0
  random = 0

  repeat{

    temp1 = try(solve(HL_function(par_old)),silent = TRUE)
    if (class(temp1) == "try-error") {

      random = random+1
      count = 0
      par_old = c(Theta.0+runif(1,-r.1,r.1),
                  log(Alpha.0*exp(runif(1,-r.2,r.2))),
                  log(Gamma.0*exp(runif(1,-r.3,r.3))))
      next

    }

    par_new = par_old-solve(HL_function(par_old))%*%SL_function(par_old)

    if (is.na(sum(par_new)) | count > 20) {

      random = random+1
      count = 0
      par_old = c(Theta.0+runif(1,-r.1,r.1),
                  log(Alpha.0*exp(runif(1,-r.2,r.2))),
                  log(Gamma.0*exp(runif(1,-r.3,r.3))))
      next

    }
    count = count+1

    if (max(abs(c(par_old[1],exp(par_old[2]),exp(par_old[3]))
               -c(par_new[1],exp(par_new[2]),exp(par_new[3])))) < epsilon) {

      temp2 = try(solve(HL_function(par_new)),silent = TRUE)
      temp3 = try(solve(H_function(c(par_new[1],exp(par_new[2]),exp(par_new[3])))),silent = TRUE)
      if (class(temp2) == "try-error" |
          class(temp3) == "try-error" |
          det(HL_function(par_new)) > 0) {

        random = random+1
        count = 0
        par_old = c(Theta.0+runif(1,-r.1,r.1),
                    log(Alpha.0*exp(runif(1,-r.2,r.2))),
                    log(Gamma.0*exp(runif(1,-r.3,r.3))))
        next

      }
      break

    }
    par_old = par_new

  }

  Theta_hat = par_new[1]
  Alpha_hat = exp(par_new[2])
  Gamma_hat = exp(par_new[3])

  Info = solve(-H_function(c(par_new[1],exp(par_new[2]),exp(par_new[3]))))
  Theta_se = sqrt(Info[1,1])
  Alpha_se = sqrt(Info[2,2])
  Gamma_se = sqrt(Info[3,3])

  InfoL = solve(-HL_function(par_new))
  CI_Theta = c(Theta_hat-qnorm(0.975)*sqrt(InfoL[1,1]),
               Theta_hat+qnorm(0.975)*sqrt(InfoL[1,1]))
  CI_Alpha = c(Alpha_hat*exp(-qnorm(0.975)*sqrt(InfoL[2,2])),
               Alpha_hat*exp(+qnorm(0.975)*sqrt(InfoL[2,2])))
  CI_Gamma = c(Gamma_hat*exp(-qnorm(0.975)*sqrt(InfoL[3,3])),
               Gamma_hat*exp(+qnorm(0.975)*sqrt(InfoL[3,3])))

  Theta.res = c(Estimate = Theta_hat,SE = Theta_se,CI.lower = CI_Theta[1],CI.upper = CI_Theta[2])
  Alpha.res = c(Estimate = Alpha_hat,SE = Alpha_se,CI.lower = CI_Alpha[1],CI.upper = CI_Alpha[2])
  Gamma.res = c(Estimate = Gamma_hat,SE = Gamma_se,CI.lower = CI_Gamma[1],CI.upper = CI_Gamma[2])

  if (bootstrap) {

    Theta_B = rep(0,B)
    Alpha_B = rep(0,B)
    Gamma_B = rep(0,B)
    cen_max = max(t.event*(1-event1-event2))

    for (b in 1:B) {

      data_B = Frank.Pareto(n,Theta_hat,Alpha_hat,Alpha_hat,Gamma_hat,Gamma_hat)
      if (cen_max == 0) {

        t.event = pmin(data_B[,1],data_B[,2])
        event1 = rep(1,n)*(data_B[,1] < data_B[,2])
        event2 = rep(1,n)*(data_B[,2] < data_B[,1])

      } else {

        cen_B = runif(n,0,cen_max)
        t.event = pmin(data_B[,1],data_B[,2],cen_B)
        event1 = rep(1,n)*(data_B[,1] < data_B[,2] & data_B[,1] < cen_B)
        event2 = rep(1,n)*(data_B[,2] < data_B[,1] & data_B[,2] < cen_B)

      }

      par_old_B = c(Theta.0,log(Alpha.0),log(Gamma.0))
      count  = 0
      random = 0

      repeat{

        temp1 = try(solve(HL_function(par_old_B)),silent = TRUE)
        if (class(temp1) == "try-error") {

          random = random+1
          count = 0
          par_old_B = c(Theta.0+runif(1,-r.1,r.1),
                        log(Alpha.0*exp(runif(1,-r.2,r.2))),
                        log(Gamma.0*exp(runif(1,-r.3,r.3))))
          next

        }

        par_new_B = par_old_B-solve(HL_function(par_old_B))%*%SL_function(par_old_B)

        if (is.na(sum(par_new_B)) | count > 20) {

          random = random+1
          count = 0
          par_old_B = c(Theta.0+runif(1,-r.1,r.1),
                        log(Alpha.0*exp(runif(1,-r.2,r.2))),
                        log(Gamma.0*exp(runif(1,-r.3,r.3))))
          next

        }
        count = count+1

        if (max(abs(c(par_old_B[1],exp(par_old_B[2]),exp(par_old_B[3]))
                   -c(par_new_B[1],exp(par_new_B[2]),exp(par_new_B[3])))) < epsilon) {

          temp2 = try(solve(HL_function(par_new_B)),silent = TRUE)
          temp3 = try(solve(H_function(c(par_new_B[1],exp(par_new_B[2]),exp(par_new_B[3])))),silent = TRUE)
          if (class(temp2) == "try-error" |
              class(temp3) == "try-error" |
              det(HL_function(par_new_B)) > 0) {

            random = random+1
            count = 0
            par_old_B = c(Theta.0+runif(1,-r.1,r.1),
                          log(Alpha.0*exp(runif(1,-r.2,r.2))),
                          log(Gamma.0*exp(runif(1,-r.3,r.3))))
            next

          }
          break

        }
        par_old_B = par_new_B

      }
      Theta_B[b] = par_new_B[1]
      Alpha_B[b] = exp(par_new_B[2])
      Gamma_B[b] = exp(par_new_B[3])

    }

    Theta_se_B = sd(Theta_B)
    Alpha_se_B = sd(Alpha_B)
    Gamma_se_B = sd(Gamma_B)

    CI_Theta_lower_B = sort(Theta_B)[B*0.025]
    CI_Alpha_lower_B = sort(Alpha_B)[B*0.025]
    CI_Gamma_lower_B = sort(Gamma_B)[B*0.025]

    CI_Theta_upper_B = sort(Theta_B)[B*0.975]
    CI_Alpha_upper_B = sort(Alpha_B)[B*0.975]
    CI_Gamma_upper_B = sort(Gamma_B)[B*0.975]

    Theta.res_B = c(Estimate = Theta_hat,SE = Theta_se_B,
                    CI.lower = CI_Theta_lower_B,CI.upper = CI_Theta_upper_B)
    Alpha.res_B = c(Estimate = Alpha_hat,SE = Alpha_se_B,
                    CI.lower = CI_Alpha_lower_B,CI.upper = CI_Alpha_upper_B)
    Gamma.res_B = c(Estimate = Gamma_hat,SE = Gamma_se_B,
                    CI.lower = CI_Gamma_lower_B,CI.upper = CI_Gamma_upper_B)

    return(list(n = n,Iteration = count,Randomization = random,
                Theta = Theta.res,Theta.B = Theta.res_B,
                Alpha = Alpha.res,Alpha.B = Alpha.res_B,
                Gamma = Gamma.res,Gamma.B = Gamma.res_B,
                logL = log_L(par_new),AIC = 2*length(par_new)-2*log_L(par_new),
                BIC = length(par_new)*log(length(t.event))-2*log_L(par_new)))

  } else {

    return(list(n = n,Iteration = count,Randomization = random,
                Theta = Theta.res,Alpha = Alpha.res,Gamma = Gamma.res,
                logL = log_L(par_new),AIC = 2*length(par_new)-2*log_L(par_new),
                BIC = length(par_new)*log(length(t.event))-2*log_L(par_new)))

  }

}
