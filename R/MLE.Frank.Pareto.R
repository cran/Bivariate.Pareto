#' Maximum likelihood estimation for bivariate dependent competing risks data under the Frank copula with the Pareto margins and fixed \eqn{\theta}
#'
#' @param t.event Vector of the observed failure times.
#' @param event1 Vector of the indicators for the failure cause 1.
#' @param event2 Vector of the indicators for the failure cause 2.
#' @param Theta Copula parameter \eqn{\theta}.
#' @param Alpha1.0 Initial guess for the scale parameter \eqn{\alpha_{1}} with default value 1.
#' @param Alpha2.0 Initial guess for the scale parameter \eqn{\alpha_{2}} with default value 1.
#' @param Gamma1.0 Initial guess for the shape parameter \eqn{\gamma_{1}} with default value 1.
#' @param Gamma2.0 Initial guess for the shape parameter \eqn{\gamma_{2}} with default value 1.
#' @param epsilon Positive tunning parameter in the NR algorithm with default value \eqn{10^{-5}}.
#' @param d Positive tunning parameter in the NR algorithm with default value \eqn{e^{10}}.
#' @param r.1 Positive tunning parameter in the NR algorithm with default value 1.
#' @param r.2 Positive tunning parameter in the NR algorithm with default value 1.
#' @param r.3 Positive tunning parameter in the NR algorithm with default value 1.
#' @param r.4 Positive tunning parameter in the NR algorithm with default value 1.
#' @description Maximum likelihood estimation for bivariate dependent competing risks data under the Frank copula with the Pareto margins and fixed \eqn{\theta}.
#'
#' @return \item{n}{Sample size.}
#' \item{count}{Iteration number.}
#' \item{random}{Randomization number.}
#' \item{Alpha1}{Positive scale parameter for the Pareto margin (failure cause 1).}
#' \item{Alpha2}{Positive scale parameter for the Pareto margin (failure cause 2).}
#' \item{Gamma1}{Positive shape parameter for the Pareto margin (failure cause 1).}
#' \item{Gamma2}{Positive shape parameter for the Pareto margin (failure cause 2).}
#' \item{MedX}{Median lifetime due to failure cause 1.}
#' \item{MedY}{Median lifetime due to failure cause 2.}
#' \item{MeanX}{Mean lifetime due to failure cause 1.}
#' \item{MeanY}{Mean lifetime due to failure cause 2.}
#' \item{logL}{Log-likelihood value under the fitted model.}
#' \item{AIC}{AIC value under the fitted model.}
#' \item{BIC}{BIC value under the fitted model.}
#'
#' @references Shih J-H, Lee W, Sun L-H, Emura T (2018), Fitting competing risks data to bivariate Pareto models, Communications in Statistics - Theory and Methods, doi: 10.1080/03610926.2018.1425450.
#' @importFrom stats qnorm runif
#' @importFrom utils globalVariables
#' @importFrom methods is
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
#' MLE.Frank.Pareto(t.event,event1,event2,Theta = -5)

MLE.Frank.Pareto = function(t.event,event1,event2,Theta,Alpha1.0 = 1,Alpha2.0 = 1,
                            Gamma1.0 = 1,Gamma2.0 = 1,epsilon = 1e-5,d = exp(10),
                            r.1 = 6,r.2 = 6,r.3 = 6,r.4 = 6) {

  ### checking inputs ###
  n = length(t.event)
  if (length(t.event[t.event < 0]) != 0) {stop("t.event must be non-negative")}
  if (length(event1) != n) {stop("the length of event1 is different from t.event")}
  if (length(event2) != n) {stop("the length of event2 is different from t.event")}
  if (length(event1[event1 == 0 | event1 == 1]) != n) {stop("elements in event1 must be either 0 or 1")}
  if (length(event2[event2 == 0 | event2 == 1]) != n) {stop("elements in event2 must be either 0 or 1")}
  temp.event = event1+event2
  if (length(temp.event[temp.event == 2]) != 0) {stop("event1 and event2 cannot be 1 simultaneously")}
  if (Theta == 0) {stop("Theta cannot be zero")}
  if (Alpha1.0 <= 0) {stop("Alpha1.0 must be positive")}
  if (Alpha2.0 <= 0) {stop("Alpha2.0 must be positive")}
  if (Gamma1.0 <= 0) {stop("Alpha1.0 must be positive")}
  if (Gamma2.0 <= 0) {stop("Alpha2.0 must be positive")}
  if (epsilon <= 0) {stop("epsilon must be positive")}
  if (d <= 0) {stop("d must be positive")}
  if (r.1 <= 0) {stop("r.1 must be positive")}
  if (r.2 <= 0) {stop("r.2 must be positive")}
  if (r.3 <= 0) {stop("r.3 must be positive")}
  if (r.4 <= 0) {stop("r.3 must be positive")}

  ### functions ###
  log_L = function(par){

    Alpha1 = exp(par[1])
    Alpha2 = exp(par[2])
    Gamma1 = exp(par[3])
    Gamma2 = exp(par[4])

    h1 = Alpha1*Gamma1/(1+Alpha1*t.event)
    h2 = Alpha2*Gamma2/(1+Alpha2*t.event)
    S1 = (1+Alpha1*t.event)^(-Gamma1)
    S2 = (1+Alpha2*t.event)^(-Gamma2)
    ST = -(1/Theta)*log(1+(exp(-Theta*S1)-1)*(exp(-Theta*S2)-1)/(exp(-Theta)-1))
    f1 = h1*S1*exp(-Theta*S1)*(exp(-Theta*S2)-1)/((exp(-Theta)-1)*exp(-Theta*ST))
    f2 = h2*S2*exp(-Theta*S2)*(exp(-Theta*S1)-1)/((exp(-Theta)-1)*exp(-Theta*ST))

    sum((1-event1-event2)*log(ST))+sum(event1*log(f1))+sum(event2*log(f2))

  }
  SL_function = function(par){

    Alpha1 = exp(par[1])
    Alpha2 = exp(par[2])
    Gamma1 = exp(par[3])
    Gamma2 = exp(par[4])

    h1 = Alpha1*Gamma1/(1+Alpha1*t.event)
    h2 = Alpha2*Gamma2/(1+Alpha2*t.event)
    S1 = (1+Alpha1*t.event)^(-Gamma1)
    S2 = (1+Alpha2*t.event)^(-Gamma2)
    p0 = exp(-Theta)
    p1 = exp(-Theta*S1)
    p2 = exp(-Theta*S2)

    ST = -(1/Theta)*log(1+(exp(-Theta*S1)-1)*(exp(-Theta*S2)-1)/(exp(-Theta)-1))

    der_h1_Alpha1 = Gamma1/(1+Alpha1*t.event)^2
    der_h2_Alpha2 = Gamma2/(1+Alpha2*t.event)^2

    der_S1_Alpha1 = -Gamma1*t.event*(1+Alpha1*t.event)^(-Gamma1-1)
    der_S2_Alpha2 = -Gamma2*t.event*(1+Alpha2*t.event)^(-Gamma2-1)

    der_h1_Gamma1 = Alpha1/(1+Alpha1*t.event)
    der_h2_Gamma2 = Alpha2/(1+Alpha2*t.event)
    der_S1_Gamma1 = -(1+Alpha1*t.event)^(-Gamma1)*log(1+Alpha1*t.event)
    der_S2_Gamma2 = -(1+Alpha2*t.event)^(-Gamma2)*log(1+Alpha2*t.event)

    der_ST_Alpha1 = der_S1_Alpha1*p1*(p2-1)/(p0-1+(p1-1)*(p2-1))
    der_ST_Alpha2 = der_S2_Alpha2*p2*(p1-1)/(p0-1+(p1-1)*(p2-1))

    der_ST_Gamma1 = der_S1_Gamma1*p1*(p2-1)/(p0-1+(p1-1)*(p2-1))
    der_ST_Gamma2 = der_S2_Gamma2*p2*(p1-1)/(p0-1+(p1-1)*(p2-1))

    d11 = sum(event1*(der_h1_Alpha1/h1+der_S1_Alpha1/S1-Theta*der_S1_Alpha1+Theta*der_ST_Alpha1))
    d12 = sum(event2*(-Theta*der_S1_Alpha1*p1/(p1-1)+Theta*der_ST_Alpha1))
    d13 = sum((1-event1-event2)*(der_ST_Alpha1/ST))
    d1  = d11+d12+d13

    d21 = sum(event1*(-Theta*der_S2_Alpha2*p2/(p2-1)+Theta*der_ST_Alpha2))
    d22 = sum(event2*(der_h2_Alpha2/h2+der_S2_Alpha2/S2-Theta*der_S2_Alpha2+Theta*der_ST_Alpha2))
    d23 = sum((1-event1-event2)*(der_ST_Alpha2/ST))
    d2  = d21+d22+d23

    d31 = sum(event1*(der_h1_Gamma1/h1+der_S1_Gamma1/S1-Theta*der_S1_Gamma1+Theta*der_ST_Gamma1))
    d32 = sum(event2*(-Theta*der_S1_Gamma1*p1/(p1-1)+Theta*der_ST_Gamma1))
    d33 = sum((1-event1-event2)*(der_ST_Gamma1/ST))
    d3  = d31+d32+d33

    d41 = sum(event1*(-Theta*der_S2_Gamma2*p2/(p2-1)+Theta*der_ST_Gamma2))
    d42 = sum(event2*(der_h2_Gamma2/h2+der_S2_Gamma2/S2-Theta*der_S2_Gamma2+Theta*der_ST_Gamma2))
    d43 = sum((1-event1-event2)*(der_ST_Gamma2/ST))
    d4  = d41+d42+d43

    c(exp(par[1])*d1,exp(par[2])*d2,exp(par[3])*d3,exp(par[4])*d4)

  }
  HL_function = function(par){

    Alpha1 = exp(par[1])
    Alpha2 = exp(par[2])
    Gamma1 = exp(par[3])
    Gamma2 = exp(par[4])

    h1 = Alpha1*Gamma1/(1+Alpha1*t.event)
    h2 = Alpha2*Gamma2/(1+Alpha2*t.event)
    S1 = (1+Alpha1*t.event)^(-Gamma1)
    S2 = (1+Alpha2*t.event)^(-Gamma2)
    p0 = exp(-Theta)
    p1 = exp(-Theta*S1)
    p2 = exp(-Theta*S2)
    ST = -(1/Theta)*log(1+(exp(-Theta*S1)-1)*(exp(-Theta*S2)-1)/(exp(-Theta)-1))

    der_h1_Alpha1 = Gamma1/(1+Alpha1*t.event)^2
    der_h2_Alpha2 = Gamma2/(1+Alpha2*t.event)^2
    der_S1_Alpha1 = -Gamma1*t.event*(1+Alpha1*t.event)^(-Gamma1-1)
    der_S2_Alpha2 = -Gamma2*t.event*(1+Alpha2*t.event)^(-Gamma2-1)
    der_h1_Gamma1 = Alpha1/(1+Alpha1*t.event)
    der_h2_Gamma2 = Alpha2/(1+Alpha2*t.event)
    der_S1_Gamma1 = -(1+Alpha1*t.event)^(-Gamma1)*log(1+Alpha1*t.event)
    der_S2_Gamma2 = -(1+Alpha2*t.event)^(-Gamma2)*log(1+Alpha2*t.event)
    der_ST_Alpha1 = der_S1_Alpha1*p1*(p2-1)/(p0-1+(p1-1)*(p2-1))
    der_ST_Alpha2 = der_S2_Alpha2*p2*(p1-1)/(p0-1+(p1-1)*(p2-1))
    der_ST_Gamma1 = der_S1_Gamma1*p1*(p2-1)/(p0-1+(p1-1)*(p2-1))
    der_ST_Gamma2 = der_S2_Gamma2*p2*(p1-1)/(p0-1+(p1-1)*(p2-1))

    der_h1_Alpha1_Alpha1 = -2*Gamma1*t.event/(1+Alpha1*t.event)^3
    der_h2_Alpha2_Alpha2 = -2*Gamma2*t.event/(1+Alpha2*t.event)^3
    der_S1_Alpha1_Alpha1 = Gamma1*(Gamma1+1)*t.event^2*(1+Alpha1*t.event)^(-Gamma1-2)
    der_S2_Alpha2_Alpha2 = Gamma2*(Gamma2+1)*t.event^2*(1+Alpha2*t.event)^(-Gamma2-2)
    der_h1_Gamma1_Gamma1 = 0
    der_h2_Gamma2_Gamma2 = 0
    der_S1_Gamma1_Gamma1 = (1+Alpha1*t.event)^(-Gamma1)*(log(1+Alpha1*t.event))^2
    der_S2_Gamma2_Gamma2 = (1+Alpha2*t.event)^(-Gamma2)*(log(1+Alpha2*t.event))^2
    der_h1_Alpha1_Gamma1 = (1+Alpha1*t.event)^(-2)
    der_h2_Alpha2_Gamma2 = (1+Alpha2*t.event)^(-2)
    der_S1_Alpha1_Gamma1 = t.event*(Gamma1*log(1+Alpha1*t.event)-1)/(1+Alpha1*t.event)^(Gamma1+1)
    der_S2_Alpha2_Gamma2 = t.event*(Gamma2*log(1+Alpha2*t.event)-1)/(1+Alpha2*t.event)^(Gamma2+1)

    der_ST_Alpha1_Alpha1 = ((p0-1+(p1-1)*(p2-1))*(p2-1)*(der_S1_Alpha1_Alpha1*p1-Theta*der_S1_Alpha1^2*p1)-der_S1_Alpha1*p1*(p2-1)*(-Theta*der_S1_Alpha1*p1*(p2-1)))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha2_Alpha2 = ((p0-1+(p1-1)*(p2-1))*(p1-1)*(der_S2_Alpha2_Alpha2*p2-Theta*der_S2_Alpha2^2*p2)-der_S2_Alpha2*p2*(p1-1)*(-Theta*der_S2_Alpha2*p2*(p1-1)))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha1_Alpha2 = ((p0-1+(p1-1)*(p2-1))*p1*p2*-Theta*der_S1_Alpha1*der_S2_Alpha2+Theta*p1*p2*(p1-1)*(p2-1)*der_S1_Alpha1*der_S2_Alpha2)/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Gamma1_Gamma1 = ((p0-1+(p1-1)*(p2-1))*(p2-1)*(der_S1_Gamma1_Gamma1*p1-Theta*der_S1_Gamma1^2*p1)-der_S1_Gamma1*p1*(p2-1)*(-Theta*der_S1_Gamma1*p1*(p2-1)))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Gamma2_Gamma2 = ((p0-1+(p1-1)*(p2-1))*(p1-1)*(der_S2_Gamma2_Gamma2*p2-Theta*der_S2_Gamma2^2*p2)-der_S2_Gamma2*p2*(p1-1)*(-Theta*der_S2_Gamma2*p2*(p1-1)))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Gamma1_Gamma2 = ((p0-1+(p1-1)*(p2-1))*p1*p2*-Theta*der_S1_Gamma1*der_S2_Gamma2+Theta*p1*p2*(p1-1)*(p2-1)*der_S1_Gamma1*der_S2_Gamma2)/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha1_Gamma1 = ((p0-1+(p1-1)*(p2-1))*(p2-1)*(der_S1_Alpha1_Gamma1*p1-Theta*der_S1_Alpha1*der_S1_Gamma1*p1)+Theta*der_S1_Alpha1*der_S1_Gamma1*p1^2*(p2-1)^2)/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha1_Gamma2 = ((p0-1+(p1-1)*(p2-1))*p1*der_S1_Alpha1*-Theta*der_S2_Gamma2*p2+Theta*der_S1_Alpha1*der_S2_Gamma2*p1*p2*(p1-1)*(p2-1))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha2_Gamma1 = ((p0-1+(p1-1)*(p2-1))*p1*p2*-Theta*der_S1_Gamma1*der_S2_Alpha2+Theta*p1*p2*(p1-1)*(p2-1)*der_S1_Gamma1*der_S2_Alpha2)/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha2_Gamma2 = ((p0-1+(p1-1)*(p2-1))*(p1-1)*(der_S2_Alpha2_Gamma2*p2-Theta*der_S2_Alpha2*der_S2_Gamma2*p2)+Theta*der_S2_Alpha2*der_S2_Gamma2*p2^2*(p1-1)^2)/(p0-1+(p1-1)*(p2-1))^2

    d11 = sum(event1*(der_h1_Alpha1/h1+der_S1_Alpha1/S1-Theta*der_S1_Alpha1+Theta*der_ST_Alpha1))
    d12 = sum(event2*(-Theta*der_S1_Alpha1*p1/(p1-1)+Theta*der_ST_Alpha1))
    d13 = sum((1-event1-event2)*(der_ST_Alpha1/ST))
    d1  = d11+d12+d13

    d21 = sum(event1*(-Theta*der_S2_Alpha2*p2/(p2-1)+Theta*der_ST_Alpha2))
    d22 = sum(event2*(der_h2_Alpha2/h2+der_S2_Alpha2/S2-Theta*der_S2_Alpha2+Theta*der_ST_Alpha2))
    d23 = sum((1-event1-event2)*(der_ST_Alpha2/ST))
    d2  = d21+d22+d23

    d31 = sum(event1*(der_h1_Gamma1/h1+der_S1_Gamma1/S1-Theta*der_S1_Gamma1+Theta*der_ST_Gamma1))
    d32 = sum(event2*(-Theta*der_S1_Gamma1*p1/(p1-1)+Theta*der_ST_Gamma1))
    d33 = sum((1-event1-event2)*(der_ST_Gamma1/ST))
    d3  = d31+d32+d33

    d41 = sum(event1*(-Theta*der_S2_Gamma2*p2/(p2-1)+Theta*der_ST_Gamma2))
    d42 = sum(event2*(der_h2_Gamma2/h2+der_S2_Gamma2/S2-Theta*der_S2_Gamma2+Theta*der_ST_Gamma2))
    d43 = sum((1-event1-event2)*(der_ST_Gamma2/ST))
    d4  = d41+d42+d43

    D111 = sum(event1*((der_h1_Alpha1_Alpha1*h1-der_h1_Alpha1^2)/h1^2+(der_S1_Alpha1_Alpha1*S1-der_S1_Alpha1^2)/S1^2-Theta*der_S1_Alpha1_Alpha1+Theta*der_ST_Alpha1_Alpha1))
    D112 = sum(event2*(((p1-1)*(-Theta*der_S1_Alpha1_Alpha1*p1+Theta^2*der_S1_Alpha1^2*p1)-Theta^2*der_S1_Alpha1^2*p1^2)/(p1-1)^2+Theta*der_ST_Alpha1_Alpha1))
    D113 = sum((1-event1-event2)*((ST*der_ST_Alpha1_Alpha1-der_ST_Alpha1^2)/ST^2))
    D11  = D111+D112+D113

    D221 = sum(event1*(((p2-1)*(-Theta*der_S2_Alpha2_Alpha2*p2+Theta^2*der_S2_Alpha2^2*p2)-Theta^2*der_S2_Alpha2^2*p2^2)/(p2-1)^2+Theta*der_ST_Alpha2_Alpha2))
    D222 = sum(event2*((der_h2_Alpha2_Alpha2*h2-der_h2_Alpha2^2)/h2^2+(der_S2_Alpha2_Alpha2*S2-der_S2_Alpha2^2)/S2^2-Theta*der_S2_Alpha2_Alpha2+Theta*der_ST_Alpha2_Alpha2))
    D223 = sum((1-event1-event2)*((ST*der_ST_Alpha2_Alpha2-der_ST_Alpha2^2)/ST^2))
    D22  = D221+D222+D223

    D121 = sum(event1*(Theta*der_ST_Alpha1_Alpha2))
    D122 = sum(event2*(Theta*der_ST_Alpha1_Alpha2))
    D123 = sum((1-event1-event2)*((ST*der_ST_Alpha1_Alpha2-der_ST_Alpha1*der_ST_Alpha2)/ST^2))
    D12  = D121+D122+D123

    D331 = sum(event1*((der_h1_Gamma1_Gamma1*h1-der_h1_Gamma1^2)/h1^2+(der_S1_Gamma1_Gamma1*S1-der_S1_Gamma1^2)/S1^2-Theta*der_S1_Gamma1_Gamma1+Theta*der_ST_Gamma1_Gamma1))
    D332 = sum(event2*(((p1-1)*(-Theta*der_S1_Gamma1_Gamma1*p1+Theta^2*der_S1_Gamma1^2*p1)-Theta^2*der_S1_Gamma1^2*p1^2)/(p1-1)^2+Theta*der_ST_Gamma1_Gamma1))
    D333 = sum((1-event1-event2)*((ST*der_ST_Gamma1_Gamma1-der_ST_Gamma1^2)/ST^2))
    D33  = D331+D332+D333

    D441 = sum(event1*(((p2-1)*(-Theta*der_S2_Gamma2_Gamma2*p2+Theta^2*der_S2_Gamma2^2*p2)-Theta^2*der_S2_Gamma2^2*p2^2)/(p2-1)^2+Theta*der_ST_Gamma2_Gamma2))
    D442 = sum(event2*((der_h2_Gamma2_Gamma2*h2-der_h2_Gamma2^2)/h2^2+(der_S2_Gamma2_Gamma2*S2-der_S2_Gamma2^2)/S2^2-Theta*der_S2_Gamma2_Gamma2+Theta*der_ST_Gamma2_Gamma2))
    D443 = sum((1-event1-event2)*((ST*der_ST_Gamma2_Gamma2-der_ST_Gamma2^2)/ST^2))
    D44  = D441+D442+D443

    D341 = sum(event1*(Theta*der_ST_Gamma1_Gamma2))
    D342 = sum(event2*(Theta*der_ST_Gamma1_Gamma2))
    D343 = sum((1-event1-event2)*((ST*der_ST_Gamma1_Gamma2-der_ST_Gamma1*der_ST_Gamma2)/ST^2))
    D34  = D341+D342+D343

    D131 = sum(event1*((der_h1_Alpha1_Gamma1*h1-der_h1_Alpha1*der_h1_Gamma1)/h1^2  +(der_S1_Alpha1_Gamma1*S1-der_S1_Alpha1*der_S1_Gamma1)/S1^2-Theta*der_S1_Alpha1_Gamma1+Theta*der_ST_Alpha1_Gamma1))
    D132 = sum(event2*(((p1-1)*(-Theta*der_S1_Alpha1_Gamma1*p1)-Theta^2*der_S1_Alpha1*der_S1_Gamma1*p1)/(p1-1)^2+Theta*der_ST_Alpha1_Gamma1))
    D133 = sum((1-event1-event2)*((ST*der_ST_Alpha1_Gamma1-der_ST_Alpha1*der_ST_Gamma1)/ST^2))
    D13  = D131+D132+D133

    D141 = sum(event1*(Theta*der_ST_Alpha1_Gamma2))
    D142 = sum(event2*(Theta*der_ST_Alpha1_Gamma2))
    D143 = sum((1-event1-event2)*((ST*der_ST_Alpha1_Gamma2-der_ST_Alpha1*der_ST_Gamma2)/ST^2))
    D14  = D141+D142+D143

    D231 = sum(event1*(Theta*der_ST_Alpha2_Gamma1))
    D232 = sum(event2*(Theta*der_ST_Alpha2_Gamma1))
    D233 = sum((1-event1-event2)*((ST*der_ST_Alpha2_Gamma1-der_ST_Alpha2*der_ST_Gamma1)/ST^2))
    D23  = D231+D232+D233

    D241 = sum(event1*(((p2-1)*(-Theta*der_S2_Alpha2_Gamma2*p2+Theta^2*der_S2_Alpha2*der_S2_Gamma2*p2)-Theta^2*der_S2_Alpha2*der_S2_Gamma2*p2^2)/(p2-1)^2+Theta*der_ST_Alpha2_Gamma2))
    D242 = sum(event2*((der_h2_Alpha2_Gamma2*h2-der_h2_Alpha2*der_h2_Gamma2)/h2^2+(der_S2_Alpha2_Gamma2*S2-der_S2_Alpha2*der_S2_Gamma2)/S2^2-Theta*der_S2_Alpha2_Gamma2+Theta*der_ST_Alpha2_Gamma2))
    D243 = sum((1-event1-event2)*((ST*der_ST_Alpha2_Gamma2-der_ST_Alpha2*der_ST_Gamma2)/ST^2))
    D24  = D241+D242+D243

    DD11 = exp(2*par[1])*D11+exp(par[1])*d1
    DD12 = exp(par[1])*exp(par[2])*D12
    DD13 = exp(par[1])*exp(par[3])*D13
    DD14 = exp(par[1])*exp(par[4])*D14
    DD22 = exp(2*par[2])*D22+exp(par[2])*d2
    DD23 = exp(par[2])*exp(par[3])*D23
    DD24 = exp(par[2])*exp(par[4])*D24
    DD33 = exp(2*par[3])*D33+exp(par[3])*d3
    DD34 = exp(par[3])*exp(par[4])*D34
    DD44 = exp(2*par[4])*D44+exp(par[4])*d4

    matrix(c(DD11,DD12,DD13,DD14,DD12,DD22,DD23,DD24,DD13,DD23,DD33,DD34,DD14,DD24,DD34,DD44),4,4)

  }
  H_function = function(par){

    Alpha1 = par[1]
    Alpha2 = par[2]
    Gamma1 = par[3]
    Gamma2 = par[4]

    h1 = Alpha1*Gamma1/(1+Alpha1*t.event)
    h2 = Alpha2*Gamma2/(1+Alpha2*t.event)
    S1 = (1+Alpha1*t.event)^(-Gamma1)
    S2 = (1+Alpha2*t.event)^(-Gamma2)
    p0 = exp(-Theta)
    p1 = exp(-Theta*S1)
    p2 = exp(-Theta*S2)
    ST = -(1/Theta)*log(1+(exp(-Theta*S1)-1)*(exp(-Theta*S2)-1)/(exp(-Theta)-1))

    der_h1_Alpha1 = Gamma1/(1+Alpha1*t.event)^2
    der_h2_Alpha2 = Gamma2/(1+Alpha2*t.event)^2
    der_S1_Alpha1 = -Gamma1*t.event*(1+Alpha1*t.event)^(-Gamma1-1)
    der_S2_Alpha2 = -Gamma2*t.event*(1+Alpha2*t.event)^(-Gamma2-1)
    der_h1_Gamma1 = Alpha1/(1+Alpha1*t.event)
    der_h2_Gamma2 = Alpha2/(1+Alpha2*t.event)
    der_S1_Gamma1 = -(1+Alpha1*t.event)^(-Gamma1)*log(1+Alpha1*t.event)
    der_S2_Gamma2 = -(1+Alpha2*t.event)^(-Gamma2)*log(1+Alpha2*t.event)
    der_ST_Alpha1 = der_S1_Alpha1*p1*(p2-1)/(p0-1+(p1-1)*(p2-1))
    der_ST_Alpha2 = der_S2_Alpha2*p2*(p1-1)/(p0-1+(p1-1)*(p2-1))
    der_ST_Gamma1 = der_S1_Gamma1*p1*(p2-1)/(p0-1+(p1-1)*(p2-1))
    der_ST_Gamma2 = der_S2_Gamma2*p2*(p1-1)/(p0-1+(p1-1)*(p2-1))

    der_h1_Alpha1_Alpha1 = -2*Gamma1*t.event/(1+Alpha1*t.event)^3
    der_h2_Alpha2_Alpha2 = -2*Gamma2*t.event/(1+Alpha2*t.event)^3
    der_S1_Alpha1_Alpha1 = Gamma1*(Gamma1+1)*t.event^2*(1+Alpha1*t.event)^(-Gamma1-2)
    der_S2_Alpha2_Alpha2 = Gamma2*(Gamma2+1)*t.event^2*(1+Alpha2*t.event)^(-Gamma2-2)
    der_h1_Gamma1_Gamma1 = 0
    der_h2_Gamma2_Gamma2 = 0
    der_S1_Gamma1_Gamma1 = (1+Alpha1*t.event)^(-Gamma1)*(log(1+Alpha1*t.event))^2
    der_S2_Gamma2_Gamma2 = (1+Alpha2*t.event)^(-Gamma2)*(log(1+Alpha2*t.event))^2
    der_h1_Alpha1_Gamma1 = (1+Alpha1*t.event)^(-2)
    der_h2_Alpha2_Gamma2 = (1+Alpha2*t.event)^(-2)
    der_S1_Alpha1_Gamma1 = t.event*(Gamma1*log(1+Alpha1*t.event)-1)/(1+Alpha1*t.event)^(Gamma1+1)
    der_S2_Alpha2_Gamma2 = t.event*(Gamma2*log(1+Alpha2*t.event)-1)/(1+Alpha2*t.event)^(Gamma2+1)

    der_ST_Alpha1_Alpha1 = ((p0-1+(p1-1)*(p2-1))*(p2-1)*(der_S1_Alpha1_Alpha1*p1-Theta*der_S1_Alpha1^2*p1)-der_S1_Alpha1*p1*(p2-1)*(-Theta*der_S1_Alpha1*p1*(p2-1)))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha2_Alpha2 = ((p0-1+(p1-1)*(p2-1))*(p1-1)*(der_S2_Alpha2_Alpha2*p2-Theta*der_S2_Alpha2^2*p2)-der_S2_Alpha2*p2*(p1-1)*(-Theta*der_S2_Alpha2*p2*(p1-1)))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha1_Alpha2 = ((p0-1+(p1-1)*(p2-1))*p1*p2*-Theta*der_S1_Alpha1*der_S2_Alpha2+Theta*p1*p2*(p1-1)*(p2-1)*der_S1_Alpha1*der_S2_Alpha2)/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Gamma1_Gamma1 = ((p0-1+(p1-1)*(p2-1))*(p2-1)*(der_S1_Gamma1_Gamma1*p1-Theta*der_S1_Gamma1^2*p1)-der_S1_Gamma1*p1*(p2-1)*(-Theta*der_S1_Gamma1*p1*(p2-1)))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Gamma2_Gamma2 = ((p0-1+(p1-1)*(p2-1))*(p1-1)*(der_S2_Gamma2_Gamma2*p2-Theta*der_S2_Gamma2^2*p2)-der_S2_Gamma2*p2*(p1-1)*(-Theta*der_S2_Gamma2*p2*(p1-1)))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Gamma1_Gamma2 = ((p0-1+(p1-1)*(p2-1))*p1*p2*-Theta*der_S1_Gamma1*der_S2_Gamma2+Theta*p1*p2*(p1-1)*(p2-1)*der_S1_Gamma1*der_S2_Gamma2)/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha1_Gamma1 = ((p0-1+(p1-1)*(p2-1))*(p2-1)*(der_S1_Alpha1_Gamma1*p1-Theta*der_S1_Alpha1*der_S1_Gamma1*p1)+Theta*der_S1_Alpha1*der_S1_Gamma1*p1^2*(p2-1)^2)/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha1_Gamma2 = ((p0-1+(p1-1)*(p2-1))*p1*der_S1_Alpha1*-Theta*der_S2_Gamma2*p2+Theta*der_S1_Alpha1*der_S2_Gamma2*p1*p2*(p1-1)*(p2-1))/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha2_Gamma1 = ((p0-1+(p1-1)*(p2-1))*p1*p2*-Theta*der_S1_Gamma1*der_S2_Alpha2+Theta*p1*p2*(p1-1)*(p2-1)*der_S1_Gamma1*der_S2_Alpha2)/(p0-1+(p1-1)*(p2-1))^2
    der_ST_Alpha2_Gamma2 = ((p0-1+(p1-1)*(p2-1))*(p1-1)*(der_S2_Alpha2_Gamma2*p2-Theta*der_S2_Alpha2*der_S2_Gamma2*p2)+Theta*der_S2_Alpha2*der_S2_Gamma2*p2^2*(p1-1)^2)/(p0-1+(p1-1)*(p2-1))^2

    d11 = sum(event1*(der_h1_Alpha1/h1+der_S1_Alpha1/S1-Theta*der_S1_Alpha1+Theta*der_ST_Alpha1))
    d12 = sum(event2*(-Theta*der_S1_Alpha1*p1/(p1-1)+Theta*der_ST_Alpha1))
    d13 = sum((1-event1-event2)*(der_ST_Alpha1/ST))
    d1  = d11+d12+d13

    d21 = sum(event1*(-Theta*der_S2_Alpha2*p2/(p2-1)+Theta*der_ST_Alpha2))
    d22 = sum(event2*(der_h2_Alpha2/h2+der_S2_Alpha2/S2-Theta*der_S2_Alpha2+Theta*der_ST_Alpha2))
    d23 = sum((1-event1-event2)*(der_ST_Alpha2/ST))
    d2  = d21+d22+d23

    d31 = sum(event1*(der_h1_Gamma1/h1+der_S1_Gamma1/S1-Theta*der_S1_Gamma1+Theta*der_ST_Gamma1))
    d32 = sum(event2*(-Theta*der_S1_Gamma1*p1/(p1-1)+Theta*der_ST_Gamma1))
    d33 = sum((1-event1-event2)*(der_ST_Gamma1/ST))
    d3  = d31+d32+d33

    d41 = sum(event1*(-Theta*der_S2_Gamma2*p2/(p2-1)+Theta*der_ST_Gamma2))
    d42 = sum(event2*(der_h2_Gamma2/h2+der_S2_Gamma2/S2-Theta*der_S2_Gamma2+Theta*der_ST_Gamma2))
    d43 = sum((1-event1-event2)*(der_ST_Gamma2/ST))
    d4  = d41+d42+d43

    D111 = sum(event1*((der_h1_Alpha1_Alpha1*h1-der_h1_Alpha1^2)/h1^2+(der_S1_Alpha1_Alpha1*S1-der_S1_Alpha1^2)/S1^2-Theta*der_S1_Alpha1_Alpha1+Theta*der_ST_Alpha1_Alpha1))
    D112 = sum(event2*(((p1-1)*(-Theta*der_S1_Alpha1_Alpha1*p1+Theta^2*der_S1_Alpha1^2*p1)-Theta^2*der_S1_Alpha1^2*p1^2)/(p1-1)^2+Theta*der_ST_Alpha1_Alpha1))
    D113 = sum((1-event1-event2)*((ST*der_ST_Alpha1_Alpha1-der_ST_Alpha1^2)/ST^2))
    D11  = D111+D112+D113

    D221 = sum(event1*(((p2-1)*(-Theta*der_S2_Alpha2_Alpha2*p2+Theta^2*der_S2_Alpha2^2*p2)-Theta^2*der_S2_Alpha2^2*p2^2)/(p2-1)^2+Theta*der_ST_Alpha2_Alpha2))
    D222 = sum(event2*((der_h2_Alpha2_Alpha2*h2-der_h2_Alpha2^2)/h2^2+(der_S2_Alpha2_Alpha2*S2-der_S2_Alpha2^2)/S2^2-Theta*der_S2_Alpha2_Alpha2+Theta*der_ST_Alpha2_Alpha2))
    D223 = sum((1-event1-event2)*((ST*der_ST_Alpha2_Alpha2-der_ST_Alpha2^2)/ST^2))
    D22  = D221+D222+D223

    D121 = sum(event1*(Theta*der_ST_Alpha1_Alpha2))
    D122 = sum(event2*(Theta*der_ST_Alpha1_Alpha2))
    D123 = sum((1-event1-event2)*((ST*der_ST_Alpha1_Alpha2-der_ST_Alpha1*der_ST_Alpha2)/ST^2))
    D12  = D121+D122+D123

    D331 = sum(event1*((der_h1_Gamma1_Gamma1*h1-der_h1_Gamma1^2)/h1^2+(der_S1_Gamma1_Gamma1*S1-der_S1_Gamma1^2)/S1^2-Theta*der_S1_Gamma1_Gamma1+Theta*der_ST_Gamma1_Gamma1))
    D332 = sum(event2*(((p1-1)*(-Theta*der_S1_Gamma1_Gamma1*p1+Theta^2*der_S1_Gamma1^2*p1)-Theta^2*der_S1_Gamma1^2*p1^2)/(p1-1)^2+Theta*der_ST_Gamma1_Gamma1))
    D333 = sum((1-event1-event2)*((ST*der_ST_Gamma1_Gamma1-der_ST_Gamma1^2)/ST^2))
    D33  = D331+D332+D333

    D441 = sum(event1*(((p2-1)*(-Theta*der_S2_Gamma2_Gamma2*p2+Theta^2*der_S2_Gamma2^2*p2)-Theta^2*der_S2_Gamma2^2*p2^2)/(p2-1)^2+Theta*der_ST_Gamma2_Gamma2))
    D442 = sum(event2*((der_h2_Gamma2_Gamma2*h2-der_h2_Gamma2^2)/h2^2+(der_S2_Gamma2_Gamma2*S2-der_S2_Gamma2^2)/S2^2-Theta*der_S2_Gamma2_Gamma2+Theta*der_ST_Gamma2_Gamma2))
    D443 = sum((1-event1-event2)*((ST*der_ST_Gamma2_Gamma2-der_ST_Gamma2^2)/ST^2))
    D44  = D441+D442+D443

    D341 = sum(event1*(Theta*der_ST_Gamma1_Gamma2))
    D342 = sum(event2*(Theta*der_ST_Gamma1_Gamma2))
    D343 = sum((1-event1-event2)*((ST*der_ST_Gamma1_Gamma2-der_ST_Gamma1*der_ST_Gamma2)/ST^2))
    D34  = D341+D342+D343

    D131 = sum(event1*((der_h1_Alpha1_Gamma1*h1-der_h1_Alpha1*der_h1_Gamma1)/h1^2  +(der_S1_Alpha1_Gamma1*S1-der_S1_Alpha1*der_S1_Gamma1)/S1^2-Theta*der_S1_Alpha1_Gamma1+Theta*der_ST_Alpha1_Gamma1))
    D132 = sum(event2*(((p1-1)*(-Theta*der_S1_Alpha1_Gamma1*p1)-Theta^2*der_S1_Alpha1*der_S1_Gamma1*p1)/(p1-1)^2+Theta*der_ST_Alpha1_Gamma1))
    D133 = sum((1-event1-event2)*((ST*der_ST_Alpha1_Gamma1-der_ST_Alpha1*der_ST_Gamma1)/ST^2))
    D13  = D131+D132+D133

    D141 = sum(event1*(Theta*der_ST_Alpha1_Gamma2))
    D142 = sum(event2*(Theta*der_ST_Alpha1_Gamma2))
    D143 = sum((1-event1-event2)*((ST*der_ST_Alpha1_Gamma2-der_ST_Alpha1*der_ST_Gamma2)/ST^2))
    D14  = D141+D142+D143

    D231 = sum(event1*(Theta*der_ST_Alpha2_Gamma1))
    D232 = sum(event2*(Theta*der_ST_Alpha2_Gamma1))
    D233 = sum((1-event1-event2)*((ST*der_ST_Alpha2_Gamma1-der_ST_Alpha2*der_ST_Gamma1)/ST^2))
    D23  = D231+D232+D233

    D241 = sum(event1*(((p2-1)*(-Theta*der_S2_Alpha2_Gamma2*p2+Theta^2*der_S2_Alpha2*der_S2_Gamma2*p2)-Theta^2*der_S2_Alpha2*der_S2_Gamma2*p2^2)/(p2-1)^2+Theta*der_ST_Alpha2_Gamma2))
    D242 = sum(event2*((der_h2_Alpha2_Gamma2*h2-der_h2_Alpha2*der_h2_Gamma2)/h2^2+(der_S2_Alpha2_Gamma2*S2-der_S2_Alpha2*der_S2_Gamma2)/S2^2-Theta*der_S2_Alpha2_Gamma2+Theta*der_ST_Alpha2_Gamma2))
    D243 = sum((1-event1-event2)*((ST*der_ST_Alpha2_Gamma2-der_ST_Alpha2*der_ST_Gamma2)/ST^2))
    D24  = D241+D242+D243

    matrix(c(D11,D12,D13,D14,D12,D22,D23,D24,D13,D23,D33,D34,D14,D24,D34,D44),4,4)

  }

  par_old = c(log(Alpha1.0),log(Alpha2.0),log(Gamma1.0),log(Gamma2.0))
  count  = 0
  random = 0

  repeat{

    temp = try(solve(HL_function(par_old),silent = TRUE))
    if (is(temp,"try-error")){

      random = random+1
      count = 0
      par_old = c(log(Alpha1.0*exp(runif(1,-r.1,r.1))),
                  log(Alpha2.0*exp(runif(1,-r.2,r.2))),
                  log(Gamma1.0*exp(runif(1,-r.3,r.3))),
                  log(Gamma2.0*exp(runif(1,-r.4,r.4))))
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
                  log(Gamma1.0*exp(runif(1,-r.3,r.3))),
                  log(Gamma2.0*exp(runif(1,-r.4,r.4))))
      next

    }

    if (max(abs(exp(par_old)-exp(par_new))) < epsilon) {break}
    par_old = par_new

  }

  Alpha1_hat = exp(par_new[1])
  Alpha2_hat = exp(par_new[2])
  Gamma1_hat = exp(par_new[3])
  Gamma2_hat = exp(par_new[4])

  Info = solve(-H_function(exp(par_new)))
  Alpha1_se = sqrt(Info[1,1])
  Alpha2_se = sqrt(Info[2,2])
  Gamma1_se = sqrt(Info[3,3])
  Gamma2_se = sqrt(Info[4,4])

  InfoL = solve(-HL_function(par_new))
  CI_Alpha1 = c(Alpha1_hat*exp(-qnorm(0.975)*sqrt(InfoL[1,1])),
                Alpha1_hat*exp(+qnorm(0.975)*sqrt(InfoL[1,1])))
  CI_Alpha2 = c(Alpha2_hat*exp(-qnorm(0.975)*sqrt(InfoL[2,2])),
                Alpha2_hat*exp(+qnorm(0.975)*sqrt(InfoL[2,2])))
  CI_Gamma1 = c(Gamma1_hat*exp(-qnorm(0.975)*sqrt(InfoL[3,3])),
                Gamma1_hat*exp(+qnorm(0.975)*sqrt(InfoL[3,3])))
  CI_Gamma2 = c(Gamma2_hat*exp(-qnorm(0.975)*sqrt(InfoL[4,4])),
                Gamma2_hat*exp(+qnorm(0.975)*sqrt(InfoL[4,4])))

  MedX_hat = (2^(1/Gamma1_hat)-1)/Alpha1_hat
  MedY_hat = (2^(1/Gamma2_hat)-1)/Alpha2_hat
  transX = c((1-2^(1/Gamma1_hat))/Alpha1_hat^2,0,-2^(1/Gamma1_hat)*log(2)/(Alpha1_hat*Gamma1_hat^2),0)
  transY = c(0,(1-2^(1/Gamma2_hat))/Alpha2_hat^2,0,-2^(1/Gamma2_hat)*log(2)/(Alpha2_hat*Gamma2_hat^2))
  MedX_se = sqrt(t(transX)%*%Info%*%transX)
  MedY_se = sqrt(t(transY)%*%Info%*%transY)

  temp_transX = c(-1,0,-2^(1/Gamma1_hat)*log(2)/((2^(1/Gamma1_hat)-1)*Gamma1_hat),0)
  temp_transY = c(0,-1,0,-2^(1/Gamma2_hat)*log(2)/((2^(1/Gamma2_hat)-1)*Gamma2_hat))
  temp_MedX_se = sqrt(t(temp_transX)%*%InfoL%*%temp_transX)
  temp_MedY_se = sqrt(t(temp_transY)%*%InfoL%*%temp_transY)

  CI_MedX = c(MedX_hat*exp(-qnorm(0.975)*temp_MedX_se),
              MedX_hat*exp(+qnorm(0.975)*temp_MedX_se))
  CI_MedY = c(MedY_hat*exp(-qnorm(0.975)*temp_MedY_se),
              MedY_hat*exp(+qnorm(0.975)*temp_MedY_se))

  Alpha1.res = c(Estimate = Alpha1_hat,SE = Alpha1_se,CI.lower = CI_Alpha1[1],CI.upper = CI_Alpha1[2])
  Alpha2.res = c(Estimate = Alpha2_hat,SE = Alpha2_se,CI.lower = CI_Alpha2[1],CI.upper = CI_Alpha2[2])
  Gamma1.res = c(Estimate = Gamma1_hat,SE = Gamma1_se,CI.lower = CI_Gamma1[1],CI.upper = CI_Gamma1[2])
  Gamma2.res = c(Estimate = Gamma2_hat,SE = Gamma2_se,CI.lower = CI_Gamma2[1],CI.upper = CI_Gamma2[2])

  MedX.res = c(Estimate = MedX_hat,SE = MedX_se,CI.lower = CI_MedX[1],CI.upper = CI_MedX[2])
  MedY.res = c(Estimate = MedY_hat,SE = MedY_se,CI.lower = CI_MedY[1],CI.upper = CI_MedY[2])

  if (Gamma1_hat < 1 & Gamma2_hat < 1) {

    return(list(n = n,Iteration = count,Randomization = random,
                Alpha1 = Alpha1.res,Alpha2 = Alpha2.res,Gamma1 = Gamma1.res,Gamma2 = Gamma2.res,
                MedX = MedX.res,MedY = MedY.res,MeanX = "Unavaliable",MeanY = "Unavaliable",
                logL = log_L(par_new),AIC = 2*length(par_new)-2*log_L(par_new),
                BIC = length(par_new)*log(length(t.event))-2*log_L(par_new)))

  } else if (Gamma1_hat >= 1 & Gamma2_hat >= 1) {

    MeanX_hat = 1/(Alpha1_hat*(Gamma1_hat-1))
    MeanY_hat = 1/(Alpha2_hat*(Gamma2_hat-1))
    trans2X = c(-1/(Alpha1_hat^2*(Gamma1_hat-1)),0,-1/(Alpha1_hat*(Gamma1_hat-1)^2),0)
    trans2Y = c(0,-1/(Alpha2_hat^2*(Gamma2_hat-1)),0,-1/(Alpha2_hat*(Gamma2_hat-1)^2))
    MeanX_se = sqrt(t(trans2X)%*%Info%*%trans2X)
    MeanY_se = sqrt(t(trans2Y)%*%Info%*%trans2Y)

    temp_trans2X = c(-1,0,-Gamma1_hat/(Gamma1_hat-1),0)
    temp_trans2Y = c(0,-1,0,-Gamma2_hat/(Gamma2_hat-1))
    temp_MeanX_se = sqrt(t(temp_trans2X)%*%InfoL%*%temp_trans2X)
    temp_MeanY_se = sqrt(t(temp_trans2Y)%*%InfoL%*%temp_trans2Y)

    CI_MeanX = c(MeanX_hat*exp(-qnorm(0.975)*temp_MeanX_se),
                 MeanX_hat*exp(+qnorm(0.975)*temp_MeanX_se))
    CI_MeanY = c(MeanY_hat*exp(-qnorm(0.975)*temp_MeanY_se),
                 MeanY_hat*exp(+qnorm(0.975)*temp_MeanY_se))

    MeanX.res = c(Estimate = MeanX_hat,SE = MeanX_se,CI.lower = CI_MeanX[1],CI.upper = CI_MeanX[2])
    MeanY.res = c(Estimate = MeanY_hat,SE = MeanY_se,CI.lower = CI_MeanY[1],CI.upper = CI_MeanY[2])

    return(list(n = n,Iteration = count,Randomization = random,
                Alpha1 = Alpha1.res,Alpha2 = Alpha2.res,Gamma1 = Gamma1.res,Gamma2 = Gamma2.res,
                MedX = MedX.res,MedY = MedY.res,MeanX = MeanX.res,MeanY = MeanY.res,
                logL = log_L(par_new),AIC = 2*length(par_new)-2*log_L(par_new),
                BIC = length(par_new)*log(length(t.event))-2*log_L(par_new)))

  } else if (Gamma1_hat >= 1 & Gamma2_hat < 1) {

    MeanX_hat = 1/(Alpha1_hat*(Gamma1_hat-1))
    trans2X = c(-1/(Alpha1_hat^2*(Gamma1_hat-1)),0,-1/(Alpha1_hat*(Gamma1_hat-1)^2),0)
    MeanX_se = sqrt(t(trans2X)%*%Info%*%trans2X)

    temp_trans2X = c(-1,0,-Gamma1_hat/(Gamma1_hat-1),0)
    temp_MeanX_se = sqrt(t(temp_trans2X)%*%InfoL%*%temp_trans2X)

    CI_MeanX = c(MeanX_hat*exp(-qnorm(0.975)*temp_MeanX_se),
                 MeanX_hat*exp(+qnorm(0.975)*temp_MeanX_se))

    MeanX.res = c(Estimate = MeanX_hat,SE = MeanX_se,CI.lower = CI_MeanX[1],CI.upper = CI_MeanX[2])

    return(list(n = n,Iteration = count,Randomization = random,
                Alpha1 = Alpha1.res,Alpha2 = Alpha2.res,Gamma1 = Gamma1.res,Gamma2 = Gamma2.res,
                MedX = MedX.res,MedY = MedY.res,MeanX = MeanX.res,MeanY = "Unavaliable",
                logL = log_L(par_new),AIC = 2*length(par_new)-2*log_L(par_new),
                BIC = length(par_new)*log(length(t.event))-2*log_L(par_new)))

  } else {

    MeanY_hat = 1/(Alpha2_hat*(Gamma2_hat-1))
    trans2Y = c(0,-1/(Alpha2_hat^2*(Gamma2_hat-1)),0,-1/(Alpha2_hat*(Gamma2_hat-1)^2))
    MeanY_se = sqrt(t(trans2Y)%*%Info%*%trans2Y)

    temp_trans2Y = c(0,-1,0,-Gamma2_hat/(Gamma2_hat-1))
    temp_MeanY_se = sqrt(t(temp_trans2Y)%*%InfoL%*%temp_trans2Y)

    CI_MeanY = c(MeanY_hat*exp(-qnorm(0.975)*temp_MeanY_se),
                 MeanY_hat*exp(+qnorm(0.975)*temp_MeanY_se))

    MeanY.res = c(Estimate = MeanY_hat,SE = MeanY_se,CI.lower = CI_MeanY[1],CI.upper = CI_MeanY[2])

    return(list(n = n,Iteration = count,Randomization = random,
                Alpha1 = Alpha1.res,Alpha2 = Alpha2.res,Gamma1 = Gamma1.res,Gamma2 = Gamma2.res,
                MedX = MedX.res,MedY = MedY.res,MeanX = "Unavaliable",MeanY = MeanY.res,
                logL = log_L(par_new),AIC = 2*length(par_new)-2*log_L(par_new),
                BIC = length(par_new)*log(length(t.event))-2*log_L(par_new)))

  }

}
