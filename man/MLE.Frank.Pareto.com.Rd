% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLE.Frank.Pareto.com.R
\name{MLE.Frank.Pareto.com}
\alias{MLE.Frank.Pareto.com}
\title{Maximum likelihood estimation for bivariate dependent competing risks data under the Frank copula with the common Pareto margins}
\usage{
MLE.Frank.Pareto.com(
  t.event,
  event1,
  event2,
  Theta.0 = 1,
  Alpha.0 = 1,
  Gamma.0 = 1,
  epsilon = 1e-05,
  r.1 = 13,
  r.2 = 3,
  r.3 = 3,
  bootstrap = FALSE,
  B = 200
)
}
\arguments{
\item{t.event}{Vector of the observed failure times.}

\item{event1}{Vector of the indicators for the failure cause 1.}

\item{event2}{Vector of the indicators for the failure cause 2.}

\item{Theta.0}{Initial guess for the copula parameter \eqn{\theta}.}

\item{Alpha.0}{Initial guess for the common scale parameter \eqn{\alpha} with default value 1.}

\item{Gamma.0}{Initial guess for the common shape parameter \eqn{\gamma} with default value 1.}

\item{epsilon}{Positive tunning parameter in the NR algorithm with default value \eqn{10^{-5}}.}

\item{r.1}{Positive tunning parameter in the NR algorithm with default value 1.}

\item{r.2}{Positive tunning parameter in the NR algorithm with default value 1.}

\item{r.3}{Positive tunning parameter in the NR algorithm with default value 1.}

\item{bootstrap}{Perform parametric bootstrap if \code{TRUE}.}

\item{B}{Number of bootstrap replications.}
}
\value{
\item{n}{Sample size.}
\item{count}{Iteration number.}
\item{random}{Randomization number.}
\item{Theta}{Copula parameter.}
\item{Theta.B}{Copula parameter (SE and CI are calculated by parametric bootstrap method).}
\item{Alpha}{Common positive scale parameter for the Pareto margin.}
\item{Alpha.B}{Common positive scale parameter for the Pareto margin (SE and CI are calculated by parametric bootstrap method).}
\item{Gamma}{Common positive shape parameter for the Pareto margin.}
\item{Gamma.B}{Common positive shape parameter for the Pareto margin (SE and CI are calculated by parametric bootstrap method).}
\item{logL}{Log-likelihood value under the fitted model.}
\item{AIC}{AIC value under the fitted model.}
\item{BIC}{BIC value under the fitted model.}
}
\description{
Maximum likelihood estimation for bivariate dependent competing risks data under the Frank copula with the common Pareto margins.
}
\details{
The parametric bootstrap method requires the assumption of the uniform censoring distribution. One must notice that such assumption is not always true in real data analysis.
}
\examples{
t.event = c(72,40,20,65,24,46,62,61,60,60,59,59,49,20, 3,58,29,26,52,20,
            51,51,31,42,38,69,39,33, 8,13,33, 9,21,66, 5,27, 2,20,19,60,
            32,53,53,43,21,74,72,14,33, 8,10,51, 7,33, 3,43,37, 5, 6, 2,
            5,64, 1,21,16,21,12,75,74,54,73,36,59, 6,58,16,19,39,26,60,
            43, 7, 9,67,62,17,25, 0, 5,34,59,31,58,30,57, 5,55,55,52, 0,
            51,17,70,74,74,20, 2, 8,27,23, 1,52,51, 6, 0,26,65,26, 6, 6,
            68,33,67,23, 6,11, 6,57,57,29, 9,53,51, 8, 0,21,27,22,12,68,
            21,68, 0, 2,14,18, 5,60,40,51,50,46,65, 9,21,27,54,52,75,30,
            70,14, 0,42,12,40, 2,12,53,11,18,13,45, 8,28,67,67,24,64,26,
            57,32,42,20,71,54,64,51, 1, 2, 0,54,69,68,67,66,64,63,35,62,
            7,35,24,57, 1, 4,74, 0,51,36,16,32,68,17,66,65,19,41,28, 0,
            46,63,60,59,46,63, 8,74,18,33,12, 1,66,28,30,57,50,39,40,24,
            6,30,58,68,24,33,65, 2,64,19,15,10,12,53,51, 1,40,40,66, 2,
            21,35,29,54,37,10,29,71,12,13,27,66,28,31,12, 9,21,19,51,71,
            76,46,47,75,75,49,75,75,31,69,74,25,72,28,36, 8,71,60,14,22,
            67,62,68,68,27,68,68,67,67, 3,49,12,30,67, 5,65,24,66,36,66,
            40,13,40, 0,14,45,64,13,24,15,26, 5,63,35,61,61,50,57,21,26,
            11,59,42,27,50,57,57, 0, 1,54,53,23, 8,51,27,52,52,52,45,48,
            18, 2, 2,35,75,75, 9,39, 0,26,17,43,53,47,11,65,16,21,64, 7,
            38,55, 5,28,38,20,24,27,31, 9, 9,11,56,36,56,15,51,33,70,32,
            5,23,63,30,53,12,58,54,36,20,74,34,70,25,65, 4,10,58,37,56,
            6, 0,70,70,28,40,67,36,23,23,62,62,62, 2,34, 4,12,56, 1, 7,
            4,70,65, 7,30,40,13,22, 0,18,64,13,26, 1,16,33,22,30,53,53,
            7,61,40, 9,59, 7,12,46,50, 0,52,19,52,51,51,14,27,51, 5, 0,
            41,53,19)

event1 = c(0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
           0,0,1,0,0,0,1,0,1,1,0,1,1,1,1,0,0,1,1,0,
           1,0,0,1,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,1,
           1,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,1,0,0,
           0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,
           0,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,0,0,0,
           1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,1,0,0,1,1,0,1,0,0,1,1,0,0,
           1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
           0,0,1,0,1,0,0,0,0,1,1,1,1,0,0,0,1,1,0,0,
           1,1,1,1,0,0,1,0,1,1,1,1,1,1,1,0,1,1,0,1,
           0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
           0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,
           0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
           1,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,0,1,
           1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,0,0,0,0,
           0,0,0,1,0,0,0,0,1,0,0,1,0,1,0,1,1,0,1,0,
           1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,
           1,0,0,1,0,0,0,1,0,1,0,0,1,0,0,0,1,1,0,1,
           1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,
           0,0,1)

event2 = c(0,1,1,0,0,1,0,0,0,0,0,0,0,1,1,0,1,1,0,1,
           0,0,0,1,1,0,0,1,0,0,1,0,0,0,0,1,1,0,0,0,
           0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,
           0,0,1,0,1,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1,
           1,1,1,0,1,1,1,1,1,1,0,1,0,1,0,1,0,0,0,1,
           0,1,1,0,0,1,0,0,1,1,1,0,0,0,0,1,1,0,1,1,
           0,1,0,0,1,1,0,0,0,1,1,0,0,1,1,1,0,1,0,0,
           1,0,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,0,0,1,
           0,1,1,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,0,1,
           0,0,1,1,0,1,0,1,1,1,0,1,0,0,0,0,0,0,1,0,
           1,1,1,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,1,
           0,0,0,0,1,0,1,0,1,1,1,1,0,1,1,1,0,1,1,1,
           1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,
           0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
           0,0,1,0,0,1,0,0,1,0,0,1,0,1,1,0,0,1,1,1,
           1,1,0,0,1,0,0,0,0,1,1,1,1,0,1,1,1,0,1,0,
           1,1,1,1,1,1,0,1,1,1,1,0,0,1,0,0,1,1,1,0,
           1,0,0,1,1,0,0,1,1,0,0,1,1,1,1,0,0,0,1,1,
           0,1,1,1,0,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,
           0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0,1,
           1,1,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
           0,1,0,0,1,1,0,1,1,1,0,0,0,1,0,1,0,0,1,1,
           0,0,0,0,1,1,1,0,1,0,1,1,0,1,1,1,0,0,1,0,
           0,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,1,1,
           1,0,0)

library(Bivariate.Pareto)
set.seed(10)
MLE.Frank.Pareto.com(t.event,event1,event2,bootstrap = FALSE)
}
\references{
Shih J-H, Lee W, Sun L-H, Emura T (2019), Fitting competing risks data to bivariate Pareto models, Communications in Statistics - Theory and Methods, 48:1193-1220.
}
