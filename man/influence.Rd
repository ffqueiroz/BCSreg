% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5_diagnostic.R
\name{influence}
\alias{influence}
\title{Influence Diagnostics for BCSreg Objects}
\usage{
influence(object, plot = TRUE, ask = grDevices::dev.interactive(), ...)
}
\arguments{
\item{object}{an object of class \code{"BCSreg"}.}

\item{plot}{logical. If \code{plot = TRUE} (default), the plots are shown.}

\item{ask}{logical; if \code{TRUE}, the user is asked before each plot, if \code{plot = TRUE}.}

\item{...}{currently not used.}
}
\value{
\code{influence} returns a list with two objects:
    \item{case.weights}{The values of \eqn{d_{max}} eigenvector based on case
    weights perturbation scheme (see Medeiros and Queiroz (2025)).}
    \item{totalLI}{The total local influence (see Lesaffre and Verbeke (1998)).}
}
\description{
The \code{influence} function provides two influence measures for a Box-Cox symmetric or a
    zero-adusted Box-Cox symmetric regression fit.
}
\examples{
## Data set: fishery (for description, run ?fishery)
hist(fishery$cpue, xlab = "Catch per unit effort")
plot(cpue ~ tide_phase, fishery, pch = 16,
   xlab = "Tide phase", ylab = "Catch per unit effort")
plot(cpue ~ location, fishery, pch = 16,
   xlab = "Location", ylab = "Catch per unit effort")
plot(cpue ~ max_temp, fishery, pch = 16,
   xlab = "Maximum temperature", ylab = "Catch per unit effort")

## Fit a double Box-Cox normal regression model:
fit <- BCSreg(cpue ~ location + tide_phase |
               location + tide_phase + max_temp, fishery)

## Influence measures under case-weights perturbation scheme:
cw <- influence(fit) ## two index plots are shown
str(cw)
}
\references{
Lesaffre, E., and Verbeke, G. (1998). Local influence in linear mixed models.
 \emph{Biometrics}, 570--582.

 Medeiros, R. M. R., and Queiroz, F. F. (2025). Modeling positive continuous data:
    Box-Cox symmetric regression models and their extensions
}
\seealso{
\code{\link{BCSreg}} for parameter estimation in the class of the Box-Cox
    symmetric or zero-adjusted Box-Cox symmetric regression models,
    \code{\link{residuals.BCSreg}} for extracting residuals for a fitted model, and
    \code{\link{plot.BCSreg}} for diagnostic plots.
}
\author{
Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>

Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
}
