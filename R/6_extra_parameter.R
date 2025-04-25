#' @name extra.parameter
#'
#' @title Select the Extra Parameter of a Box-Cox Symmetric Regression
#'
#' @description Estimation of the extra parameter in a Box-Cox regression fit based on the Upsilon
#'     goodness-of-fit statistic and profiled likelihood.
#'
#' @param object an object of class \code{"BCSreg"}, a result of a call to \code{\link{BCSreg}}.
#' @param family a character that specifies the symmetric generating family of the
#'     BCS distribution. The options available are: \code{"ST"}, \code{"PE"},
#'     \code{"SN"}, \code{"HP"}, and \code{"SL"}. The \code{"NO"}, \code{"LOI"},
#'     and \code{"LOII"} generating families do not depend on additional parameters.
#' @param grid grid of values that will be used to evaluate the Upsilon statistics and the profiled
#'     log-likelihood function.
#' @param trace logical; if \code{TRUE}, a summary with the profiled log-likelihood value and the
#'     Upsilon statistics is displayed.
#' @param plot logical; if \code{TRUE}, a graph of the Upsilon statistics evaluated in the
#'     considered grid of values is shown.
#' @param control a list of control arguments specified via \code{\link{BCSreg.control}}.
#' @param ... further arguments passed to \code{\link{BCSreg.control}}.
#'
#' @return An object of class \code{"extra.parameter"}. More specifically, it returns a list in which
#'     each element consists of the fit of the BCS regression with each value of the extra
#'     parameter zeta specified in \code{grid}. In addition, it has the elements \code{"logLik"}
#'     with the vector of log-likelihood values for each fit, \code{"Upsilon"} with the values
#'     of the Upsilon statistics for each fit, and \code{"grid"} with the specified grid of
#'     values.
#'
#'     The \code{print} function summarizes the fits by displaying, for each value in \code{grid},
#'     the log-likelihood value and the Upsilon statistic. The
#'     \code{plot} function returns a graph of the Upsilon statistics, highlighting its
#'     minimum.
#'
#' @references Medeiros, R. M. R., and Queiroz, F. F. (2025). Modeling positive continuous data:
#'     Box-Cox symmetric regression models and their extensions
#'
#' @examples
#' # Examples
#'
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @export
#'
extra.parameter <- function(object, family, grid = seq(1, 30, 2), trace = TRUE, plot = TRUE,
         control = BCSreg.control(...), ...) {

  if(family == "NO" | family == "LOI" | family == "LOII")
    stop("This model does not depend on extra parameters.")

  n <- object$nobs
  fit_update <- lapply(grid, function(zeta) {

    opt <- try(stats::update(object, family = family, zeta = zeta, control = control), silent = TRUE)

    if (trace) {
      cat(
        "\nzeta:", zeta,
        "|",
        "logLik:", if (unique(grepl("Error", opt))) NA else sprintf("%.3f", stats::logLik(opt)),
        "|",
        "Upsilon:", if (unique(grepl("Error", opt))) NA else sprintf("%.3f", summary(opt)$Upsilon)
      )
    }

    opt
  })

  ll <- Upsilon <- vector("numeric", length(grid))
  for (i in 1:length(grid)) {
    ll[i] <- if (unique(grepl("Error", fit_update[[i]]))) NA else stats::logLik(fit_update[[i]])
    Upsilon[i] <- if (unique(grepl("Error", fit_update[[i]]))) NA else summary(fit_update[[i]])$Upsilon
  }



  if (plot) {

    plot(grid, Upsilon, type = "o", pch = 16, cex = 0.6,
         xlab = expression(zeta), ylab = expression(Upsilon(zeta)), las = 1)
    graphics::abline(h = Upsilon[which.min(Upsilon)], lty = 3, col = "grey", lwd = 2)
    graphics::points(c(grid[which.min(Upsilon)], grid[which.min(Upsilon)]),
                     c(Upsilon[which.min(Upsilon)], Upsilon[which.min(Upsilon)]),
                     col = c("dodgerblue", 1), pch = c(16, 1))


  }

  fit_update <- stats::setNames(fit_update, paste0("zeta = ", grid))
  fit_update$logLik <- ll
  fit_update$Upsilon <- Upsilon
  fit_update$grid <- grid
  class(fit_update) <- "extra.parameter"
  fit_update

}


# Print
#' @rdname extra.parameter
#' @param x an object of class \code{"extra.parameter"}.
#'
#' @export
print.extra.parameter <- function(x, ...) {

  #cat("BCS regression fit with mode", x$zeta, "\n")

  grid <- x$grid
  i <- 1
  ll <- Upsilon <- vector("numeric", length(grid))
  for (zeta in grid) {
    ll[i] <- if (unique(grepl("Error", x[[i]]))) NA else x$logLik[i]
    Upsilon[i] <- if (unique(grepl("Error", x[[i]]))) NA else x$Upsilon[i]

    cat(
      "zeta:", zeta,
      "|",
      "logLik:", ll[i],
      "|",
      "Upsilon:", Upsilon[i], "\n"
    )

    i <- i + 1
  }

  cat("\n\nBest value for zeta according to Upsilon:", grid[which.min(Upsilon)],
      "and logLik:", grid[which.max(ll)])

  invisible(x)
}


# Plot
#' @rdname extra.parameter
#' @param x an object of class \code{"extra.parameter"}.
#'
#' @export
plot.extra.parameter <- function(x, ...) {

  grid <- x$grid
  ll <- Upsilon <- vector("numeric", length(grid))
  for (i in 1:length(grid)) {
    ll[i] <- if (unique(grepl("Error", x[[i]]))) NA else x$logLik[i]
    Upsilon[i] <- if (unique(grepl("Error", x[[i]]))) NA else x$Upsilon[i]
  }

  plot(grid, Upsilon, type = "o", pch = 16, cex = 0.6,
       xlab = expression(zeta), ylab = expression(Upsilon(zeta)), las = 1)
  graphics::abline(h = Upsilon[which.min(Upsilon)], lty = 3, col = "grey", lwd = 2)
  graphics::points(c(grid[which.min(Upsilon)], grid[which.min(Upsilon)]),
                   c(Upsilon[which.min(Upsilon)], Upsilon[which.min(Upsilon)]),
                   col = c("dodgerblue", 1), pch = c(16, 1))

  invisible(x)
}
