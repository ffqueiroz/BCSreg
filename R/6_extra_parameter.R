#' @name extra.parameter
#'
#' @title Select the Extra Parameter of a Box-Cox Symmetric Regression Model
#'
#' @description Estimation of the extra parameter in a Box-Cox symmetric or zero-adjusted
#'     Box-Cox symmetric regression model based on the Upsilon goodness-of-fit statistic and
#'     the profile log-likelihood.
#'
#' @param object an object of class \code{"BCSreg"}, resulting from a call to \code{\link{BCSreg}}.
#' @param family a character string specifying the symmetric generating family of the
#'     BCS distribution. The available options are: \code{"ST"}, \code{"PE"},
#'     \code{"SN"}, \code{"HP"}, and \code{"SL"}. The families \code{"NO"}, \code{"LOI"},
#'     and \code{"LOII"} do not depend on an additional parameter.
#' @param grid a numeric vector of positive values at which the Upsilon statistic and
#'     the profile log-likelihood function will be evaluated.
#' @param trace logical; if \code{TRUE}, a summary displaying the profile log-likelihood
#'     values and the Upsilon statistics is shown.
#' @param plot logical; if \code{TRUE}, a plot of the Upsilon statistics evaluated over
#'     the specified grid of values is displayed.
#' @param which numeric; if \code{which = 1}, the plot shows the Upsilon statistic versus \code{zeta};
#'     if \code{which = 2}, the plot shows the profile log-likelihood versus \code{zeta}.
#' @param control a list of control arguments specified via \code{\link{BCSreg.control}}.
#' @param ... further arguments passed to \code{\link{BCSreg.control}}.
#'
#' @return An object of class \code{"extra.parameter"}, which is a list containing the fits of the
#'     BCS regression model for each value of the extra parameter \code{zeta} specified in \code{grid}.
#'     The object also includes:
#'     \itemize{
#'      \item{\code{logLik}: a vector with the log-likelihood values for each fit;}
#'      \item{\code{Upsilon}: a vector with the Upsilon statistic values for each fit;}
#'      \item{\code{grid}: the specified grid of values.}
#'     }
#'
#'     The value of the extra parameter (\code{zeta}) can be selected using two alternative
#'     approaches:
#'     \itemize{
#'      \item{the value that minimizes the Upsilon goodness-of-fit statistic;}
#'      \item{the value that maximizes the log-likelihood.}
#'     }
#'
#'     The \code{print} method summarizes the fits by displaying, for each value in \code{grid},
#'     the corresponding log-likelihood value and Upsilon statistic.
#'     The \code{plot} method returns a graph of the Upsilon statistics, highlighting its minimum.
#'
#' @references
#'
#'     Medeiros, R. M. R., and Queiroz, F. F. (2025). Flexible modeling of nonnegative continuous
#'     data: Box-Cox symmetric regression and its zero-adjusted extension.
#'
#' @examples
#' ## Data set: raycatch (for description, run ?raycatch)
#' hist(raycatch$cpue, xlab = "Catch per unit effort")
#' plot(cpue ~ tide_phase, raycatch, pch = 16,
#'    xlab = "Tide phase", ylab = "Catch per unit effort")
#' plot(cpue ~ location, raycatch, pch = 16,
#'    xlab = "Location", ylab = "Catch per unit effort")
#' plot(cpue ~ max_temp, raycatch, pch = 16,
#'    xlab = "Maximum temperature", ylab = "Catch per unit effort")
#'
#' ## Fit the Box-Cox normal regression as a reference model
#' fit_bcno <- BCSreg(cpue ~ location + tide_phase |
#'                 location + tide_phase + max_temp, raycatch)
#'
#' ## Use the specifications of the reference model to change the distribution
#' ## to, for example, Box-Cox t, and select the value of the extra parameter:
#' select_bct <- extra.parameter(fit_bcno, family = "ST", grid = 1:20)
#'
#' ## Class
#' class(select_bct)
#'
#' ## It is possible to recover the plots:
#' plot(select_bct)
#' plot(select_bct, which = 2)
#'
#' ## and the trace:
#' select_bct
#'
#' ## Selected fit based on the Upsilon statistic
#' fit_bct <- select_bct[["zeta = 19"]]
#' summary(fit_bct)
#' @author Francisco F. de Queiroz <\email{felipeq@ime.usp.br}>
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@ufrn.br}>
#'
#' @export
extra.parameter <- function(object, family, grid = seq(1, 30, 2), trace = TRUE, plot = TRUE,
         control = BCSreg.control(...), ...) {

  if(family == "NO" | family == "LOI" | family == "LOII")
    stop("This model does not depend on extra parameters.")

  n <- object$nobs
  fit_update <- lapply(grid, function(zeta) {

    opt <- try(stats::update(object, family = family, zeta = zeta, control = control), silent = TRUE)

    if (trace) {
      cat(
        "zeta:", zeta,
        "|",
        "logLik:", if (unique(grepl("Error", opt))) NA else sprintf("%.3f", stats::logLik(opt)),
        "|",
        "Upsilon:", if (unique(grepl("Error", opt))) NA else sprintf("%.3f", suppressWarnings(summary(opt)$Upsilon)),
        "\n"
      )
    }

    opt
  })

  ll <- Upsilon <- vector("numeric", length(grid))
  for (i in 1:length(grid)) {
    ll[i] <- if (unique(grepl("Error", fit_update[[i]]))) NA else stats::logLik(fit_update[[i]])
    Upsilon[i] <- if (unique(grepl("Error", fit_update[[i]]))) NA else suppressWarnings(summary(fit_update[[i]])$Upsilon)
  }

  if (trace) {
    cat("\nBest value for zeta according to Upsilon:", grid[which.min(Upsilon)],
        "and Profile log-lik.:", grid[which.max(ll)], "\n")
  }

  if (plot) {

    op <- graphics::par(mfrow = c(1, 2))
    on.exit(graphics::par(op))

    plot(grid, Upsilon, type = "o", pch = 16, cex = 0.6,
         xlab = expression(zeta), ylab = expression(Upsilon(zeta)), las = 1)
    graphics::abline(h = Upsilon[which.min(Upsilon)], lty = 3, col = "grey", lwd = 2)
    graphics::points(c(grid[which.min(Upsilon)], grid[which.min(Upsilon)]),
                     c(Upsilon[which.min(Upsilon)], Upsilon[which.min(Upsilon)]),
                     col = c("dodgerblue", 1), pch = c(16, 1))


    plot(grid, ll, type = "o", pch = 16, cex = 0.6,
         xlab = expression(zeta), ylab = "Profile log-likelihood", las = 1)
    graphics::abline(h = ll[which.max(ll)], lty = 3, col = "grey", lwd = 2)
    graphics::points(c(grid[which.max(ll)], grid[which.max(ll)]),
                     c(ll[which.max(ll)], ll[which.max(ll)]),
                     col = c("dodgerblue", 1), pch = c(16, 1))

  }

  fit_update <- stats::setNames(fit_update, paste0("zeta = ", grid))
  fit_update$logLik <- ll
  fit_update$Upsilon <- Upsilon
  fit_update$grid <- grid
  class(fit_update) <- "extra.parameter"
  invisible(fit_update)

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

  cat("\nBest value for zeta according to Upsilon:", grid[which.min(Upsilon)],
      "and logLik:", grid[which.max(ll)])

  invisible(x)
}


# Plot
#' @rdname extra.parameter
#' @param x an object of class \code{"extra.parameter"}.
#'
#' @export
plot.extra.parameter <- function(x, which = 1, ...) {

  if(!is.numeric(which) || any(which < 1) || any(which > 2))
    stop("`which' must be in 1:2")

  grid <- x$grid
  ll <- Upsilon <- vector("numeric", length(grid))
  for (i in 1:length(grid)) {
    ll[i] <- if (unique(grepl("Error", x[[i]]))) NA else x$logLik[i]
    Upsilon[i] <- if (unique(grepl("Error", x[[i]]))) NA else x$Upsilon[i]
  }

  if (which == 1L) {
    plot(grid, Upsilon, type = "o", pch = 16, cex = 0.6,
         xlab = expression(zeta), ylab = expression(Upsilon(zeta)), las = 1)
    graphics::abline(h = Upsilon[which.min(Upsilon)], lty = 3, col = "grey", lwd = 2)
    graphics::points(c(grid[which.min(Upsilon)], grid[which.min(Upsilon)]),
                     c(Upsilon[which.min(Upsilon)], Upsilon[which.min(Upsilon)]),
                     col = c("dodgerblue", 1), pch = c(16, 1))
  } else if (which == 2L) {
    plot(grid, ll, type = "o", pch = 16, cex = 0.6,
         xlab = expression(zeta), ylab = "Profile log-likelihood", las = 1)
    graphics::abline(h = ll[which.max(ll)], lty = 3, col = "grey", lwd = 2)
    graphics::points(c(grid[which.max(ll)], grid[which.max(ll)]),
                     c(ll[which.max(ll)], ll[which.max(ll)]),
                     col = c("dodgerblue", 1), pch = c(16, 1))
  }

  invisible(x)
}
