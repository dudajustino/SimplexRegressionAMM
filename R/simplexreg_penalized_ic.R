################################################################################
#               PENALIZED INFORMATION CRITERIA FOR SIMPLEX REGRESSION          #
################################################################################

#' @title Penalized Information Criteria for Parametric Simplex Regression
#' @description Computes AIC, BIC, and HQ information criteria with a penalty term
#' depending on the estimated parameter \eqn{\lambda} of the parametric link
#' function ("plogit1" or "plogit2").
#'
#' @param model An object of class \code{"simplexregression"} fitted with a
#' parametric link ("plogit1" or "plogit2").
#' @param kappa A numeric constant controlling the penalization strength
#'   (default is 0.1).
#'
#' @details
#' The penalized criteria are computed as:
#' \deqn{AIC_c = -2 \ell + (2 + c \, |\log(\lambda)|)(n - df_{res})}
#' \deqn{BIC_c = -2 \ell + (\log n + c \, |\log(\lambda)|)(n - df_{res})}
#' \deqn{HQ_c  = -2 \ell + (2 \log(\log n) + c \, |\log(\lambda)|)(n - df_{res})}
#' where:
#' \itemize{
#'   \item \eqn{\ell} is the log-likelihood at the maximum;
#'   \item \eqn{\lambda} is the power parameter of the parametric link function;
#'   \item \eqn{df_{res}} is the number of residual degrees of freedom.
#' }
#'
#' These penalized versions add a smooth penalty on the magnitude of
#' \eqn{\lambda}, encouraging simpler link structures.
#'
#' @return A named vector with components \code{AICc}, \code{BICc}, and \code{HQc}.
#'
#' @examples
#' \dontrun{
#' # Fit two models with parametric links
#' fit1 <- simplexreg(y ~ x1 + x2, link = "plogit1")
#' fit2 <- simplexreg(y ~ x1 + x2, link = "plogit2")
#'
#' # Compute penalized criteria
#' penalized_ic(fit1)
#' penalized_ic(fit2, kappa = 0.2)
#' }
#'
#' @export
penalized_ic <- function(model, kappa = 0.1) {

  if (is.null(model$lambda.fv)) {
    stop("Model does not have a parametric link (lambda.fv missing).")
  }

  n <- length(model$fitted.values)
  r <- n - model$df.residual
  lambda <- model$lambda.fv
  ll <- model$loglik
  c <- kappa

  # penalized terms
  aic_c <- -2 * ll + (2 + c * abs(log(lambda))) * r
  bic_c <- -2 * ll + (log(n) + c * abs(log(lambda))) * r
  hq_c  <- -2 * ll + (2 * log(log(n)) + c * abs(log(lambda))) * r

  out <- c(AICc = aic_c, BICc = bic_c, HQc = hq_c)
  return(out)
}
