################################################################################
#              SIMPLEX REGRESSION - LIKELIHOOD RATIO TEST                      #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: Likelihood ratio test for comparing nested simplex models       #
################################################################################

# ==============================================================================
# LIKELIHOOD RATIO TEST
# ==============================================================================

#' @title Likelihood Ratio Test for Simplex Regression
#' @description Performs a likelihood ratio test to compare two nested simplex
#' regression models.
#'
#' @param model An object of class \code{"simplexregression"} representing the
#' restricted (null) model
#' @param model_amp An object of class \code{"simplexregression"} representing the
#' augmented (alternative) model
#'
#' @details
#' The likelihood ratio test compares two nested models: a restricted model
#' (under the null hypothesis) and an augmented model (under the alternative
#' hypothesis). The test statistic is:
#'
#' \deqn{LR = 2(\ell_1 - \ell_0)}
#'
#' where \eqn{\ell_0} and \eqn{\ell_1} are the log-likelihoods of the
#' restricted and augmented models, respectively.
#'
#' Under the null hypothesis, the test statistic follows a chi-squared
#' distribution with degrees of freedom equal to the difference in the number
#' of parameters between the two models.
#'
#' The models must be nested, meaning that the restricted model is a special
#' case of the augmented model obtained by setting some parameters to specific
#' values (typically zero).
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: The likelihood ratio test statistic
#'   \item \code{parameter}: Degrees of freedom
#'   \item \code{p.value}: The p-value of the test
#'   \item \code{method}: Description of the test
#'   \item \code{data.name}: Names of both model objects being compared
#' }
#'
#' @references
#' Wilks, S. S. (1938). The large-sample distribution of the likelihood ratio
#' for testing composite hypotheses. \emph{The Annals of Mathematical Statistics},
#' 9(1), 60-62.
#'
#' @examples
#' \dontrun{
#' # Fit restricted model (null hypothesis)
#' model0 <- simplexreg(y ~ x1, data = mydata)
#'
#' # Fit augmented model (alternative hypothesis)
#' model1 <- simplexreg(y ~ x1 + x2 + x3, data = mydata)
#'
#' # Perform likelihood ratio test
#' lrtest.simplexreg(model0, model1)
#' }
#'
#' @importFrom stats pchisq
#'
#' @export
lrtest.simplexreg <- function(model, model_amp) {
  METHOD = "Likelihood Ratio test"
  DNAME  = deparse(substitute(modelh0))
  DNAME  = paste(DNAME, "vs", deparse(substitute(modelh1)))

  lH0 <- model$loglik
  lH1 <- model_amp$loglik
  LR <- 2*(lH1 - lH0)
  gl = model$df.residual - model_amp$df.residual
  PVAL = pchisq(LR, gl, lower.tail= F)

  names(gl) = "df"
  names(LR) = "LR"

  RVAL <- list(statistic = LR, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"

  return(RVAL)
}
