################################################################################
#                     SIMPLEX REGRESSION - RESET TEST                          #
# Author: Maria Eduarda da Cruz Justino and Francisco Cribari-Neto             #
# Date: 2025-11-08                                                             #
# Description: RESET (Regression Equation Specification Error Test) for        #
#              testing functional form misspecification                        #
################################################################################

# ==============================================================================
# RESET TEST
# ==============================================================================

#' @title RESET Test for Simplex Regression
#' @description Performs the RESET (Regression Equation Specification Error Test)
#' to detect functional form misspecification in simplex regression models.
#'
#' @param model An object of class \code{"simplexregression"}
#' @param dispersion Logical. If \code{TRUE}, includes the squared linear predictor
#' in the dispersion submodel as well. Default is \code{TRUE}
#'
#' @details
#' The RESET test augments the original model by adding the squared linear
#' predictor as an additional covariate. Under the null hypothesis of correct
#' functional form, this additional term should not be significant.
#'
#' The test statistic follows a chi-squared distribution with degrees of freedom
#' equal to the difference in the number of parameters between the augmented
#' and original models.
#'
#' If \code{dispersion = TRUE}, the squared linear predictor is added to both
#' the mean and dispersion submodels. If \code{FALSE}, it is only added to the
#' mean submodel.
#'
#' @return An object of class \code{"htest"} containing:
#' \itemize{
#'   \item \code{statistic}: The likelihood ratio test statistic
#'   \item \code{parameter}: Degrees of freedom
#'   \item \code{p.value}: The p-value of the test
#'   \item \code{method}: Description of the test
#'   \item \code{data.name}: Name of the model object
#' }
#'
#' @references
#' Ramsey, J. B. (1969). Tests for specification errors in classical linear
#' least-squares regression analysis. \emph{Journal of the Royal Statistical
#' Society: Series B}, 31(2), 350-371.
#'
#' @examples
#' \dontrun{
#' # Fit a simplex regression model
#' model <- simplexreg(y ~ x1 + x2, data = mydata)
#'
#' # Perform RESET test
#' reset.simplexregression(model)
#'
#' # RESET test only for mean submodel
#' reset.simplexregression(model, dispersion = FALSE)
#' }
#'
#' @importFrom stats pchisq
#'
#' @export
reset_simplexreg <- function(model, dispersion = TRUE){
  METHOD = "RESET: Regression Equation Specification Error Test"
  DNAME  = deparse(substitute(model))

  y <- as.vector(model$y)
  mu_lp_squared <- model$mu.lp^2

  x <- cbind(model$mu.x, mu_lp_squared)
  if(dispersion == TRUE){
    z <- cbind(model$sigma2.x, mu_lp_squared)
  } else {
    z <- model$sigma2.x
  }

  modelh1 <- simplexreg_fit(y, x[,-1, drop=FALSE], z[,-1, drop=FALSE], diag = 0,
                               link.mu = model$mu.link, link.sigma2 = model$sigma2.link,
                               x_names = model$x_names, z_names = model$z_names)

  lH0 <- model$loglik
  lH1 <- modelh1$loglik
  LR <- 2*(lH1 - lH0)
  gl = model$df.residual - modelh1$df.residual
  PVAL = pchisq(LR, gl, lower.tail= F)

  names(gl) = "df"
  names(LR) = "LR"

  RVAL <- list(statistic = LR, parameter = gl, p.value = PVAL,
               method = METHOD, data.name = DNAME)
  class(RVAL) <- "htest"

  return(RVAL)
}
