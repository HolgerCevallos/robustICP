#' @title Jackknife Conformal Method (Calibration)
#' 
#' @description Computes the calibration for a jackknife conformal
#' prediction method for detecting potentially anomalous trajectories
#' at the individual level. This approach is suitable for multilevel,
#' temporally structured (functional) data, such as those
#' commonly encountered in free-living digital health studies.
#' The proposed jackknife conformal method
#' relies exclusively on within-individual data
#' and uses a robust functional PCA based on a sieve procedure.
#' The latter consists of projecting each trajectory 
#' of the target subject onto an
#' orthonormalized B-spline basis. 
#' Next, leave-one-out robust PCA fits are obtained
#' from the corresponding scores using \code{\link{DsubLTS}},
#' a robust and fast PCA method
#' for multivariate data.
#' Standardized jackknife residuals are then computed
#' for each left-out score, and the empirical distribution
#' of these residuals is used to obtain the \eqn{(1-\alpha)}
#' conformal calibration quantile. This quantile can then be used
#' to flag new trajectories from the same target individual
#' as regular or potentially anomalous.
#'
#' @param x A numeric matrix or data frame whose rows correspond to the observed
#' trajectories of the target subject and whose columns correspond to the
#' time points specified by \code{mesh}. Missing values are not allowed.
#' @param dimension.Bspline Number of B-spline basis functions to be used 
#' in the sieve procedure.
#' @param mesh Grid points at which the B-spline basis is evaluated.
#' These also correspond to the observation time points of the trajectories.
#' @param q Chosen dimension of the functional PCA.
#' @param indexTargetSubject Index of the list element in \code{x}
#' corresponding to the target subject.
#' @param alpha Miscoverage level of the prediction region. 
#' Thus, \eqn{(1-\alpha)} is the nominal coverage level.
#' @param alpha.DsubLTS Trimming fraction used in the \code{\link{DsubLTS}} algorithm. This parameter
#' determines the fraction of outliers the algorithm should resist. Allowed values are
#' between \eqn{0} and \eqn{0.5}, with a default value of \eqn{0.5}.
#' 
#' @return 
#' \item{mesh}{same as input parameter \code{mesh}.}
#' \item{dimension.Bspline}{same as input parameter \code{dimension.Bspline}.}
#' \item{quantile.calib}{Calibration threshold that determines the prediction region.}
#' \item{mu.calib}{Robust location estimate of the 
#' calibration scores obtained by applying DsubLTS 
#' to the target subject's calibration scores.}
#' \item{B.calib}{Estimated first \eqn{q} orthonormal loadings of the 
#' calibration scores obtained by applying DsubLTS to the
#' target subject's calibration scores.}
#' \item{var.weights}{Importance weights for each of the
#' \code{dimension.Bspline} components, used to compute the standardized
#' residuals. The weights are normalized to sum to one.}
#' \item{fitted.values.jackknife}{Orthogonal projections of the
#' left-out scores onto their corresponding
#' leave-one-out robust subspaces estimated with DsubLTS.}
#' \item{std.residuals.jackknife}{Standardized residuals of the
#' left-out scores relative to their corresponding fitted scores
#' in \code{fitted.values.jackknife}.}
#' \item{fitted.values.targetsubj}{Orthogonal projections of the
#' available (calibration) scores of the target subject onto its
#' robust subspace estimated using DsubLTS.}
#' \item{std.residuals}{Standardized residuals of the
#' available (calibration) scores of the target subject relative to their
#' corresponding fitted scores in \code{fitted.values.targetsubj}.}
#' \item{predRegIndices}{Indices of the available (calibration)
#' trajectories of the target subject included in the
#' prediction region.}
#' \item{predRegCurves}{Available (calibration) trajectories
#' of the target subject included in the 
#' prediction region.}
#' \item{coverage}{Empirical proportion of the available 
#' (calibration) trajectories
#' of the target subject included in the prediction region.}
#' 
#' @details
#' Our jackknife conformal prediction method uses robust functional PCA
#' based on a sieve procedure to obtain a robust representation of
#' the target subject's calibration data, which serves as a reference
#' for detecting potentially anomalous future curves.
#' The sieve procedure consists
#' of projecting the trajectories of the target subject onto an
#' orthonormalized B-spline basis of dimension \code{dimension.Bspline},
#' yielding a matrix of scores, with each row corresponding to a trajectory.
#' This allows the functional PCA problem to be treated 
#' in a finite-dimensional space of dimension \code{dimension.Bspline}. 
#' The conformal calibration quantile is obtained using a jackknife
#' scheme on the scores, whereas the robust PCA fit used for prediction
#' is obtained by applying \code{\link{DsubLTS}}, 
#' a fast and robust PCA method for multivariate data,
#' to the complete calibration score matrix.
#' In particular, to obtain the calibration quantile, \code{\link{DsubLTS}}
#' is applied to the scores matrix following a jackknife scheme.
#' For each left-out score, a standardized residual is computed as the
#' weighted Euclidean distance between that score and its projection
#' onto the corresponding leave-one-out robust subspace estimate.
#' The weights estimate the contribution of each of the
#' \code{dimension.Bspline} components in explaining the variance of
#' the target subject's scores. They are normalized to sum to one.
#' This weighting scheme therefore gives greater importance to components
#' that explain a larger proportion of the variance when computing
#' the standardized residuals.
#' The conformal calibration quantile is then calculated as the
#' \eqn{(1-\alpha)}-quantile of the empirical distribution of the
#' resulting standardized jackknife residuals, where \code{alpha}
#' specifies the conformal miscoverage level.
#' To obtain the robust PCA fit used for prediction,
#' \code{\link{DsubLTS}} is applied to the complete
#' calibration score matrix of the target subject.
#' The estimated location vector and loading matrix from this fit,
#' together with the calibration quantile, are used to detect
#' potentially anomalous future trajectories from the target subject.
#' Note that the subspace dimension used in DsubLTS
#' is given by the parameter \code{q}.
#' 
#' @seealso
#' \code{\link{jackknife.pred}},
#' \code{\link{pooled.calib}},
#' \code{\link{pooled.pred}},
#' \code{\link{DsubLTS}}
#'
#' @examples
#' # Generate the calibration data
#' calibration.data <- robustICP::calibration.data(Npatients = 10,
#'                                       sigma = 1,
#'                                       tm = (0:20)/100,
#'                                       curves.per.patient = 50,
#'                                       seed = 1234)
#' 
#' # Application of the Jackknife Conformal Method (Calibration)
#' fit.jackknife.method <- jackknife.calib(x = calibration.data[["Yij.todos"]],
#'                                         mesh = calibration.data[["tm"]],
#'                                         dimension.Bspline = 20,
#'                                         q = 2,
#'                                         indexTargetSubject = 1,
#'                                         alpha = 0.5)
#'
#' # Inspect the calibration quantile, residuals and prediction region indices
#' fit.jackknife.method$quantile.calib
#' fit.jackknife.method$var.weights
#' fit.jackknife.method$std.residuals.jackknife
#' 
#' @export
jackknife.calib <- function(x,
                            mesh,
                            dimension.Bspline,
                            q,
                            indexTargetSubject,
                            alpha,
                            alpha.DsubLTS=0.5){
  
  # Input validation: x
  if (!is.list(x) || length(x) == 0L) {
    stop("'x' must be a non-empty list of numeric matrices.", call. = FALSE)
  }
  
  if (!all(vapply(x, is.matrix, logical(1)))) {
    stop("All elements of 'x' must be matrices.", call. = FALSE)
  }
  
  if (!all(vapply(x, is.numeric, logical(1)))) {
    stop("All elements of 'x' must be numeric matrices.", call. = FALSE)
  }
  
  # Validate indexTargetSubject
  if (!is.numeric(indexTargetSubject) ||
      length(indexTargetSubject) != 1L ||
      !is.finite(indexTargetSubject) ||
      indexTargetSubject != floor(indexTargetSubject) ||
      indexTargetSubject < 1 ||
      indexTargetSubject > length(x)) {
    stop( "'indexTargetSubject' must be an integer between 1 and length(x).", call. = FALSE )}  
  
# R functions to generate B-spline basis
  Bspline.manual <- function(mesh, dimension.Bspline){
    if (dimension.Bspline == 5)
      part <- 2
    else
      part <- dimension.Bspline/2
    aa <- ((-part):part)/part
    dimension.Bsplinenots <- (aa+1)*1.0/2
    base.estim.disc <- mgcv::cSplineDes(mesh,dimension.Bsplinenots)
    base.ortonorm.disc <- qr.Q(qr(base.estim.disc))
    return(base.ortonorm.disc)
  }
  
# Generating B-spline basis functions
  base.estim.Bspline <- Bspline.manual(mesh, dimension.Bspline) 
  
# Check compatibility between x and the B-spline basis
  if (ncol(x[[indexTargetSubject]]) != nrow(base.estim.Bspline)) {
    stop("The target subject and the B-spline basis must have the same number of evaluation points.", call. = FALSE)
  }

# Function to compute standardized residuals as a weighted Euclidean distance
  distancescore <- function(y,z,w){
    out <- sqrt(sum(((y-z)^2) * w))
    return(out)
  }
  
  x<-data.matrix(x[[indexTargetSubject]])
  #x<-data.matrix(x)
  
# Projecting subject data onto the B-spline basis
  y <- x %*% base.estim.Bspline
  n <- nrow(y)
  p <- ncol(y)
  
# Variance estimation for the standardization of residuals
var.explained <- sapply(1:p, function(k) {
  y.DsubLTS <- robustPCA::DsubLTS(x = y, q = k, alpha = alpha.DsubLTS, scale=FALSE)
    mu.y.DsubLTS <- y.DsubLTS$mu
    bon <- y.DsubLTS$B
    yc <- y - matrix(mu.y.DsubLTS, nrow = n, ncol = p, byrow = TRUE)
    y.rulo.DsubLTS <- ((yc %*% bon) %*% t(bon)) + 
      matrix(mu.y.DsubLTS, nrow = n, ncol = p, byrow = TRUE)

    1 - (mean(rowSums((y- y.rulo.DsubLTS)^2)) /
           mean(rowSums((y - mu.y.DsubLTS)^2)))
  })
  
# Standardization weights
  pesos <- pmax(0, c(var.explained[1], diff(var.explained)))
  pesos <- pesos / sum(pesos);pesos
  

  fitted.values.jackknife <- matrix(NA_real_, nrow = n, ncol = p)
  resid.jackknife.pred <- numeric(n)
  
  for (i in seq_len(n)) {
    jackknife.fit.DsubLTS <- robustPCA::DsubLTS(
      x = y[-i, , drop = FALSE],
      q = q,
      alpha = alpha.DsubLTS,
      scale = FALSE
    )
    yc <- y[i, ] - jackknife.fit.DsubLTS$mu
    bon <- jackknife.fit.DsubLTS$B
    jackknife.DsubLTS.i.prediction <-
      ((yc %*% bon) %*% t(bon)) + jackknife.fit.DsubLTS$mu
    fitted.values.jackknife[i, ] <- jackknife.DsubLTS.i.prediction
    resid.jackknife.pred[i] <- distancescore(
      y[i, ], jackknife.DsubLTS.i.prediction, w = pesos
    )
  }
  
# Computing the $(1-\alpha)$-quantile of
# the empirical distribution of the calibration residuals.
  cuantil <- quantile(resid.jackknife.pred, probs= (1-alpha))
  
# Calibration using robust PCA (DsubLTS algorithm)
  y.DsubLTS <- robustPCA::DsubLTS(x=y, q=q, alpha=alpha.DsubLTS, scale=FALSE)
  mu.y.DsubLTS <- y.DsubLTS$mu
  bon <- y.DsubLTS$B
  
  # Prediction on the observed curves using the fitted robust PCA
  
  yc <- y - matrix(
    mu.y.DsubLTS,
    nrow = n,
    ncol = p,
    byrow = TRUE
  )
  
  y.pred <- (yc %*% bon) %*% t(bon) +
    matrix(
      mu.y.DsubLTS,
      nrow = n,
      ncol = p,
      byrow = TRUE
    )
  
  # Standardized residuals on the observed curves
  
  resid.pred <- sqrt(
    rowSums(
      (y - y.pred)^2 *
        matrix(pesos, nrow = n, ncol = p, byrow = TRUE)
    )
  )
  
  # Identify curves included in the prediction region
  
  index.good.curves <- which(resid.pred <= as.numeric(cuantil))
  
  # Empirical proportion of included curves
  
  coverage <- length(index.good.curves) / n
  
  res <- list(
    mesh = mesh,
    dimension.Bspline = dimension.Bspline,
    quantile.calib = cuantil,
    mu.calib = mu.y.DsubLTS,
    B.calib = bon,
    var.weights = pesos,
    fitted.values.jackknife = fitted.values.jackknife,
    std.residuals.jackknife = resid.jackknife.pred,
    fitted.values.targetsubj = y.pred,
    std.residuals = resid.pred,
    predRegIndices = index.good.curves,
    predRegCurves = y[index.good.curves, , drop = FALSE],
    coverage = coverage
    )
  
  class(res) <- "jackknife.calib"
  return(res)
}
