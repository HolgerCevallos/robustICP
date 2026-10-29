#' @title Jackknife Conformal Method (Prediction)
#' 
#' @description
#' Performs subject-specific conformal prediction on 
#' one or more new trajectories from a target individual
#' using a previously calibrated jackknife conformal model.
#' The prediction procedure requires the jackknife conformal fit
#' and calibration threshold obtained from a calibration set,
#' both returned in a \code{\link{jackknife.calib}} object.
#' These are used to determine
#' whether new trajectories can be considered regular or potentially anomalous.
#' The jackknife conformal method aims to achieve the nominal coverage
#' \eqn{(1-\alpha)} for future repeated trajectories of a subject,
#' with empirical results supporting this coverage property.
#' 
#' @param xnew A numeric vector, matrix, or data frame containing one or more
#' new (future) trajectories for the target subject. For a matrix or data frame,
#' rows correspond to trajectories and columns correspond to the time points
#' specified by \code{mesh} during calibration. For a vector, its length
#' must equal the number of time points in \code{mesh}.
#' @param object An object inheriting from class \code{"jackknife.calib"}.
#'
#' @return 
#' \item{quantile.calib}{Calibration threshold that determines the prediction region.}
#' \item{var.weights}{Importance weights for each of the
#' \code{dimension.Bspline} components during calibration.
#' They are used to compute the standardized residuals and 
#' are normalized to sum to one.}
#' \item{fitted.values.new}{Orthogonal projections of the
#' new scores of the target subject onto the
#' robust subspace estimated during calibration using DsubLTS.}
#' \item{std.residuals.new}{Standardized residuals of the
#' new scores of the target subject relative to their
#' corresponding fitted scores in \code{fitted.values.new}.}
#' \item{predRegIndices}{Indices of the new trajectories of the
#' target subject included in the prediction region.}
#' \item{predRegCurves}{New trajectories of the target subject
#' included in the prediction region.}
#' \item{coverage}{Empirical proportion of new trajectories
#' of the target subject included in the prediction region.}
#'
#' @details
#' \code{jackknife.pred} determines whether one or more new trajectories
#' from a target individual belong to the calibrated
#' conformal prediction region and are thus considered regular
#' or otherwise potentially anomalous.
#' The function applies a previously calibrated jackknife conformal model
#' returned by \code{\link{jackknife.calib}} to new trajectories
#' from the target individual.
#' The calibration of the jackknife model is based on previously observed
#' trajectories from the target individual, which constitute the calibration set.
#' The new trajectories are first projected onto the same orthonormalized
#' B-spline basis used during calibration to obtain corresponding scores.
#' A prediction for each new score is then obtained by projecting it
#' onto the robust principal subspace estimated by DsubLTS
#' from the calibration set.
#' Next, standardized residuals are computed 
#' as weighted Euclidean distances between
#' each new score and its prediction.
#' Note that the weights correspond to the importance weights
#' obtained during calibration.
#' The standardized residuals are then compared with the jackknife calibration quantile.
#' A new trajectory is included in the conformal prediction region when its
#' residual is smaller than or equal to this calibration threshold.
#' The function returns the calibration threshold, the indices and
#' trajectories included in the prediction region, and an
#' empirical coverage estimate, along with additional prediction results.
#' 
#' @seealso
#' \code{\link{jackknife.calib}},
#' \code{\link{DsubLTS}}
#'
#' @examples
#' # Generate the calibration data
#' calibration.train <- robustICP::calibration.data(Npatients = 10,
#'                                       sigma = 1,
#'                                       tm = (0:20)/100,
#'                                       curves.per.patient = 50,
#'                                       seed = 1234)
#' 
#' # Application of the Jackknife Conformal Method (Calibration)
#' fit.jackknife.method <- jackknife.calib(x = calibration.train$Y,
#'                                         mesh = calibration.train$timepoints,
#'                                         dimension.Bspline = 20,
#'                                         q = 2,
#'                                         indexTargetSubject = 1,
#'                                         alpha = 0.5)
#'
#' # Application of the Jackknife Conformal Method (Prediction)
#' pred.jackknife.method <- jackknife.pred(xnew = test.data(indexTargetSubject = 1,
#'                                                          curves.per.patient = 50,
#'                                                          sigma=calibration.train$sigma,
#'                                                          Nt=calibration.train$Nt,
#'                                                          seed=calibration.train$seed),
#'                                         object = fit.jackknife.method)
#'
#' # Calibration quantile
#' fit.jackknife.method$quantile.calib
#'
#' # Prediction region indices
#' pred.jackknife.method$predRegIndices
#'
#' # Prediction region curves
#' pred.jackknife.method$predRegCurves
#'
#' # Empirical coverage
#' pred.jackknife.method$coverage
#'
#' @export
jackknife.pred <- function(xnew, object){
  if (!inherits(object, "jackknife.calib")) {
    stop("The 'object' argument must be a model object returned by jackknife.calib().", call. = FALSE)
  }
  
  # xnew validation
  if (is.vector(xnew)) {
    if (!is.numeric(xnew)) {
      stop("'xnew' must be a numeric vector, matrix, or data frame.", call. = FALSE)
    }
    xnew <- matrix(xnew, nrow = 1)
    
  } else if (is.matrix(xnew)) {
    if (!is.numeric(xnew)) {
      stop("'xnew' must contain only numeric values.", call. = FALSE)
    }
    
  } else if (is.data.frame(xnew)) {
    if (!all(vapply(xnew, is.numeric, logical(1)))) {
      stop("'xnew' must contain only numeric values.", call. = FALSE)
    }
    xnew <- data.matrix(xnew)
    
  } else {
    stop("'xnew' must be a numeric vector, matrix, or data frame.", call. = FALSE)
  }
  
  #Extract calibration measurements from model
  mesh <- object$mesh
  dimension.Bspline <- object$dimension.Bspline
  quantile.calib <- object$quantile.calib
  mu.calib <- object$mu.calib
  B.calib <- object$B.calib
  var.weights <- object$var.weights
  
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
  
  # Function to compute standardized residuals as a weighted Euclidean distance
  distancescore <- function(y,z,w){
    out <- sqrt(sum(((y-z)^2) * w))
    return(out)
  }
  
  # Check compatibility between xnew and the B-spline basis
  if (ncol(xnew) != nrow(base.estim.Bspline)) {
    stop("'xnew' and the B-spline basis must have the same number of evaluation points.", call. = FALSE)
  }
  
  # Projecting new data of the given target subject onto the B-spline basis
  y <- xnew %*% base.estim.Bspline
  n.y.nuevo <- nrow(y)
  p <- ncol(y)
  
  # Computing predictions on the projection data for the new trajectories of the given target subject 
  # using the robust PCA fit from the calibration set
  
  yc <- y - matrix(mu.calib, nrow = n.y.nuevo, ncol = p, byrow = TRUE)
  y.new.rulo.DsubLTS <- ((yc %*% B.calib) %*% t(B.calib)) + 
    matrix(mu.calib, nrow = n.y.nuevo, ncol = p, byrow = TRUE)
  
  resid.y.new.pred <- sqrt(rowSums(((y - y.new.rulo.DsubLTS)^2) * 
                                     matrix(var.weights, nrow = n.y.nuevo, 
                                            ncol = p, byrow = TRUE)))
  
  # Constructing the prediction region on the new curves of the target subject
  index.good.new.curves <- which(resid.y.new.pred <= quantile.calib)
  
  res = list(
    quantile.calib = quantile.calib,
    var.weights = var.weights,
    fitted.values.new = y.new.rulo.DsubLTS,
    std.residuals.new = resid.y.new.pred,
    predRegIndices = index.good.new.curves,
    predRegCurves = xnew[index.good.new.curves,  , drop = FALSE])

  multiple.curves <- nrow(xnew) > 1
  
  # Computing empirical coverage when exist multiple trajectories in the subject
  if (multiple.curves) {
    coverage.new.y <- length(index.good.new.curves) / n.y.nuevo
    res$coverage <- coverage.new.y
  }else {
    res$coverage <- as.numeric(length(index.good.new.curves) > 0)
  }
  return(res)
}
