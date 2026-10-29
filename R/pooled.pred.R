#' @title Pooled Conformal Method (Prediction)
#'
#' @description
#' Performs subject-specific conformal prediction on
#' one or more new trajectories
#' from a target individual using a previously calibrated pooled conformal model.
#' The prediction procedure requires the pooled conformal fit
#' and calibration threshold obtained from a calibration set,
#' both returned in a \code{\link{pooled.calib}} object.
#' These are used to determine
#' whether new trajectories can be considered regular or potentially anomalous.
#' The pooled conformal method aims to achieve the nominal coverage
#' \eqn{(1-\alpha)} for future repeated trajectories of a subject,
#' with empirical results supporting this coverage property.
#' 
#' @param xnew A numeric vector, matrix, or data frame containing one or more
#' new (future) trajectories for the target subject. For a matrix or data frame,
#' rows correspond to trajectories and columns correspond to the time points
#' specified by \code{mesh} during calibration. For a vector, its length
#' must equal the number of time points in \code{mesh}.
#' @param object An object inheriting from class \code{"pooled.calib"}.
#' 
#' @return
#' \item{quantile.calib}{Calibration threshold that determines the prediction region.}
#' \item{var.weights}{Importance weights for each of the
#' \code{dimension.Bspline} components during calibration.
#' They are used to compute the standardized residuals and 
#' are normalized to sum to one.}
#' \item{fitted.values.new}{Orthogonal projections of the
#' new centered scores of the target subject onto the
#' robust subspace estimated during calibration using DsubLTS.}
#' \item{std.residuals.new}{Standardized residuals of the
#' new centered scores of the target subject relative to their
#' corresponding fitted scores in \code{fitted.values.new}.}
#' \item{predRegIndices}{Indices of the new trajectories of the
#' target subject included in the prediction region.}
#' \item{predRegCurves}{New trajectories of the target subject
#' included in the prediction region.}
#' \item{coverage}{Empirical proportion of new trajectories
#' of the target subject included in the prediction region.}
#'
#' @details
#' \code{pooled.pred} determines whether one or more new trajectories
#' from a target individual belong to the calibrated
#' conformal prediction region and are thus considered regular
#' or otherwise potentially anomalous.
#' The function applies a previously calibrated pooled conformal model
#' returned by \code{\link{pooled.calib}} to new trajectories
#' from the target individual.
#' The calibration of the pooled model is based on trajectories
#' from other individuals, accounting for their distributional
#' similarity to the trajectories of the target individual during calibration.
#' The new trajectories of the target individual are first projected 
#' onto the same orthonormalized
#' B-spline basis used during calibration to obtain corresponding scores. 
#' These scores are robustly centered using the target subject's
#' location estimate obtained with DsubLTS during calibration.
#' A prediction for each centered score is then obtained by projecting it
#' onto the robust principal subspace estimated by DsubLTS
#' from the calibration set.
#' Next, standardized residuals are computed 
#' as weighted Euclidean distances between
#' each centered score and its prediction.
#' Note that the weights correspond to the importance weights
#' obtained during calibration.
#' The standardized residuals are then compared with the pooled calibration quantile.
#' A new trajectory is included in the conformal prediction region when its
#' residual is smaller than or equal to this calibration threshold.
#' The function returns the calibration threshold, the indices
#' and trajectories included in the prediction region, and an
#' empirical coverage estimate, along with additional prediction results.
#'
#' @seealso
#' \code{\link{pooled.calib}},
#' \code{\link{jackknife.calib}},
#' \code{\link{jackknife.pred}},
#' \code{\link{DsubLTS}},
#' \code{\link{.mmd_gaussian}}
#'
#' @examples
#' # Generate the calibration data
#' calibration.train <- robustICP::calibration.data(Npatients = 10,
#'                                       sigma = 1,
#'                                       tm = (0:20)/100,
#'                                       curves.per.patient = 50,
#'                                       seed = 1234)
#' 
#' # Application of the calibration for a pooled conformal prediction method
#' fit.pooled.method <- pooled.calib(x = calibration.train$Y,
#'                                         mesh = calibration.train$timepoints,
#'                                         dimension.Bspline = 20,
#'                                         q = 2,
#'                                         indexTargetSubject = 1,
#'                                         alpha = 0.5)
#'
#' # Application of the Pooled Conformal Method (Prediction)
#' pred.pooled.method <- pooled.pred(xnew = test.data(indexTargetSubject = 1,
#'                                                    curves.per.patient = 50,
#'                                                    sigma=calibration.train$sigma,
#'                                                    Nt=calibration.train$Nt,
#'                                                    seed=calibration.train$seed),
#'                                   object = fit.pooled.method)
#'
#' # Calibration quantile
#' pred.pooled.method$quantile.calib
#'
#' # Prediction region indices
#' pred.pooled.method$predRegIndices
#'
#' # Prediction region curves
#' pred.pooled.method$predRegCurves
#'
#' # Empirical coverage
#' pred.pooled.method$coverage
#' 
#' @export
pooled.pred <- function(xnew, object)
{
  if (!inherits(object, "pooled.calib")) {
    stop("The 'object' argument must be a model object returned by pooled.calib().", call. = FALSE )}

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
  
  # Extracting calibration information
  mesh <- object$mesh
  dimension.Bspline <- object$dimension.Bspline
  quantile.calib <- object$quantile.calib
  mu.calib <- object$mu.calib
  B.calib <- object$B.calib
  var.weights <- object$var.weights
  mu.target.calib <- object$mu.target.calib
  
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
  
  # Check compatibility between xnew and the B-spline basis
  if (ncol(xnew) != nrow(base.estim.Bspline)) {
    stop("'xnew' and the B-spline basis must have the same number of evaluation points.", call. = FALSE)
  }
  
  # Projecting new trajectories onto the B-spline basis to obtain scores
  y.new <- xnew %*% base.estim.Bspline
  
  n.y.nuevo <- nrow(y.new)
  p <- ncol(y.new)
   
  # Centering scores of target subject using the robust location
  # estimated with DsubLTS during calibration.
  y.nuevo <- y.new - matrix(mu.target.calib, nrow = n.y.nuevo, ncol = p, byrow = TRUE)
  
  # Projecting centered scores onto the robust principal subspace estimated by DsubLTS
  # from the calibration set.
  yc <- y.nuevo - matrix(mu.calib, nrow = n.y.nuevo, ncol = p, byrow = TRUE)
  y.new.rulo.DsubLTS <- ((yc %*% B.calib) %*% t(B.calib)) + 
    matrix(mu.calib, nrow = n.y.nuevo, ncol = p, byrow = TRUE)
  
  #Computing standardized residuals on the centered scores of the target subject
  resid.y.new.pred <- sqrt(rowSums(((y.nuevo - y.new.rulo.DsubLTS)^2) *
                                     matrix(var.weights, nrow = n.y.nuevo,
                                            ncol = p, byrow = TRUE)))
  
  #Constructing the prediction region on the new curves of the target subject
    index.good.new.curves <- which(resid.y.new.pred <= quantile.calib)
    
    # Creating output object
    res <- list(
      quantile.calib = quantile.calib,
      var.weights = var.weights,
      fitted.values.new = y.new.rulo.DsubLTS,
      std.residuals.new = resid.y.new.pred,
      predRegIndices = index.good.new.curves,
      predRegCurves = xnew[index.good.new.curves,  , drop = FALSE])
    
    multiple.curves <- n.y.nuevo > 1
    
    # Computing empirical coverage when exist multiple trajectories in the subject
    if (multiple.curves) {
      coverage.new.y <- length(index.good.new.curves) / n.y.nuevo
      res$coverage <- coverage.new.y
    }else {
      res$coverage <- as.numeric(length(index.good.new.curves) > 0)
    }
    return(res)
}
