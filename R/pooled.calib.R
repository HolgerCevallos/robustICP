#' @title Pooled Conformal Method (Calibration)
#' 
#' @description
#' Computes the calibration for a pooled conformal prediction
#' method that accounts for departures from exchangeability
#' in multilevel functional data.
#' The aim is to detect outlying
#' trajectories at the individual level.
#' The pooled conformal method borrows information from the curves of the other
#' available individuals, which form the calibration set,
#' with the goal of obtaining smaller prediction regions for the target individual.
#' As with the jackknife method, the pooled procedure uses
#' a robust functional PCA based on a sieve procedure.
#' After projecting the trajectories onto an
#' orthonormalized B-spline basis, it applies robust multivariate PCA
#' using \code{\link{DsubLTS}} to the scores of the calibration set.
#' The corresponding standardized residuals are then computed.
#' A weighted empirical distribution of those residuals is used to obtain the
#' \eqn{(1-\alpha)} conformal calibration quantile.
#' The weights reflect the distributional similarity of each calibration individual
#' to the target individual based on the maximum mean discrepancy (MMD).
#' The resulting quantile can then be used
#' to flag new trajectories from the target individual
#' as regular or potentially anomalous.
#' This approach is suitable for multilevel,
#' temporally structured (functional) data, such as those
#' commonly encountered in free-living digital health studies.
#' 
#' @param x A list where each element contains trajectories from a different
#' individual. Each element must be a numeric matrix whose rows correspond
#' to observed trajectories for that individual and whose columns correspond
#' to the time points specified by \code{mesh}. Missing values are not allowed.
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
#' \item{mesh}{Same as the input parameter \code{mesh}.}
#' \item{dimension.Bspline}{Same as the input parameter \code{dimension.Bspline}.}
#' \item{mu.target.calib}{Robust location estimate of the target subject's
#' scores obtained by applying DsubLTS to those scores.}
#' \item{quantile.calib}{Calibration threshold that determines the prediction region.}
#' \item{mu.calib}{Robust location estimate of the calibration scores obtained by applying DsubLTS to the pooled matrix of centered calibration scores.}
#' \item{B.calib}{Estimated first \eqn{q} orthonormal loadings of the calibration scores obtained by applying DsubLTS to the pooled matrix of centered calibration scores.}
#' \item{var.weights}{Importance weights for each of the
#' \code{dimension.Bspline} components, used to compute the standardized
#' residuals. The weights are normalized to sum to one.}
#' \item{MMD.weights}{Weights used to construct the weighted
#' empirical distribution of standardized residuals from which
#' the calibration threshold is derived. These weights are computed
#' from inverse MMD distances. Each subject's weight is repeated
#' according to its number of calibration trajectories.}
#' \item{fitted.values.pooled}{Orthogonal projections of the
#' pooled centered calibration scores onto the robust subspace estimated
#' from these scores using DsubLTS.}
#' \item{std.residuals.pooled}{Standardized residuals of the
#' pooled centered calibration scores relative to their corresponding
#' fitted scores in \code{fitted.values.pooled}.}
#' \item{fitted.values.targetsubj}{Orthogonal projections of the
#' available centered scores of the target subject onto the
#' robust subspace estimated during calibration using DsubLTS.}
#' \item{std.residuals}{Standardized residuals of the
#' available centered scores of the target subject relative to their
#' corresponding fitted scores in \code{fitted.values.targetsubj}.}
#' \item{predRegIndices}{Indices of the available trajectories of the
#' target subject included in the prediction region.}
#' \item{predRegCurves}{Available trajectories of the target subject
#' included in the prediction region.}
#' \item{coverage}{Empirical proportion of the available trajectories
#' of the target subject included in the prediction region.}
#'
#' @details
#' The element of \code{x} indexed by \code{indexTargetSubject}
#' corresponds to the target subject, while the remaining elements
#' correspond to the calibration subjects for the pooled method.
#' Our pooled conformal prediction method uses robust functional PCA 
#' based on a sieve procedure to obtain
#' a robust representation of the calibration data, which serves as a reference
#' for detecting potentially anomalous curves for the target subject.
#' The sieve procedure consists of
#' projecting the trajectories onto an
#' orthonormalized B-spline basis of dimension \code{dimension.Bspline},
#' yielding a matrix of scores for each subject, with each row
#' corresponding to a trajectory.
#' This allows the functional PCA problem to be treated in a finite-dimensional
#' space of dimension \code{dimension.Bspline}.
#' The scores of each subject are first centered using its corresponding robust
#' location estimate obtained with \code{\link{DsubLTS}}, a robust and fast PCA method for multivariate data.
#' The robust representation of the calibration set is then obtained
#' by applying \code{\link{DsubLTS}} to the matrix of pooled, centered scores
#' from the calibration subjects.
#' Standardized residuals are computed as weighted Euclidean distances between
#' each centered calibration score and its projection onto the estimated robust subspace.
#' The weights estimate the contribution of each of the
#' \code{dimension.Bspline} components in explaining the variance of the
#' calibration scores. They are normalized to sum to one.
#' This weighting scheme allows components that explain a larger proportion
#' of the variance to receive greater importance in the computation of the
#' standardized residuals.
#' To obtain the \eqn{(1-\alpha)} conformal calibration quantile,
#' a weighted empirical distribution of these residuals is used.
#' The weights at this stage are based on inverse MMD distances,
#' where MMD distances quantify distributional dissimilarities
#' between calibration subjects and the target subject.
#' Consequently, subjects that are more distributionally similar to the
#' target subject receive greater weight when determining the
#' conformal calibration quantile.
#' MMD distances are computed by \code{\link{.mmd_gaussian}} 
#' using a Gaussian kernel with a diagonal covariance matrix.
#' Its diagonal entries equal the squared median of the 
#' pairwise Euclidean distances
#' computed from the pooled score data of all subjects.
#' Note that the subspace dimension used in \code{\link{DsubLTS}}
#' is given by the parameter \code{q}.
#' The estimated location vector and loading matrix
#' from the robust PCA fit on the calibration data,
#' together with the calibration quantile, are used to detect
#' potentially anomalous trajectories from the target subject.
#' The function identifies potentially anomalous trajectories among
#' the available trajectories of the target subject. It also returns
#' the calibration threshold and the estimated quantities required
#' to evaluate new trajectories from the target subject.
#'
#' @seealso
#' \code{\link{pooled.pred}},
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
#' fit.pooled.method$quantile.calib
#' fit.pooled.method$var.weights
#' fit.pooled.method$std.residuals.pooled
#'                                  
#' @export
pooled.calib <- function(x,
                         mesh, dimension.Bspline,
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

# Function to compute the $(1-\alpha)$-quantile of
# the weighted empirical distribution of the calibration residuals.
umbral.train <- function(x, w, p){
  sw <- sum(w)
  o <- order(x)
  xs <- x[o]
  ws <- w[o]
  cw <- cumsum(ws) / sw
  xs[ which(cw >= p)[1] ]
}

# Function to compute standardized residuals as a weighted Euclidean distance
distancescore <- function(y,z,w){
  out <- sqrt(sum(((y-z)^2) * w))
  return(out)
} 

list.scores.todos <- list() 
n.subjects <- length(x)


# Robust centering of projection data for each subject using the corresponding robust PCA center estimate
for (subject in seq_len(n.subjects)) {
  x.obs <- x[[subject]] # Extract subject functional data
  
  # Project subject data onto the B-spline basis
  y <- x.obs %*% base.estim.Bspline
  
  # Estimate the robust center using DsubLTS
  y.DsubLTS <- robustPCA::DsubLTS(
    x = y,
    q = q,
    alpha = alpha.DsubLTS,
    scale = FALSE
  )
  
  mu.y.DsubLTS <- y.DsubLTS$mu
  
  # Save the robust center of the target subject
  if (subject == indexTargetSubject) {
    mu.target <- mu.y.DsubLTS
  }
  
  # Robust centering
  list.scores.todos[[subject]] <- sweep(y, 2, mu.y.DsubLTS, FUN = "-")
}
 
# Computing MMD matrix.
# Row index of the MMD matrix indicates target subject.
# Each row of the MMD matrix contains information of similarity of each of the other subjects
# with respect to that target subject.

dist_matrix <- matrix(0, nrow = n.subjects, ncol = n.subjects)

# Building the reference sequence
for (i.dmm in 1:n.subjects) {
  list_reference <- list.scores.todos
  tmp <- list_reference[[i.dmm]]
  list_reference[[i.dmm]] <- list_reference[[n.subjects]]
  list_reference[[n.subjects]] <- tmp
  
  matriz_A <- do.call(rbind, list_reference)
  
  # Building the swaped sequence
  for (j.dmm in 1:n.subjects) {
    
    list_j <- list_reference
    tmp <- list_j[[j.dmm]]
    list_j[[j.dmm]] <- list_j[[n.subjects]]
    list_j[[n.subjects]] <- tmp
    
    matriz_B <- do.call(rbind, list_j)
    
    xy <- rbind(t(matriz_A), t(matriz_B))
    
    # MMD distance with Gaussian kernel with a diagonal covariance matrix,
    # variance set to the squared median of the pairwise Euclidean distances 
    # between the feature vectors in the pooled data
    dists <- dist(xy)  
    median_dist <- median(as.numeric(dists))  
    var_kernel <- diag(rep(median_dist^2, ncol(t(matriz_A))))
  
    # Computing biased MMD distances
    dist_matrix[i.dmm, j.dmm] <- .mmd_gaussian(
                                                t(matriz_A), 
                                                t(matriz_B),
                                                var = var_kernel
                                                )               
  }
}

# Excluding subject indexTargetSubject for calibration
list.sin.indexTargetSubject <- list.scores.todos[-indexTargetSubject]

# Calibration set
y.sin.indexTargetSubject <- do.call(rbind, list.sin.indexTargetSubject)

################## Fitting the method on the calibration set ####################################

# The calibration set includes all curves from the other available individuals
n <- nrow(y.sin.indexTargetSubject)
p <- ncol(y.sin.indexTargetSubject)

# Variance estimation for the standardization of residuals
var.explained <- sapply(1:p, function(k) {
  y.DsubLTS <- robustPCA::DsubLTS(x = y.sin.indexTargetSubject, q = k, alpha = alpha.DsubLTS, scale=FALSE)
  mu.y.DsubLTS <- y.DsubLTS$mu
  yc <- y.sin.indexTargetSubject - matrix(mu.y.DsubLTS, nrow = n, ncol = p, byrow = TRUE)
  bon <- y.DsubLTS$B
  y.rulo.DsubLTS <- ((yc %*% bon) %*% t(bon)) + 
    matrix(mu.y.DsubLTS, nrow = n, ncol = p, byrow = TRUE)
  
  1 - (
    mean(rowSums((y.sin.indexTargetSubject - y.rulo.DsubLTS)^2)) /
      mean(rowSums((y.sin.indexTargetSubject - mu.y.DsubLTS)^2))
  )
})

# Standardization weights
pesos <- pmax(0, c(var.explained[1], diff(var.explained)))
pesos <- pesos / sum(pesos)

# Calibration using robust PCA (DsubLTS algorithm)
fit.DsubLTS <- robustPCA::DsubLTS(x = y.sin.indexTargetSubject, q = q, alpha = alpha.DsubLTS, scale=FALSE)
yc <- y.sin.indexTargetSubject - matrix(fit.DsubLTS$mu, nrow = n, ncol = p, byrow = TRUE)
bon <- fit.DsubLTS$B
DsubLTS.prediction <- ((yc %*% bon) %*% t(bon)) + 
  matrix(fit.DsubLTS$mu, nrow = n, ncol = p, byrow = TRUE)

# Residuals on the calibration set

# Creating names according to subject index
ids_originales <- seq_along(list.scores.todos) 
# Extracting index of the target individual
ids_filtrados <- ids_originales[-indexTargetSubject]  
# Repeating indexes of the other individuals according to the number of curves each of them has
paciente_id <- rep(ids_filtrados, times = sapply(list.sin.indexTargetSubject, nrow)) 

# Computing standardized residuals on the calibration set
resid.pred <- rep(0, n)
for(i.vir in 1:n){
  resid.pred[i.vir] <- distancescore(y.sin.indexTargetSubject[i.vir, ], 
                                     DsubLTS.prediction[i.vir,], 
                                     w = pesos)
  names(resid.pred)[i.vir] <- paste0("R", paciente_id[i.vir]) # assigning subject indices 
  # as names to the residuals
} 

# Computing weights for the other individuals relative to the target subject
# via the inverse-MMD rule

# The last swap leaves the sequence identical to the reference one (no-swap configuration), 
# as a result, the last column of the MMD distance matrix is equal to zero and is therefore removed.
D <- dist_matrix[ , -n.subjects] # MMD distance matrix
denominator <- colSums(1 / D)
W.subject.indexTargetSubject <- (1/D[indexTargetSubject,])/denominator

# Replicating weights according to the number of curves each of the other subjects have
# and assigning subject indices as names to the weights vector
nombres.R <- names(resid.pred) 
indices.R <- as.integer(sub("R", "", nombres.R)) 
W.rep <- W.subject.indexTargetSubject[indices.R]
indice.Rtarget <- which(nombres.R == paste0("R", n.subjects))
# Weights in the location of target subject corresponds to last
# subject and therefore are placed at the end of the list.
W.rep[indice.Rtarget] <- W.subject.indexTargetSubject[indexTargetSubject] 
names(W.rep) <- nombres.R

################## Predictions on the target subject’s projection data ####################################

# Computing predictions for the target subject’s scores data
# using the robust PCA fit from the calibration set
y.nuevo <- list.scores.todos[[indexTargetSubject]]
n.y.nuevo <- dim(y.nuevo)[1]
yc <- y.nuevo - matrix(fit.DsubLTS$mu, nrow = n.y.nuevo, ncol = p, byrow = TRUE)
y.new.rulo.DsubLTS <- ((yc %*% bon) %*% t(bon)) +
  matrix(fit.DsubLTS$mu, nrow = n.y.nuevo, ncol = p, byrow = TRUE)

# Computing standardized residuals on the projection data of the target subject
resid.y.new.pred <- sqrt(rowSums(((y.nuevo - y.new.rulo.DsubLTS)^2) *
                                   matrix(pesos, nrow = n.y.nuevo,
                                          ncol = p, byrow = TRUE)))

# Computing the $(1-\alpha)$-quantile of
# the weighted empirical distribution of the calibration residuals.
  cuantil <- umbral.train(unname(resid.pred), unname(W.rep), p=(1-alpha))

# Constructing the prediction region on the curves of the target subject
  index.good.new.curves <- which(resid.y.new.pred <= cuantil)

# Computing empirical coverage
coverage.new.y <- length(index.good.new.curves)/n.y.nuevo

res <- list(
  mesh = mesh,
  dimension.Bspline = dimension.Bspline,
  mu.target.calib = mu.target,
  quantile.calib = cuantil,
  mu.calib = fit.DsubLTS$mu,
  B.calib = bon,
  var.weights = pesos,
  MMD.weights = W.rep,
  fitted.values.pooled = DsubLTS.prediction,
  std.residuals.pooled = resid.pred,
  fitted.values.targetsubj = y.new.rulo.DsubLTS,
  std.residuals = resid.y.new.pred,
  predRegIndices = index.good.new.curves,
  predRegCurves = y.nuevo[index.good.new.curves, , drop = FALSE],
  coverage = coverage.new.y)

class(res) <- "pooled.calib"
return(res)
}
