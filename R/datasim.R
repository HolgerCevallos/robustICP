# Design 1: Generation of simulated multilevel functional data
lambda.k1.values <- function(k){
  lambda.k1 <- 0.5^(k-1)
  return(lambda.k1)
}

lambda.l2.values <- function(l){
  lambda.l2 <- 0.5^(l-1)
  return(lambda.l2)
}

xik.values <- function(lambda.k1){
  xik <- rnorm(1,0,lambda.k1)
  return(xik)
}

zetaijl.values <- function(lambda.l2){
  zetaijl <- rnorm(1,0,lambda.l2)
  return(zetaijl)
}

phik1.values <- function(t){
  phik1 <- list(sqrt(2)*sin(2*pi*t),
                sqrt(2)*cos(2*pi*t),
                sqrt(2)*sin(4*pi*t),
                sqrt(2)*cos(4*pi*t))
  return(phik1)
}

phil2.values <- function(t){
  phil2 <- list(sqrt(2)*sin(6*pi*t),
                sqrt(2)*cos(6*pi*t),
                sqrt(2)*sin(8*pi*t),
                sqrt(2)*cos(8*pi*t))
  return(phil2)
}

#' @title Simulated Calibration Dataset
#'
#' @description Simulated multilevel functional data generated
#' following the simulation framework described in
#' Di et al. (2009). The dataset contains simulated data
#' for ten subjects, each with 50 trajectories observed
#' at 21 time points as default settings.
#'
#' @param Npatients Integer specifying the number of subjects to simulate. Defaults to 10.
#' @param sigma Numeric value specifying the standard deviation of the 
#' measurement error added to each trajectory. Defaults to 1.
#' @param tm Numeric vector specifying the time points at which the
#' functional trajectories are evaluated. Defaults to `(0:20)/100`,
#' corresponding to 21 time points between 0 and 0.2.
#' @param curves.per.patient Integer specifying the number of trajectories
#' generated for each subject. Defaults to 50.
#' @param seed Integer used to initialize the random number generator
#' during data generation. Defaults to 1234.
#' 
#' @return A list containing the following components:
#' \item{Yij.todos}{A list of length \code{Npatients}. Each element is a
#' numeric matrix with \code{curves.per.patient} rows and
#' \code{length(tm)} columns, containing the simulated trajectories
#' for one subject.}
#' \item{tm}{Numeric vector of observation time points supplied through
#' the \code{tm} argument.}
#' \item{sigma}{Standard deviation of the measurement error used to
#' generate the trajectories.}
#' \item{Nt}{Integer giving the number of observation time points,
#' calculated as \code{length(tm)}.}
#' \item{seed}{Integer supplied through the \code{seed} argument.}
#' 
#' @references
#' Di, C. Z., Crainiceanu, C. M., Caffo, B. S., and Punjabi, N. M. (2009).
#' Multilevel functional principal component analysis.
#' \emph{The Annals of Applied Statistics}, 3(1), 458--488.
#' \url{https://doi.org/10.1214/08-AOAS206}
#
#' @keywords datasets
#' 
#' @export
calibration.data <- function(Npatients = 10,
                             sigma = 1,
                             tm = (0:20)/100,
                             curves.per.patient = 50,
                             seed = 1234) {

Yij.todos <- vector("list", Npatients) 
Nt <- length(tm)

for (s in 1:Npatients){
  
  fixed.subject.effect <- rep(0,Nt)
  visit.effect <- matrix(0, nrow=curves.per.patient, ncol=Nt)
  Yij <- matrix(0, nrow=curves.per.patient, ncol=Nt)
  
  
  for(t in 1:Nt){
    for(k in 1:4){
      lambda.k1 <- lambda.k1.values(k)
      set.seed(seed+40+s)
      fixed.xk <- xik.values(lambda.k1)
      fixed.subject.effect[t] <- fixed.subject.effect[t] + 
        fixed.xk*phik1.values(tm[t])[[k]]
    }
  }
  
  for(j in 1:curves.per.patient){
    for(t in 1:Nt){
      for(l in 1:4){
        lambda.l2 <- lambda.l2.values(l)
        set.seed(seed+s+j)
        zetaijl <- zetaijl.values(lambda.l2)
        visit.effect[j,t] <- visit.effect[j,t] + 
          zetaijl*phil2.values(tm[t])[[l]]
      }
    }
  }
  
  for(j in 1:curves.per.patient){
    for(t in 1:Nt){
      Yij[j,t] <- fixed.subject.effect[t] + visit.effect[j,t] +
        rnorm(1,0,sigma)
    }
  }
  Yij.todos[[s]] <- Yij
}

return(list(
  Y = Yij.todos,
  timepoints = tm,
  sigma = sigma,
  Nt = Nt,
  seed = seed
))
}