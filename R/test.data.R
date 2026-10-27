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

#' @title Simulated Test Dataset
#'
#' @description Simulated functional observations generated for the subject
#' with index \code{indexTargetSubject}, which is the target
#' subject of interest. The data were generated following the
#' simulation framework of Di et al. (2009).
#' The calibration data for this subject can be found in
#' \code{calibration.data[[indexTargetSubject]]}.
#' The dataset contains 50 new simulated trajectories
#' observed at 21 time points for that subject.
#' 
#' @param indexTargetSubject Integer identifying the target subject for which
#' the test trajectories are generated. Defaults to 1.
#' @param sigma Numeric value specifying the standard deviation of the
#' measurement error added to each trajectory. Defaults to 1.
#' @param Nt Integer specifying the number of time points at which the
#' trajectories are observed. Defaults to 21.
#' @param seed Integer used to initialize the random number generator
#' during the generation of the subject-specific and visit-specific effects.
#' Defaults to 1234.
#' @param curves.per.patient Integer specifying the number of test trajectories
#' to generate for the target subject. Defaults to 50.
#'
#' @return A numeric matrix with \code{curves.per.patient} rows and
#' \code{Nt} columns. Each row represents a simulated functional trajectory
#' for the target subject, and each column corresponds to a time point.
#' 
#' @references
#' Di, C. Z., Crainiceanu, C. M., Caffo, B. S., and Punjabi, N. M.
#' (2009). Multilevel functional principal component analysis.
#' \emph{The Annals of Applied Statistics}, 3(1), 458--488.
#' \url{https://doi.org/10.1214/08-AOAS206}
#'
#' @keywords datasets
#' 
#' @export
test.data <- function(indexTargetSubject = 1,
                      sigma=1,
                      Nt=21,
                      seed=1234,
                      curves.per.patient = 50) {

s <- as.integer(indexTargetSubject)
#sigma <- 1
#Nt <- 21
#curves.per.patient <- 50
fixed.subject.effect <- rep(0,Nt)
visit.effect <- matrix(0, nrow=curves.per.patient, ncol=Nt)
Yij <- matrix(0, nrow=curves.per.patient, ncol=Nt)
tm <- (0:20)/100

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
      set.seed(seed+100+s+j)
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
Yij
}
