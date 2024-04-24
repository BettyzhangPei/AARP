#'Assessing the Performance of Dietary Assessment Instruments
#' @details
#'
#' A package for estimating both of the correlations of a dietary assessment instrument (Q) and nutrient intake and attenuation factors and the correlations of nutrient residual.
#' The model is established as the measurement error model to evaluate the data collected by a food-frequency questionnaire (FFQ) using two non-consecutive 24-hour dietary recalls (24HRs) as the reference.
#' We employ the moment estimation method to the complete longitudinal data structure. Moreover, our method can provide the estimates of fixed effects including intake-related biases in the FFQ and overall biases of true intake for the FFQ and 24HR, respectively in the model.
#' Notably, the correlations with true intake and attenuation factors for energy-adjusted should be interpreted on the scale to which each particular nutrient was transformed.
#'
#'
#' @param data data has to be complete data.
#' @param Q.N.1: An assessment instrument for nutrient (N) at time 1.
#' @param Q.N.2: An assessment instrument for nutrient (N) at time 2.
#' @param Q.E.1: An assessment instrument for total energy (E) at time 1.
#' @param Q.E.2: An assessment instrument for total energy (E) at time 2.
#' @param F.N.1: A reference instrument for nutrient (N) at time 1.
#' @param F.N.2: A reference instrument for nutrient (N) at time 2.
#' @param F.E.1: A reference instrument for total energy (E) at time 1.
#' @param F.E.2: A reference instrument for total energy (E) at time 2.
#'
#'
#'
#' @return
#' \item{Rho.N}{correlation coefficient between reported and true intake of nutrient N.}
#' \item{Rho.E}{correlation coefficient between reported and true intake of total energy E.}
#' \item{Rho.R}{correlation coefficient between reported and true intake residuals.}
#' \item{Lambda.N}{attenuation factor of nutrient N. }
#' \item{Lambda.E}{attenuation factor of total energy E.}
#' \item{Lambda.R}{attenuation factor of nutrient residual R.}
#' \item{est.covariance.matrix}{estimated covariance matrix for observed data}
#' \item{Beta.N1}{intake-related bias of nutrient N in the FFQ}
#' \item{Beta.E1}{intake-related bias of total energy E in the FFQ}
#' \item{Mu.T.N}{mean of ture intake of nutrient N}
#' \item{Mu.T.E}{mean of ture intake of total energy E}
#' \item{Mu.N0.2}{overall bias independent of true intake of nutrient N for the 24HR}
#' \item{Mu.E0.2}{overall bias independent of true intake of total energy E for the 24HR}
#' \item{Beta.N0.1}{overall bias independent of true intake of nutrient N for the FFQ in the first measurement}
#' \item{Beta.N0.2}{overall bias independent of true intake of total energy E for the FFQ in the first measurement}
#' \item{Beta.E0.1}{overall bias independent of true intake of nutrient N for the FFQ in the first measurement}
#' \item{Beta.E0.2}{overall bias independent of true intake of total energy E for the FFQ in the first measurement}
#'
#'
#'
#'
#'
#' @importFrom stats sd
#' @examples
#' ## example using random numbers
#'
#' library(MASS)
#' set.seed(1)
#' # number of observed individuals:
#' I<- 100
#' # fixed effects in the measurement error model:
#' # beta = (beta_N0_1, beta_N0_2, beta_N1, beta_E0_1, beta_E0_2, beta_E1)
#' beta<-  c(0.4, 0.5, 0.3,0.32,0.3,0.35)
#' # mu_1 = (mu_T_N, mu_T_E, mu_N0_2, mu_E0_2)
#' mu_1<-  c(0.3,0.32,0.4,0.35)
#' # mu = first moments:  the mean vector of D_i
#' mu<- rep(0,8)
#' mu[1]<- beta[1] + beta[3]*mu_1[1]
#' mu[2]<- beta[2] + beta[3]*mu_1[1]
#' mu[3]<- beta[4] + beta[6]*mu_1[2]
#' mu[4]<- beta[5] + beta[6]*mu_1[2]
#' mu[5]<- mu_1[1]
#' mu[6]<- mu_1[3] + mu_1[1]
#' mu[7]<- mu_1[2]
#' mu[8]<- mu_1[4] + mu_1[2]
#' # second moments:
#' # theta = (theta1, theta2, ..., theta15)
#' # we notice that theta1, theta7, theta11, theta12, theta14, theta15 are variance components, then they should be set as positive numbers.
#' theta<- rep(0,15)
#' # the setting of theta:
#' theta<- c(0.5, 0, 0, 0, 0.175 , 0, 0.58, 0, 0, 0.232, 0.5, 0.35, 0 , 0.55 ,0.3)
#' theta_set<- c(theta[1], theta[2], theta[3], theta[4], theta[5], theta[5], theta[6], theta[6], theta[2], theta[1], theta[4], theta[3], theta[5], theta[5], theta[6], theta[6],
#'              theta[3], theta[4], theta[7], theta[8], theta[9], theta[9], theta[10], theta[10], theta[4], theta[3], theta[8], theta[7], theta[9], theta[9], theta[10], theta[10],
#'              theta[5], theta[5], theta[9], theta[9], theta[11], theta[12], theta[13],theta[13],
#'              theta[5], theta[5], theta[9], theta[9], theta[12], theta[11], theta[13],theta[13],
#'              theta[6], theta[6], theta[10], theta[10], theta[13], theta[13], theta[14], theta[15], theta[6], theta[6],theta[10], theta[10], theta[13], theta[13], theta[15], theta[14])
#' # variance covariance matrix of D_i
# Sigma<- matrix(theta_set, 8, 8)
#' # the observed data:
#' D<- matrix(0, nrow=8 , ncol=I)
#' D<- t(as.matrix(mvrnorm(I, mu, Sigma)))
#' Data<- data.frame(t(D))
#' colnames(Data)<- c("FFQ.N.1","FFQ.N.2", "FFQ.E.1", "FFQ.E.2", "Ref.N.1", "Ref.N.2","Ref.E.1", "Ref.E.2")
#' ex1<- NME(Q.N.1=Data$FFQ.N.1, Q.N.2=Data$FFQ.N.2, Q.E.1= Data$FFQ.E.1, Q.E.2= Data$FFQ.E.2, F.N.1= Data$Ref.N.1, F.N.2= Data$Ref.N.2, F.E.1= Data$Ref.E.1, F.E.2= Data$Ref.E.2)
#' print(ex1)
#'
#' @references Thompson, F. E., Kipnis, V., Midthune, D., Freedman, L.S., Carroll, R.J., Subar, A.F., Brown, C.C., Butcher, M.S.,  Mouw, T., Leitzmann, M. and Schatzkin, A. (2008) Performance of a food-frequency questionnaire in the us NIH--AARP (National Institutes of Health--American Association of Retired Persons) Diet and Health Study. \emph{Public Health Nutrition}, 11, 183-195.
#' @export
NME<- function(Q.N.1, Q.N.2, Q.E.1, Q.E.2, F.N.1, F.N.2, F.E.1, F.E.2)
{

  data<- data.frame(Q.N.1, Q.N.2, Q.E.1, Q.E.2, F.N.1, F.N.2, F.E.1, F.E.2 )
  data<- as.matrix(data)


  # Observed Data: D with nrow=8, ncol=I
  # if(!is.matrix(data)) warning("Forcing the data to be the matrix.")

  D <- t(data)

  if(dim(D)[1]!=8) warning("The dimension of the data matrix should be 8xI where I is the total number of individuals.")

  # Total number of individuals:
  I <- dim(D)[2]


  est.theta<- rep(0,15)

  est.theta[1] <- (sd(D[1,])^2 + sd(D[2,])^2)/2

  est.theta[2] <- sum((D[1,] - mean(D[1,]))*(D[2,] - mean(D[2,])) )/(I-1)
  est.theta[3] <- sum( (D[1,] - mean(D[1,]))*(D[3,] - mean(D[3,])) +  (D[2,] - mean(D[2,]))*(D[4,] - mean(D[4,])) )/ (2*(I-1))

  est.theta[4] <-   sum( (D[1,] - mean(D[1,]))*(D[4,] - mean(D[4,])) +  (D[2,] - mean(D[2,]))*(D[3,] - mean(D[3,])) )/ (2*(I-1))

  est.theta[5]<-   sum( (D[5,] - mean(D[5,]))*(D[1,] - mean(D[1,])) +  (D[6,] - mean(D[6,]))*(D[1,] - mean(D[1,])) +  (D[5,] - mean(D[5,]))*(D[2,] - mean(D[2,])) +  (D[6,] - mean(D[6,]))*(D[2,] - mean(D[2,]))  )/ (4*(I-1))


  est.theta[6] <- sum( (D[7,] - mean(D[7,]))*(D[1,] - mean(D[1,])) +  (D[8,] - mean(D[8,]))*(D[1,] - mean(D[1,])) +  (D[7,] - mean(D[7,]))*(D[2,] - mean(D[2,])) +  (D[8,] - mean(D[8,]))*(D[2,] - mean(D[2,]))  )/ (4*(I-1))

  est.theta[7] <- (sd(D[3,])^2 + sd(D[4,])^2)/2

  est.theta[8] <- sum((D[3,] - mean(D[3,]))*(D[4,] - mean(D[4,]))) /(I-1)

  est.theta[9]<- sum( (D[5,] - mean(D[5,]))*(D[3,] - mean(D[3,])) +  (D[6,] - mean(D[6,]))*(D[3,] - mean(D[3,])) +  (D[5,] - mean(D[5,]))*(D[4,] - mean(D[4,])) +  (D[6,] - mean(D[6,]))*(D[4,] - mean(D[4,]))  )/ (4*(I-1))

  est.theta[10]<-  sum( (D[7,] - mean(D[7,]))*(D[3,] - mean(D[3,])) +  (D[8,] - mean(D[8,]))*(D[3,] - mean(D[3,])) +  (D[7,] - mean(D[7,]))*(D[4,] - mean(D[4,])) +  (D[8,] - mean(D[8,]))*(D[4,] - mean(D[4,]))  )/ (4*(I-1))

  est.theta[11] <- (sd(D[5,])^2 + sd(D[6,])^2)/2


  est.theta[12]<-  (sum( ((D[5,]- mean(D[5,])) * (D[6,]- mean(D[6,])) )  )  / ((I-1)))

  est.theta[13] <- sum( (D[5,] - mean(D[5,]))*(D[7,] - mean(D[7,])) +  (D[5,] - mean(D[5,]))*(D[8,] - mean(D[8,])) +  (D[6,] - mean(D[6,]))*(D[7,] - mean(D[7,])) +  (D[6,] - mean(D[6,]))*(D[8,] - mean(D[8,]))  )/ (4*(I-1))

  est.theta[14] <-  (sd(D[7,])^2 + sd(D[8,])^2)/2

  est.theta[15] <-  (sum( ((D[7,]- mean(D[7,])) * (D[8,]- mean(D[8,])) )) /((I-1)))


  est.theta.set<- c(est.theta[1], est.theta[2], est.theta[3], est.theta[4], est.theta[5], est.theta[5], est.theta[6], est.theta[6], est.theta[2], est.theta[1], est.theta[4], est.theta[3], est.theta[5], est.theta[5], est.theta[6], est.theta[6],
                    est.theta[3], est.theta[4], est.theta[7], est.theta[8], est.theta[9], est.theta[9], est.theta[10], est.theta[10], est.theta[4], est.theta[3], est.theta[8], est.theta[7], est.theta[9], est.theta[9], est.theta[10], est.theta[10],
                    est.theta[5], est.theta[5], est.theta[9], est.theta[9], est.theta[11], est.theta[12], est.theta[13],est.theta[13],
                    est.theta[5], est.theta[5], est.theta[9], est.theta[9], est.theta[12], est.theta[11], est.theta[13], est.theta[13],
                    est.theta[6], est.theta[6], est.theta[10], est.theta[10], est.theta[13], est.theta[13], est.theta[14], est.theta[15], est.theta[6], est.theta[6], est.theta[10], est.theta[10], est.theta[13], est.theta[13], est.theta[15], est.theta[14])



  # estimated covariance matrix
  est.cov.mat <- matrix(est.theta.set , nrow= 8, ncol=8)

  # the estimates of parameters (correlations and attenuation factors):
  est.coeff<- rep(0,6)

  # correlation coefficient for nutrient N:
  est.coeff[1]<- est.theta[5]/sqrt(est.theta[1]* est.theta[12])
  # attenuation factor for nutrient N:
  est.coeff[2]<- est.theta[5]/est.theta[1]

  # correlation coefficient for nutrient E:
  est.coeff[3]<- est.theta[10]/sqrt(est.theta[7]*est.theta[15])
  # attenuation factor for nutrient E:
  est.coeff[4]<- est.theta[10]/est.theta[7]


  # correlation coefficient for residual R:
  est.coeff[5]<- (est.theta[5] * est.theta[7]* est.theta[15] - est.theta[3] * est.theta[9]* est.theta[15]  - est.theta[6] * est.theta[7]* est.theta[13] + est.theta[3] * est.theta[10]* est.theta[13]) / sqrt(est.theta[7]* est.theta[15]*(est.theta[12]* est.theta[15]- (est.theta[13])^2)*(est.theta[1]*est.theta[7] - (est.theta[3])^2))
  # attenuation factor for residual R:
  est.coeff[6]<- (est.theta[5] * est.theta[7]* est.theta[15] - est.theta[3] * est.theta[9]* est.theta[15]  - est.theta[6] * est.theta[7]* est.theta[13] + est.theta[3] * est.theta[10]* est.theta[13]) / (est.theta[15]*(est.theta[1]* est.theta[7]- (est.theta[3])^2))

  est.coeff<- data.frame(est.coeff)
  rownames(est.coeff) <- c("Rho.N", "Lambda.N", "Rho.E",  "Lambda.E", "Rho.R",  "Lambda.R")



  # other unknown fixed parameters in the model:
  est.fix.par<- rep(0,10)

  # estimated beta_N1:
  est.fix.par[1]<- est.theta[5]/est.theta[15]
  # estimated beta_E1:
  est.fix.par[2]<- est.theta[10]/est.theta[15]
  # estimated mu_T_N:
  est.fix.par[3]<- mean(D[5,])
  # estimated mu_T_E:
  est.fix.par[4]<- mean(D[7,])
  # estimated mu_N0_2:
  est.fix.par[5]<- mean(D[6,]) - mean(D[5,])
  # estimated mu_E0_2:
  est.fix.par[6]<- mean(D[8,]) - mean(D[7,])
  # estimated beta_N0_1:
  est.fix.par[7]<- mean(D[1,]) - est.fix.par[1]*est.fix.par[3]
  # estimated beta_N0_2:
  est.fix.par[8]<- mean(D[2,]) - est.fix.par[1]*est.fix.par[3]
  # estimated beta_E0_1:
  est.fix.par[9]<- mean(D[3,]) - est.fix.par[2]*est.fix.par[4]
  # estimated beta_E0_2:
  est.fix.par[10]<- mean(D[4,]) -est.fix.par[2]*est.fix.par[4]

  est.fix.par<- data.frame(est.fix.par)
  rownames(est.fix.par) <- c("Beta.N1", "Beta.E1", "Mu.T.N",  "Mu.T.E", "Mu.N0.2", "Mu.E0.2", "Beta.N0.1", "Beta.N0.2", "Beta.E0.1","Beta.E0.2" )

  return(list(est.coefficients= as.matrix(est.coeff)))
  #return(list(est.coefficients= as.matrix(est.coeff), est.covariance.matrix= est.cov.mat, est.fix.parameters= as.matrix(est.fix.par)))

}


