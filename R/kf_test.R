############################ New Kalman filtering test script
####
library("pracma")
#'@title Simulation script for a nonlinear spring
#'@description Uses the 4th order Runge-Kutta method to solve the differential 
#'equation
#'\eqn{mx''+gx'+kx+\epsilon x^3=0}
#'corresponding to a nonlinear spring-mass-damper system with mass m. This is also known 
#'as the unforced Duffing oscillator. The model solves the differential equation subject to random forcing on the 
#'velocity and outputs the location and velocity of the mass as a function of time.
#'@param nsteps Scalar, the number of time steps to run the system for. Default 100
#'@param dt time step. Default 0.01
#'@param x0 array, vector, or similar. The initial position and velocity of the spring. Default c(1,0)
#'@param k Scalar stiffness parameter of the spring. Default 2.
#'@param m Scalar, mass of the object. Default 1
#'@param g Scalar damping parameter. Default 0.3 
#'@param epsilon Scalar nonlinearity term. Default 0(linear spring)
#'@param sigma_w Scalar, standard deviation of the random noise term added to the velocity. Default \eqn{10^-3}.
#'@returns The state vector (position and velocity) of the mass.
#'@export
sim_spring = function(nsteps = 100, dt=0.01,
                      x0=c(1,0), k=2,m=1,g=0.3,
                      epsilon=0, sigma_w = 0.001){
  state_vector = array(0,dim=c(2,nsteps))
  state_vector[,1] = x0
  ### Assuming that wth model matches the dynamics?
  h=dt
  a=d=1/6
  b=c=1/3
  k=k/m
  g=g/m
  epsilon=epsilon/m
  ## RK4
  for(ii in 2:nsteps){
    k1 = c(x0[2], -g*x0[2]-k*x0[1]-epsilon*x0[1]^3)
    x1 = x0+h*k1
    x1p = 0.5*(x0+x1)
    k2 =  c(x1p[2],-g*x1p[2]-k*x1p[1]-epsilon*x1p[1]^3 )
    x2 = x0+h*k2
    x2p = 0.5*(x0+x2)
    k3 = c(x2p[2],-g*x2p[2]-k*x2p[1]-epsilon*x2p[1]^3 )
    x3 = x0+h*k3
    k4 = c( x3[2], -g*x3[2]-k*x3[1]-epsilon*x3[1]^3 )
    state_vector[,ii] = state_vector[,ii-1] + h*(a*k1+b*k2+c*k3+d*k4) +c(0,sqrt(sigma_w))*rnorm(1)
    x0 = state_vector[,ii]
  }
  return(state_vector)
}
#'@title Kalman filter for the spring-mass damper system assuming perfect linearity.
#'@description Using the 4th order Runge-Kutta method on the ODE
#'\eqn{mx''+gx'+kx=0}, performs Kalman filtering of the position and velocity of the spring
#'mass damper system under observation error covariance matrix R. The method assumes that the spring is
#'linear, introducing model form error that prevents accurate model calibration. This is implemented as a  
#'demonstration of the package's calibration checking methodology.
#'@param truth 2\times nsteps array containing the mass position along the first row and velocity along the second. 
#'@param p0 2\times 2 matrix. The initial guess for the covariance matrix of the predictions. Default is \eqn{P0=[[0.01,0];[0,0.01]]}.
#'@param Q A covariance matrix for the perturbations applied to the spring-mass system. Default is the degenerate covariance matrix
#'[[0,0];[0,0.001^2]] corresponding to no perturbations to position and Gaussian perturbations to velocity with mean 0 and sstandard deviation 0.001.
#'@param R 2\times 2 observation error covariance. Default [[0.001,0];[0,0.001]]. This represents the observation errors through which the true state is 
#'observed. If you wish to simulate a system where velocity is unobserved, set the second diagonal entry to a large value, i.e. 1e6. 
#'@param k,m,g The spring stiffness, mass, and damping. Defaults are 2,1, 0.3 respectively.
#'@param dt The time step. Default is 0.01.
#'@returns a list:
#' $filt: The Kalman filtered state vector
#' $pred: The predicted state, used for the calibration estimation.
#' $predcov: the predictive covariance matrices.
#' $postcov: Posterior covariance matrices
#' $truth: The true state, same as the corresponding input parameter
#'@export
kalman_filter_spring = function( truth, P0=diag(c(0.01, 0.01)), 
                                 Q=diag(c(0,0.001^2)), 
                                 R=diag(c(1e-3,1e-3)), 
                                 k=2,m=1,g=0.3, dt=0.01 ){
#####
##Kalman filtering for the spring mass damper. Assumes linear, corrupted fwd model
## Note: We assume that we observe the state data[,i] + v, v\sim N(0,R)
I2 = diag(2)
data = truth
A = t(matrix(c( 0,1,-k/m,-g/m ), nrow=2,ncol=2))*dt#+diag(2)
FF_ = array(0,dim=c(2, 2, 5))
FF_[,,1] = diag(2)
for(ii in 2:5){
  FF_[,,ii] = A%*%FF_[,,ii-1]/(ii-1)
}
FF = rowSums(FF_[,,], dims=2) #Matrix form of RK4 for a linear system 
filtered = array(0,dim=dim(data))
filtered[,1] = data[,1]## Assuming we know this up to P0
predicted = filtered
observed = filtered
Pkk = array(0,dim=c(2,2,dim(data)[2]))
Pkk[,,1] = P0
Pkkn1 = array(0,dim=c(2,2,dim(data)[2]))
##
obs_err = t(MASS::mvrnorm(dim(data)[2],mu=c(0,0), Sigma = R))
for(ii in 2:dim(data)[2]){
  predicted[,ii] = FF%*%filtered[,ii-1]
  Pkkn1[,,ii] = FF%*%Pkk[,,ii-1]%*%t(FF) + Q
  observed[,ii] = data[,ii] + obs_err[,ii]
  innov = observed[,ii]-predicted[,ii]
  SS = Pkkn1[,,ii] + R
  Kgain = Pkkn1[,,ii]%*%solve(SS)
  filtered[,ii] = predicted[,ii] + Kgain%*%innov
  Pkk[,,ii] =
    (I2 - Kgain) %*% Pkkn1[,,ii] %*% t(I2 - Kgain) +
    Kgain %*% R %*% t(Kgain)
}
return(list(filt =filtered,pred=predicted,predcov=Pkkn1,postcov=Pkk, truth = data ))
}
#### Reliability plot building function
#'@title Reliability plot estimation
#'@description Estimates coverage rates and generates reliability plots for the 
#'probabilistic model based on inputs. Currently under development, so just returns the 
#'empirical coverages
#'@param epsvals The spring nonlinearity parameters
#'@param ground_truth The true state (position OR velocity)
#'@param preds A time series of predicted positions OR velocities
#'@param predstds The predictive standard deviations
#'@param ptest A vector of the nominal coverage levels.
#'@returns Empirical coverage estimates
#'@export 
rel_plot_spring = function(epsvals, ground_truth, preds, predstds, ptest ){
  ####
  #Want to biuld a data frame that covers everything.
  #Note: This is just the plot for one of the variables, position and
  #velocity should be done separately
  ####
  #Results should be 2 x ptest x epsvals
  cols <- hcl.colors(length(epsvals), palette = "viridis")
  ntrials=dim(preds)[2]
  # cvg_res = array(0, dim=c(3,length(ptest), length(epsvals)))
  # Prior: DIR(1,...,1)
  coverage = array(0,dim=c(length(epsvals), length(ptest)))
  dir_params = array(1,dim=c(length(epsvals), length(ptest) ))
  # colnames(coverage) <- as.character(nominal_levels)
  for (j in seq_along(ptest)) {
    alpha    <- 1 - ptest[j]
    z        <- qnorm(1 - alpha / 2)
    for(ee in seq_along(epsvals)){
      for(nn in 1:ntrials){
      true_q = ground_truth[ee,nn,]
      # pred_sd  <- sqrt(sapply(ekf_res$P_pred, function(P) P[1, 1]))
      lower    <- preds[ee,nn,] - z * predstds[ee,nn,]
      upper    <- preds[ee,nn,] + z * predstds[ee,nn,]
      coverage[ee, j] = coverage[ee, j]+sum( (true_q >= lower) & (true_q <= upper))
    }
    }
  }
  ###
  return(coverage)
}





