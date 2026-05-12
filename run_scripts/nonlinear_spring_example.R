#### Now we try to do some analysis
## Let's use those functions that we have. Should clean out probably the first 100
## observations at least to keep the filter converged
## Then we can go from there
if(TRUE){
  ptest = c(0.05,0.25,0.5,0.75,0.95)
  ##
  k = 4; g = 0.3; m = 1
  epsvals = seq(0,2,by=0.2)
  ntrials = 100
  dt = 0.01
  nsteps = 500
  system_noise_var = 0.0001 #Q
  obs_err_cov = matrix(c(0.01,0,0,10000), nrow=2)
  ## Discard first 200 observations
  nburn = 200
  ## Here we will need to generate matrices for the coverages and such
  preds = array(0,dim=c(length(epsvals), ntrials, nsteps))
  predstds = array(0,dim=c(length(epsvals), ntrials,nsteps))
  ground_truth = array(0,dim=c(length(epsvals), ntrials,2,nsteps))
  ## First: Let's set up the loop
  for(ee in 1:length(epsvals)){
    
    
    for(ii in 1:ntrials){
      ground_truth[ee,ii,,] = sim_spring(
        nsteps=nsteps,dt=dt,x0=c(1,0),
        k=k,m=m,g=g,epsilon = epsvals[ee], sigma_w = sqrt(system_noise_var)
      )
      kalman_est = kalman_filter_spring(ground_truth[ee,ii,,],
                                 Q=diag(c(0.000,system_noise_var)),
                                 k=k,g=g,m=m,R=obs_err_cov,dt=dt)
      preds[ee,ii,] = kalman_est$pred[1,]
      predstds[ee,ii,] = sqrt(kalman_est$predcov[1,1,])
      if(ii==1){
        plot((1:nsteps)*dt-dt, ground_truth[ee,ii,1,], col="black", type="l", lty=2,
             xlab="Time", ylab="Position",main=paste("Kalman filtered and true states, eps=",epsvals[ee]))
        lines((1:nsteps)*dt-dt,kalman_est$pred[1,],col="red")
        lines((1:nsteps)*dt-dt,kalman_est$pred[1,]+1.96*predstds[ee,ii,],col="red", lty=2)
        lines((1:nsteps)*dt-dt,kalman_est$pred[1,]-1.96*predstds[ee,ii,],col="red", lty=2)
      }
    }
  }
}
## Now proceed

# reliability_plot(epsvals, ground_truth[,,1,], preds, predstds, ptest)
eval_ndcs = seq(nburn,500,by=20)
cvgs = rel_plot_spring(epsvals, ground_truth[,,1,eval_ndcs],
                        preds[,,eval_ndcs],
                        predstds[,,eval_ndcs], ptest)


