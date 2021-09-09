#' @title The Gibbs sampler for RFID with fixed K
#' @description MCMC main function for RFID with fixed number of subclones
#' @param Y observed count matrix
#' @param Init an Init object used for initializing the Markov chain
#' @param hyper the hyperparameters for the MCMC
#' @param B number of burn-in iterations
#' @param nmc number of post-burn-in iterations
#' @param K number of cell types
#' @param Print_Iter logical variable indicating whether or not to print the number of current iteration
#' @param lab_switch.check logical variable indicating whether label switching phenomenon is checked
#' @return mcmc an mcmc object for RFID

RFID_MCMC<-function(Y, Init, hyper, B, nmc, K, Print_Iter = FALSE, lab_switch.check = FALSE){
    
    purity=0.2
    n=apply(Y,2,sum)
    
    
    G = dim(Y)[1]
    S = dim(Y)[2]
    
  mcmc<-NULL
  # Latent counts are not necessary for the final inference;
  # only store the current values
  ## mcmc$ell=array(0,dim=c(G,S,K))
  mcmc$ell=array(0,dim=c(G,S,K))
  mcmc$ELL=mcmc$z=mcmc$x=mcmc$m=array(0,dim=c(G,S))
  mcmc$ell_tilde=array(0,dim=c(S,K))
  # Store the rest parameters for all iterations
  mcmc$Phi=array(NA,dim=c(G,K,nmc+B))
  mcmc$Theta=array(NA,dim=c(K,S,nmc+B))
  mcmc$p=mcmc$c=array(NA, dim=c(S,nmc+B))
  mcmc$gamma0 = mcmc$lambda = mcmc$c0=rep(NA,nmc+B)
   mcmc$r=matrix(NA, nrow = K, ncol = nmc+B)
 
  mcmc$alpha=array(NA,dim=c(G,nmc+B))
  # Initialization for the latent counts
  ##mcmc$x=Y-mcmc$z
  for (i in 1:G){
    for (j in 1:S){
      mcmc$ell[i,j,]=rmultinom(1,mcmc$ELL[i,j],rep(1/K,K))
    }
  }
  # Initialization for the rest parameters
  ## comment out by Jun Li Oct-15-2019
  mcmc$alpha[,1] = mcmc$lambda[1] = 1
  ## add By Jun Li Oct-15-2019
  ##mcmc$alpha[,1]=1/G
  ##mcmc$lambda[1] = 15
  
  mcmc$r[,1]=1
  mcmc$Phi[,,1]=Init$Phi
  mcmc$Theta[,,1]=Init$Theta
  mcmc$p[,1]=rep(1/2,S)
  mcmc$c[,1]=rep(1,S)
  mcmc$gamma0[1]=mcmc$c0[1]=1
  # Importing updating functions from cplusplus files
  # Rcpp::sourceCpp('update_pos_NBFA_RF_cpp.cpp')
  # Rcpp::sourceCpp('update_pos_cube_NBFA_RF_cpp.cpp')
  #MCMC sampling procedure
  #Time Computing Point 1
  # if (Print_Iter == TRUE){
  log_likelihood = rep(0, B + nmc)
  ptm<-proc.time()
    # }
  ##Initialize
  for (iter in 2:(B+nmc)){
    # Update PhiTheta
    print(iter)
    PhiTheta = matrix(mcmc$Phi[,,iter-1], ncol = K) %*% matrix(mcmc$Theta[,,iter-1], nrow = K)
    # PhiTheta = mcmc$Phi[,,iter-1]%*%mcmc$Theta[,,iter-1]
    # Sampling z
    ##comment out by Jun Li Oct-14_2019
    ## mcmc$z = update_z(mcmc$alpha[,iter-1], PhiTheta, Y)
   
      mcmc$z = update_z(mcmc$alpha[,iter-1], PhiTheta, Y,mcmc$lambda[iter-1])
   
   ## mcmc$z = update_z_purity(mcmc$alpha[,iter-1], PhiTheta, Y,mcmc$lambda[iter-1],purity)
     
    # Update x
     mcmc$x = Y - mcmc$z
    ## mcmc$x = Y
    # Sampling m
    mcmc$m = update_m(mcmc$alpha[,iter-1], mcmc$z, mcmc$lambda[iter-1])
    # Sampling lambda
    mcmc$lambda[iter] = rgamma(1, shape = hyper$u0 + sum(mcmc$m),
                             rate = hyper$v0 - sum(log(1-mcmc$p[, iter - 1])))
                               
    
    # Sampling alpha
    mcmc$alpha[,iter] = update_alpha(hyper$g0, hyper$h0, mcmc$m, mcmc$p[,iter-1])
    
    # Sampling ELL
    mcmc$ELL = update_ELL(PhiTheta, mcmc$x)
    # Sampling ell
    mcmc$ell = update_ell(matrix(mcmc$Phi[,,iter-1], ncol = K),matrix(mcmc$Theta[,,iter-1], nrow = K),mcmc$ELL)
    # Computing sum_j ell_ijk
    if (K == 1){
      ell_ik = matrix(apply(mcmc$ell,c(1),sum),ncol=1)
    }else{
      ell_ik = apply(mcmc$ell,c(1,3),sum)
    }
    # Sampling Phi
    
    ##tmp=update_Phi(hyper$eta,ell_ik,matrix(mcmc$Theta[,,iter-1], nrow = K),mcmc$p[,iter-1])
     ## tmp=tmp*0.5+Init$Seeded_Phi*0.5
     ##   p=rbeta(1,1,1)
     ##for(i in 1:length(tmp[1,]))
     {
        
        ##  indicator=rbinom(1,1,p)
        ## if(indicator==1)
         {
             ## mcmc$Phi[,i,iter]=tmp[,i]
             ## }else{
             
             ##  mcmc$Phi[,i,iter]=Init$Seeded_Phi[,i]
         }
         
         
     }
      mcmc$Phi[,,iter]=update_Phi(hyper$eta,ell_ik,matrix(mcmc$Theta[,,iter-1], nrow = K),mcmc$p[,iter-1])
     ##   mcmc$Phi[,,iter]= Init$Seeded_Phi[,2]
    ## Fixed phi, comment out by Jun Li 09-06-2019
    ## tmp<-update_Phi(hyper$eta,ell_ik,matrix(mcmc$Theta[,,iter-1], nrow = K),mcmc$p[,iter-1])
    ##tmp1<-tmp[,K]
    ##tmp2<-mcmc$Phi[,1:(K-1),iter-1]
    ##total<-cbind(tmp2,tmp1)
    ## mcmc$Phi[,,iter]= mcmc$Phi[,,iter-1]
     
     
    # Computing sum_i ell_ijk
    if (K == 1){
      ell_jk = matrix(apply(mcmc$ell,c(2),sum), ncol = 1)
    }else{
      ell_jk = apply(mcmc$ell,c(2,3),sum)
    }
    # Sampling Theta
    mcmc$Theta[,,iter]=update_Theta(mcmc$r[,iter-1],mcmc$c[,iter-1],mcmc$p[,iter-1],ell_jk,matrix(mcmc$Phi[,,iter], ncol = K))
    # Sampling p
    ## mcmc$p[,iter]=update_p(hyper$a0,hyper$b0,n,matrix(mcmc$Theta[,,iter], nrow = K),mcmc$alpha[,iter])
    mcmc$p[,iter]=update_p(hyper$a0,hyper$b0,n,matrix(mcmc$Theta[,,iter], nrow = K))
    
    # Sampling c
    mcmc$c[,iter]=update_c(hyper$e0,hyper$f0,mcmc$r[,iter-1],matrix(mcmc$Theta[,,iter], nrow = K))
    #Sampling ell_tilde
    mcmc$ell_tilde=update_ell_tilde(ell_jk,mcmc$r[,iter-1])
    # Computing sum_j ell_tilde_jk
    if (K == 1){
      ell_tilde_k = sum(mcmc$ell_tilde)
    }else{
      ell_tilde_k = apply(mcmc$ell_tilde,c(2),sum)
    }
    # Computing K_plus
    # K_plus = sum(ell_tilde_k > 0)
    # Computing p_tilde_j and p_tilde_tilde
    p_tilde_j = -log(1-mcmc$p[,iter])/(mcmc$c[,iter]-log(1-mcmc$p[,iter]))
    p_tilde = -sum(log(1-p_tilde_j))/(mcmc$c0[iter-1]-sum(log(1-p_tilde_j)))
    print("p_tilde_j")
    print(p_tilde_j)
    print("ell_tilde_k")
    print(ell_tilde_k)
    # Sampling gamma0
    mcmc$gamma0[iter]=update_gamma0(hyper$a0,hyper$b0,p_tilde,mcmc$gamma0[iter-1],ell_tilde_k)
    print("Gamma0")
    print(mcmc$gamma0[iter])
    # Sampling c0
    mcmc$c0[iter]=update_c0(hyper$e0,hyper$f0,mcmc$gamma0[iter],mcmc$r[,iter-1])
    # Sampling r
    # mcmc$r[,iter]=update_r(mcmc$c0[iter-1],mcmc$gamma0[iter-1],mcmc$ell_tilde[,,iter],mcmc$p[,iter],mcmc$c[,iter],ell_k)
    # mcmc$r[,iter]=c(6.4,8.78,4.7778)
    print("c0")
    print (mcmc$c0[iter])
   

    mcmc$r[,iter]=update_r(mcmc$c0[iter], mcmc$gamma0[iter], p_tilde_j, ell_tilde_k)
    print("after")
    print("crt ell_jk")
    print(ell_tilde_k)
    print("ell_jk")
    print(ell_jk)
    print("Theta")
    print(mcmc$Theta[,,iter])
    print("r")
    print(mcmc$r[,iter])
    print("P")
    print(mcmc$p[,iter])
    if ((Print_Iter == TRUE) && (floor(iter/100) == iter/100)){
      print(paste("Iteration: ",iter))
    }
    PhiTheta = matrix(mcmc$Phi[,,iter], ncol = K) %*%
      matrix(mcmc$Theta[,,iter], nrow = K)
    log_likelihood[iter] =log_likelihood_Y(Y, PhiTheta, alpha = mcmc$alpha[, iter], p = mcmc$p[, iter], lambda = mcmc$lambda[iter])
    
    
  }
  #Time computing Point 2
  runtime = proc.time() - ptm
  colnames(runtime) = NULL
  if (Print_Iter == TRUE){
    print(paste("Total runtime: ", ceiling(runtime[3]), "s", sep = ""))
  }
  # Normalize Theta such that sum_k theta_kj = 1 for all j
  sample_ind = (B + 1):(B + nmc)
  iter_ref = which.max(log_likelihood[sample_ind])
  perm_all = combinat::permn(1:K)
  mcmc$Theta_norm=array(NA,dim=c(K, S, (B + nmc)))
  for (iter in 1:(B + nmc)){
    if (lab_switch.check == TRUE){
      mse = rep(NA, factorial(K))
      for (t in 1:factorial(K)){
        perm_tmp = perm_all[[t]]
        mse[t] = sum((mcmc$Theta[perm_tmp, , iter] - mcmc$Theta[, , iter_ref])^2)
      }
      ind_perm = perm_all[[which.min(mse)]]
      mcmc$Phi[, , iter] = mcmc$Phi[, ind_perm, iter]
      mcmc$Theta[, , iter] = mcmc$Theta[ind_perm, , iter]
    }
    mcmc$Theta_norm[,,iter]=mcmc$Theta[,,iter]/
      matrix(1,nrow=K,ncol=1)%*%matrix(apply(matrix(mcmc$Theta[,,iter], nrow = K),MARGIN = c(2),sum),nrow=1,ncol=S)
  }
  # Label switching detection
  return(mcmc)
}
