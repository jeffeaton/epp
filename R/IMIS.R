## Modified from IMIS package by Le Bao (http://cran.r-project.org/web/packages/IMIS/)

IMIS <- function(B0, B, B.re, number_k){
  ## B0 = B*100
  ## X_all = X_k = sample.prior(B0)                                # Draw initial samples from the prior distribution
  X_k <- sample.prior(B0)                                
  X_all = matrix(0, B0 + B*number_k, dim(X_k)[2])
  n_all <- 0
  X_all[1:B0,] <- X_k
  
  ## Sig2_global = cov(X_all)        # the prior covariance
  Sig2_global = cov(X_all[1:B0,])        # the prior covariance
  stat_all = matrix(NA, 6, number_k)                            # 6 diagnostic statistics at each iteration
  ##center_all = prior_all = like_all = NULL                      # centers of Gaussian components, prior densities, and likelihoods
  center_all = NULL
  prior_all = like_all = gaussian_sum = vector("numeric", B0 + B*number_k)
  sigma_all = list()                                            # covariance matrices of Gaussian components
  
  iter.start.time = proc.time()
  for (k in 1:number_k ){
    ##Rprof(paste("IMIS-k", k, ".out", sep=""))
    ptm.like = proc.time()
    ## prior_all = c(prior_all, prior(X_k))                # Calculate the prior densities
    ## like_all = c(like_all, likelihood(X_k))             # Calculate the likelihoods
    prior_all[n_all + 1:dim(X_k)[1]] <-  prior(X_k)
    like_all[n_all + 1:dim(X_k)[1]] <-  likelihood(X_k)
    ptm.use = (proc.time() - ptm.like)[3]
    if (k==1)   print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,2), "minutes"))
    which_pos <- which(like_all[1:(n_all + dim(X_k)[1])] > 0)
    
    ## if (k==2){
    ##   gaussian_pos = matrix(dmvnorm(X_all[which_pos,], center_all[1,], sigma_all[[1]]), 1, length(which_pos))
    ## }
    ## if (k>2){
    ##   gaussian_new = matrix(NA, k-1, length(which_pos))
    ##   gaussian_new[1:(k-2), 1:dim(gaussian_pos)[2]] = gaussian_pos
    ##   gaussian_new[k-1, ] = dmvnorm(X_all[which_pos,], center_all[k-1,], sigma_all[[k-1]])
    ##   for (j in 1:(k-2))    gaussian_new[j, (dim(gaussian_pos)[2]+1):length(which_pos) ] = dmvnorm(X_all[which_pos[(dim(gaussian_pos)[2]+1):length(which_pos)],], center_all[j,], sigma_all[[j]])
    ##   gaussian_pos = gaussian_new
    ## }

    
    ## unsure if i should really multithread this...
    if(k > 2)
      gaussian_sum[(n_pos+1):length(which_pos)] <- rowSums(sapply(1:(k-2), function(j)dmvnorm(X_all[which_pos[(n_pos+1):length(which_pos)],], center_all[j,], sigma_all[[j]])))
    if(k > 1){
      n_pos <- length(which_pos)
      gaussian_sum[1:n_pos] <- gaussian_sum[1:n_pos] + dmvnorm(X_all[which_pos,], center_all[k-1,], sigma_all[[k-1]])
    }

    
    ##if (k==1)   envelop_all = prior_all                 # envelop stores the sampling densities
    if (k==1)   envelop_pos = prior_all[which_pos]        # envelop stores the sampling densities
    ##if (k>1)    envelop_all = apply( rbind(prior_all*B0/B, gaussian_all), 2, sum) / (B0/B+D+(k-2))
    if (k>1)    envelop_pos = (prior_all[which_pos]*B0/B + gaussian_sum[1:n_pos]) / (B0/B+(k-1))
    ##Weights = prior_all*like_all / envelop_all  # importance weight is determined by the posterior density divided by the sampling density
    Weights = prior_all[which_pos]*like_all[which_pos]/ envelop_pos  # importance weight is determined by the posterior density divided by the sampling density
    ##stat_all[1,k] = log(mean(Weights))                  # the raw marginal likelihood
    stat_all[1,k] = log(mean(Weights)*length(which_pos)/(n_all+dim(X_k)[1]))                  # the raw marginal likelihood
    Weights = Weights / sum(Weights)
    stat_all[2,k] = sum(1-(1-Weights)^B.re)             # the expected number of unique points
    stat_all[3,k] = max(Weights)                                # the maximum weight
    stat_all[4,k] = 1/sum(Weights^2)                    # the effictive sample size
    stat_all[5,k] = -sum(Weights*log(Weights), na.rm = TRUE) / log(length(Weights))     # the entropy relative to uniform
    stat_all[6,k] = var(Weights/mean(Weights))  # the variance of scaled weights
    if (k==1)   print("Stage   MargLike   UniquePoint   MaxWeight   ESS   IterTime")
    iter.stop.time = proc.time()
    print(c(k, round(stat_all[1:4,k], 3), as.numeric(iter.stop.time - iter.start.time)[3]))
    iter.start.time = iter.stop.time

    ## choose the important point
    important = which(Weights == max(Weights))
    if (length(important)>1)  important = important[1]
    ##if (is.matrix(X_all))     X_imp = X_all[important,]                               # X_imp is the maximum weight input
    X_imp = X_all[which_pos[important],]
    center_all = rbind(center_all, X_imp)
    distance_all = mahalanobis(X_all[1:(n_all+dim(X_k)[1]),], X_imp, diag(diag(Sig2_global)) )
    ##label_nr = sort(distance_all, decreasing = FALSE, index=TRUE)             # Sort the distances
    label_nr = sort(distance_all, decreasing = FALSE, index=TRUE, method="quick")             # Sort the distances
    which_var = label_nr$ix[1:B]                                                              # Pick B inputs for covariance calculation

    ###########
    weight_close <- Weights[match(which_var, which_pos)]
    weight_close[!which_var %in% which_pos] <- 0
    
    ##if (is.matrix(X_all))     Sig2 = cov.wt(X_all[which_var,], wt = Weights[which_var]+1/length(Weights), cor = FALSE, center = X_imp, method = "unbias")$cov
    Sig2 = cov.wt(X_all[which_var,], wt = weight_close+1/(n_all + dim(X_k)[1]), cor = FALSE, center = X_imp, method = "unbias")$cov
    ##sigma_all[[D+k-1]] = Sig2
    sigma_all[[k]] = Sig2
    n_all <- n_all + dim(X_k)[1]
    if (is.matrix(X_all))     X_k = rmvnorm(B, X_imp, Sig2)                           # Draw new samples
    X_all[n_all + 1:B,] <- X_k
    ##if (is.matrix(X_all))     X_all = rbind(X_all, X_k)

    ##Rprof(NULL)
    ##save(list=ls(), file="IMIS-intermediate.RData")
    
    if (stat_all[2,k] > (1-exp(-1))*B.re)       break
  } # end of k

  ##nonzero = which(Weights>0)
  ##which_X = sample(nonzero, B.re, replace = TRUE, prob = Weights[nonzero])
  which_X = sample(which_pos, B.re, replace = TRUE, prob = Weights)
  if (is.matrix(X_all)) resample_X = X_all[which_X,]

  return(list(stat=t(stat_all), resample=resample_X, center=center_all))
} # end of IMIS
