## ------------------------------------------------------
## Functions for Bayesian online changepoint detection
## ------------------------------------------------------

## Classic BOCPD from Prof Fearnhead's 2007 paper

## Linear trend model for all segments
## Closed form marginal likelihood


## Trend model, Normal prior
integrate.fun5 <- function(hmm, cpt, cps, params.list, prior.list) {
  alpha <- params.list[[1]]
  sigma <- params.list[[2]]
  
  a.prior <- prior.list[[1]]
  sig.prior <- prior.list[[2]]
  
  n <- cpt - cps
  y <- matrix(hmm$x[(cps+1):cpt], nrow=n, ncol=1)
  x <- c(rep(1, n), 1:n)
  x <- matrix(x, nrow=n, ncol=2, byrow=FALSE)
  
  mu0 <- matrix(a.prior[1:2], nrow=2, ncol=1)
  D0 <- diag(a.prior[3:4])
  # V0 <- sigma * D0
  D0.inv <- diag(1/a.prior[3:4])
  
  D1 <- solve(t(x) %*% x + D0.inv)
  # V1 <- sigma * D1
  # mu1 <- D1 %*% (t(x) %*% y + D0.inv %*% mu0)
  
  P1 <- diag(n) - x %*% D1 %*% t(x)
  const.y <- t(y) %*% P1 %*% y 
  
  u0 <- sig.prior[1]           
  v0 <- sig.prior[2]
  u1 <- u0 + n/2
  v1 <- v0 + as.numeric(const.y/2)
  
  det.D0 <- sum(log(a.prior[3:4]))    # 0.5*(log(det(D1)) - log(det(D0)))
  det.D1 <- log(D1[1, 1]*D1[2, 2] - D1[1, 2]*D1[2, 1])
  Int.par <- -0.5*n*log(2*pi) + 0.5*(det.D1 - det.D0) - lgamma(u0) + u0*log(v0) + 
    lgamma(u1) - u1*log(v1) 
  return(as.numeric(Int.par))
}


## Trend model, truncated normal prior
integrate.fun5.v2 <- function(hmm, cpt, cps, params.list, prior.list) {
  alpha <- params.list[[1]]
  sigma <- params.list[[2]]
  
  a.prior <- prior.list[[1]]
  sig.prior <- prior.list[[2]]
  
  n <- cpt - cps
  y <- matrix(hmm$x[(cps+1):cpt], nrow=n, ncol=1)
  x <- c(rep(1, n), 1:n)
  x <- matrix(x, nrow=n, ncol=2, byrow=FALSE)
  
  mu0 <- matrix(a.prior[1:2], nrow=2, ncol=1)
  D0 <- diag(a.prior[3:4])
  # V0 <- sigma * D0
  D0.inv <- diag(1/a.prior[3:4])
  
  D1 <- solve(t(x) %*% x + D0.inv)
  # V1 <- sigma * D1
  mu1 <- D1 %*% (t(x) %*% y + D0.inv %*% mu0)
  
  P1 <- diag(n) - x %*% D1 %*% t(x)
  const.y <- t(y) %*% P1 %*% y 
  
  # phi0 = TruncatedNormal::pmvnorm(lb=c(0, -1), ub=c(1, 0), mu=mu0, sigma=sigma*D0)[1]
  # phi1 = TruncatedNormal::pmvnorm(lb=c(0, -1), ub=c(1, 0), mu=mu1, sigma=sigma*D1)[1]
  phi0 = mvtnorm::pmvnorm(lower=c(0, -1), upper=c(1, 0), mean=as.vector(mu0),
                          sigma=D0*sigma)[1]
  phi1 = mvtnorm::pmvnorm(lower=c(0, -1), upper=c(1, 0), mean=as.vector(mu1),
                          sigma=D1*sigma)[1]
  # for truncated normal
  # the posterior may be out of range
  if ((phi1 < 1e-14) | (is.nan(phi1))) { 
    log.phi = -100 - log(phi0)
  } else { log.phi = log(phi1) - log(phi0) }
  
  u0 <- sig.prior[1]           
  v0 <- sig.prior[2]
  u1 <- u0 + n/2
  v1 <- v0 + as.numeric(const.y/2)
  
  det.D0 <- sum(log(a.prior[3:4]))
  det.D1 <- log(D1[1, 1]*D1[2, 2] - D1[1, 2]*D1[2, 1])
  Int.par <- -0.5*n*log(2*pi) + 0.5*(det.D1 - det.D0) - lgamma(u0) + u0*log(v0) + 
    lgamma(u1) - u1*log(v1) + log.phi
  return(as.numeric(Int.par))
}


## Run-length functions
# modification when using min seg length
runl.fun1 <- function(l, ptr, msl=1) {
  # if msl is specified, use some form of a shifted distribution
  gx <- ((1-ptr)^(l-msl)) * ptr
  Gx <- 1 - (1-ptr)^(l-msl)
  # this includes the case of msl = 1
  return(list(gx=gx, Gx=Gx))
}

runl.fun2 <- function(l, ptr, msl=1) {
  # before reaching msl, g(x) = 0, we need a new normalised distribution g*(x)
  # normalising constant: 1 - sum_x=1:(msl-1) g(x)
  cdf.msl <- 1 - sum( (1-ptr)^(1:(msl-1)-1) * ptr ) 
  gx <- ((1-ptr)^(l-1)) * ptr / cdf.msl
  Gx <- (1 - (1-ptr)^(l-msl)) * (1-ptr)^(msl-1) / cdf.msl
  return(list(gx=gx, Gx=Gx))
}


## log-sum-exp for the log-likelihood
logsumexp <- function(x){
  x <- as.matrix(x)
  a <- apply(x, 1, max)
  log(rowSums(exp(x-a))) + a
}


## Forward filtering recursion of BOCPD
forward.cpt.int <- function(hmm, cpt, nmodel, params.ini.list, prior.list, 
                            intfun.list, msl) {
  # PF for whether particle filter is used 
  
  l <- length(hmm$forward[[cpt-1]]) + 1
  pos <- hmm$ParsID[[cpt-1]] 
  
  # a matrix for storing the integrated likelihood
  intPar.mat <- matrix(0, nrow=nmodel, ncol=l)
  
  ## stay: j = cpt-msl-1, ..., msl, 0
  for (i in 1:(l-1)) {   
    # used to be j in ((cpt-2):0)[pos]
    j <- ((cpt-2):0)[pos[i]]
    params.list <- params.ini.list
    # The marginal log-likelihood 
    intPar.m <- integrate.par(hmm, cpt=cpt, cps=j, nmodel, params.list,
                              prior.list, intfun.list)
    # the order is t-2, ..., 0 (closest changepoint first)
    intPar.mat[, i+1] <- unlist(intPar.m)
  }
  
  ## changepoint: j = cpt-msl 
  params.list <- params.ini.list
  intPar.m <- integrate.par(hmm, cpt=cpt, cps=cpt-msl, nmodel, params.list,
                            prior.list, intfun.list)
  intPar.mat[, 1] <- unlist(intPar.m)
  
  # the extra one
  intExtra.m <- integrate.par(hmm, cpt=cpt-1, cps=cpt-msl, nmodel, params.list,
                              prior.list, intfun.list)
  intExtra <- unlist(intExtra.m)
  
  # putting everything together
  intPar.list <- list(intPar.mat)
  intExtra.list <- list(intExtra)
  
  return(list(intPar.list=intPar.list, intExtra.list=intExtra.list, 
              hmm=hmm))
}


## A function to calculate the integrated likelihood
integrate.par <- function(hmm, cpt, cps, nmodel, params.list, prior.list, 
                          intfun.list) {
  # rfun.list calculates the function to be integrated, is NULL if no mumerical
  # integration is required
  if ((length(intfun.list) == nmodel) & (length(params.list) == nmodel) &
      (length(prior.list) == nmodel)) {
    # compute the integral
    Int <- vector(mode="list", length=nmodel)
    for (i in 1:nmodel) {
      params.list.i <- params.list[[i]]
      prior.list.i <- prior.list[[i]]
      intfun.i <- intfun.list[[i]]
      Int[[i]] <- intfun.i(hmm, cpt, cps, params.list=params.list.i, 
                           prior.list=prior.list.i)
    }
    return(Int)
    
  } else {print("Incorrect list")}
}


## Calculate the forward quantity, using re-sampling
forward.cpt.calc.resample <- function(hmm, cpt, nmodel, maxp, np, 
                                      prune=FALSE, delay=0, ptr, msl) {
  # maxp is the maximum number of changepoints to keep
  # np is the number of changepoints being pruned to
  # ptr is the transition parameter in the run-length distribution
  # msl is minimum segment length
  
  # calculate the weights W_cpt^{s}
  l <- length(hmm$forward[[cpt-1]]) + 1
  pos <- hmm$ParsID[[cpt-1]]  
  
  # last cpt < i-msl
  lparticles12 <- numeric(0)
  lsum.temp <- matrix(NA, nrow=l, ncol=nmodel)
  for (s in 2:l) {
    # compute the marginalised integration I_t^{s, m} for each s
    # intPar.lsum <- hmm$intPar[[cpt]][[1]][, s]
    # intPar.psum <- hmm$intPar[[cpt-1]][[1]][, s-1]
    intPar.lsum <- sapply(1:nmodel, FUN=function(m) {
      hmm$intPar[[cpt]][[m]][, s] })
    intPar.psum <- sapply(1:nmodel, FUN=function(m) {
      hmm$intPar[[cpt-1]][[m]][, s-1] })  
    
    # compute the weights as in the 2007 paper
    lsum <- intPar.lsum + log(hmm$pModels)
    psum <- intPar.psum + log(hmm$pModels)
    division <- logsumexp(matrix(lsum, nrow=1)) - logsumexp(matrix(psum, nrow=1))
    lparticles12 <- c(lparticles12, division)
    # for used in the MAP calculation for pruning
    lsum.temp[s, ] <- lsum 
  }
  f12 <- lparticles12 + log(hmm$transition[1:(cpt-1), 1][pos]) +
    log(hmm$forward[[cpt-1]])
  
  # last cpt = i-msl
  # intPar.lsum <- hmm$intPar[[cpt]][[1]][, 1]
  # intPar.psum <- hmm$intExtra[[cpt]][[1]]
  intPar.lsum <- sapply(1:nmodel, FUN=function(m) {
    hmm$intPar[[cpt]][[m]][, 1] })
  intPar.psum <- sapply(1:nmodel, FUN=function(m) {
    hmm$intExtra[[cpt]][[m]] })  
  
  lsum <- intPar.lsum + log(hmm$pModels)
  psum <- intPar.psum + log(hmm$pModels)
  lparticles34 <- logsumexp(matrix(lsum, nrow=1)) - logsumexp(matrix(psum, nrow=1))
  logsum.f34 <- 0
  lsum.temp[1, ] <- lsum    # for the MAP calculation
  f34 <- lparticles34 + log(hmm$transition[msl, 2]) + logsum.f34
  
  # the order here is t-msl, t-msl-1, .., msl, 0
  temp.forward <- c(f34, f12)
  temp.pos <- c(msl, pos+1)   # push the index 1 step forward
  temp.norm <- logsumexp(matrix(temp.forward, nrow=1))
  temp.weight <- exp(temp.forward - temp.norm)
  
  # Pruning where there is more than maxp candidate cpts
  if ((prune == TRUE) & (l > maxp))  {
    # re-sampling
    if (delay > 0) {
      swt <- sum(temp.weight[-(1:delay)])   
      SR <- sor.cpt(cpts.weight=temp.weight[-(1:delay)], np=np, thresh=NULL)
      keepID <- c(1:delay, SR$keepID+delay)
      keepWt <- c(temp.weight[1:delay], SR$keepWt * swt)
      
    } else if (delay == 0) {
      SR <- sor.cpt(cpts.weight=temp.weight, np=np, thresh=NULL)
      keepID <- SR$keepID
      keepWt <- SR$keepWt
    }
    
    hmm$ParsID[[cpt]] <- temp.pos[keepID]
    hmm$forward[[cpt]] <- temp.weight[keepID]
    hmm$forward[[cpt]] <- hmm$forward[[cpt]] / sum(hmm$forward[[cpt]])
    
    # update intPar here, so we don't need to trace keepID later
    intPar.prune <- list(matrix(hmm$intPar[[cpt]][[1]][, keepID], ncol=length(keepID)))
    hmm$intPar[[cpt]] <- intPar.prune
    
  } else {
    # if no re-sampling is required
    hmm$ParsID[[cpt]] <- temp.pos
    hmm$forward[[cpt]] <- temp.weight
  }
  
  return(hmm)
}


## Stratified re-sampling
sor.cpt <- function(cpts.weight, np, thresh=NULL) {
  # cpts.weight: the forward quantity hmm$forwatd[[t]] (t > maxp)
  # np: number of particles to keep (np is n in the code)
  # thresh: threshold (alpha in the paper)
  
  # first re-scale to have the sum = 1
  cpts.weight = cpts.weight / sum(cpts.weight)
  
  ## Find the threshold
  if (is.null(thresh)) {
    m <- length(cpts.weight)
    l.pre <- l <- 0   
    # l counts the particles with small weights, these are the ones to be re-sampled
    # if l* is the number of particles to keep with probability 1, 
    # then we need to sample (n - l*) particles from l particles
    thresh <- 1 / np     
    # to initialise, set thresh = 1/np, note that c.value = np = 1/thresh
    sum.weight <- 0
    # condition: sum.weight/thresh + (m-l) = np
    for (i in 1:m) {
      if (cpts.weight[i] < thresh) {
        sum.weight <- sum.weight + cpts.weight[i]
        l <- l + 1
      }
    }
    
    while ((l != l.pre) & (l != (m-np))) {
      # in the c++ code, the 2nd condition is (l != (np-m))
      # as np-m <= 0 and l >= 0, dose this just mean l != 0 ?
      # should this be l != (m-np) ?  so that np + l - m won't be 0
      
      # update the threshold
      new.n <- (np + l - m) / sum.weight  # np - (m - l) particles to resample
      thresh <- 1 / new.n     
      # sometimes the threshold is 0, why?
      
      l.pre <- l
      l <- 0
      sum.weight <- 0
      # condition: sum.weight/thresh + (m-l) = np
      for (i in 1:m) {
        if (cpts.weight[i] < thresh) {
          sum.weight <- sum.weight + cpts.weight[i]
          l <- l + 1
        }
      }
    }
  }
  
  ## re-sample the particles 
  ## As in the HMM particle filter paper (2003)
  Ln <- sum(cpts.weight >= thresh)   # the ones to keep with original weights
  weights.Ln = cpts.weight[cpts.weight >= thresh]
  
  Sn <- np - Ln                      # the number to sample from the rest
  if (Sn > 0) {
    weights <- cpts.weight[cpts.weight < thresh]
    K <- sum(weights) / Sn    # the 2003 paper uses u <- u + K
    # K <- thresh    # the 2007 paper uses u <- u + thresh
    u <- runif(1, min=0, max=K)
    i <- 1
    keep <- numeric(0)
    while (i <= length(weights)) {
      u <- u - weights[i]
      if (u < 0) {
        keep <- c(keep, i)
        u <- u + K
        # the original c++ code uses u <- u + 1
      }
      i <- i + 1
    }
    weights.Sn = rep(thresh, times=length(keep))       # assign equal weights
    
    keepID <- c(which(cpts.weight >= thresh), which(cpts.weight < thresh)[keep])
    keepWt <- c(weights.Ln, weights.Sn)
    
  } else {
    keepID <- which(cpts.weight >= thresh)
    keepWt <- weights.Ln
  }
  
  ## the final weights
  oID <- order(keepID)
  keepID <- keepID[oID]
  keepWt <- keepWt[oID]
  
  return(list(keepID=keepID, keepWt=keepWt))
}


## Backward simulation for changepoints
simulate.cpt <- function(hmm, nmodel, msl=2) { 
  # sample changepoints and model index
  nt <- hmm$Nx
  temp <- hmm$forward[[nt]]
  
  tau.new <- ((nt-1):0)[hmm$ParsID[[nt]]]
  t.new <- sample(1:length(temp), size=1, prob=temp)
  intPar.sum <- sapply(1:nmodel, FUN=function(m) {
    hmm$intPar[[nt]][[m]][, t.new] })
  pmod.post <- exp(intPar.sum)
  m.new <- sample(1:nmodel, size=1, prob=pmod.post)
  t.list <- tau.new[t.new]
  m.list <- m.new
  
  # change the following
  while (tau.new[t.new] >= 2*msl) {
    nt <- tau.new[t.new]
    temp <- hmm$forward[[nt]]
    
    # need to multiply hmm$forward by the transition prob P(C_nt+1 = nt | C_nt = s)
    # for geometric distribution, this is a constant, so we can omit that
    # otherwise, read it from the transition probability table
    
    tau.new <- ((nt-1):0)[hmm$ParsID[[nt]]]
    t.new <- sample(1:length(temp), size=1, prob=temp)
    intPar.sum <-  sapply(1:nmodel, FUN=function(m) {
      hmm$intPar[[nt]][[m]][, t.new] })
    pmod.post <- exp(intPar.sum)
    m.new <- sample(1:nmodel, size=1, prob=pmod.post)
    t.list <- c(tau.new[t.new], t.list)
    m.list <- c(m.new, m.list)
  }
  
  # the following is for the first changepoint at 0 at t=1
  # didn't manage to put it in the while loop...
  if (tau.new[t.new] > 0) {   
    t.list <- c(0, t.list)
    m.list <- c(1, m.list)
  }
  
  return(list(breaks=t.list, models=m.list))
}


## MAP estimation of changepoints
map.cpt <- function(hmm, nmodel, global.params, runl.fun, msl=2) { 
  # last changepoints and model index
  tstar <- rep(NA, hmm$Nx)
  mstar <- rep(NA, hmm$Nx)
  
  # model probability
  pm <- global.params[[1]]
  # transition parameters
  ptr <- global.params[[2]]
  
  logMAP.list <- vector(mode="list", length=hmm$Nx)
  logPosterior.list <- vector(mode="list", length=hmm$Nx)   # for diagnosing
  # P_t^MAP = Prob{changepoint at t | history, y_1:tstar},
  # where tstar >= t+msl and t=0, msl, msl+1, ..., hmm$Nx-msl
  # p_t(s, m) = Prob(Ct=s, Mt=m, history, y_1:t) = (1-G(t-s-1)) * I_t^{s, m} * pm * P_s^MAP
  # P_t^MAP = max_{s, m} { p_t(s, m) * g(t-s) / (1-G(t-s-1)) }
  
  logMAP.list[[1]] <- 0
  
  ## from t=msl to 2msl-1, tau=0
  for (i in msl:(2*msl-1)) {
    intPar.sum <-  sapply(1:nmodel, FUN=function(m) {
      hmm$intPar[[i]][[m]][, 1] })
    logPosterior.list[[i]] <- intPar.sum + log(hmm$pModels) + logMAP.list[[1]]
    logMAP.temp <- logPosterior.list[[i]] + log(runl.fun(i, ptr, msl)$gx)   
    logMAP.list[[i+1]] <- logMAP.temp[which.max(logMAP.temp)]
    tstar[i] <- 0
    mstar[i] <- which.max(logMAP.temp)
  }
  
  ## t=2*msl, ..., hmm$Nx
  for (i in (2*msl):hmm$Nx) { 
    l <- length(hmm$ParsID[[i]])
    pos <- ((i-1):0)[hmm$ParsID[[i]]]
    logMAP.temp <- logPosterior.temp <- matrix(0, nrow=l, ncol=nmodel)
    
    # for changepoint that already existed
    for (s in 2:l) {
      intPar.sum <- sapply(1:nmodel, FUN=function(m) {
        hmm$intPar[[i]][[m]][, s] })
      logPosterior.temp[s, ] <- intPar.sum + log(hmm$pModels) + logMAP.list[[pos[s]+1]]
      logMAP.temp[s, ] <- logPosterior.temp[s, ] + log(runl.fun(i-pos[s], ptr, msl)$gx)
    }
    # for the new changepoint tau=i-msl
    intPar.sum <- sapply(1:nmodel, FUN=function(m) {
      hmm$intPar[[i]][[m]][, 1] })
    logPosterior.temp[1, ] <- intPar.sum + log(hmm$pModels) + logMAP.list[[i-msl+1]]
    logMAP.temp[1, ] <- logPosterior.temp[1, ] + log(runl.fun(msl, ptr, msl)$gx)
    
    # record the maximum
    maxID <- which.max(logMAP.temp)
    maxl <- if ( maxID %% l == 0 ) { l } else { maxID %% l }
    tstar[i] <- pos[maxl]    # track back to the actual time point 
    mstar[i] <- ceiling(maxID / l)
    logPosterior.list[[i]] <- logPosterior.temp
    logMAP.list[[i+1]] <- logMAP.temp[maxID]
  }
  
  # trace back to find the MAP
  t.list <- m.list <- numeric(0)
  last <- hmm$Nx
  while (last > msl) {
    t.list <- c(tstar[last], t.list)
    m.list <- c(mstar[last], m.list)
    last <- tstar[last]
  }
  if (last > 0) {
    t.list <- c(0, t.list)
    m.list <- c(1, m.list)
  }
  
  return(list(breaks=t.list, models=m.list, logMAP.list=logMAP.list))
}

