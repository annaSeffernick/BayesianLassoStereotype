#########################################################
# Bayesian Stereotype Model Example Code 
# Anna Eames Seffernick
# February 9, 2022
########################################################

# Note you must have JAGS installed. See https://mcmc-jags.sourceforge.io

# Load required packages
library(GEOquery)
library(affy)
library(rjags)
library(dclone)
library(VGAM)

# Load the data, called small.data.scale
load("BLASSO-ST_SmallData.RData")

################
# Fit the Model
################

# write a function to perform MCMC with parallelization
# credit: Han Fu provided the original code
coda_samples <- function(dataList, model.file, parallel=T){
  if(parallel) {
    cl <- parallel::makePSOCKcluster(nChains)
    dclone::parJagsModel(cl = cl, name = "res", file = model.file, data = dataList, n.chains = nChains, 
                         inits = list(inits1, inits2, inits3), quiet=FALSE, n.adapt=adaptSteps)
    cat( "Burning in the MCMC chain...\n")
    dclone::parUpdate(cl = cl, object = "res", n.iter = burnInSteps)
    cat("Sampling final MCMC chain...\n" )
    codaSamples <- dclone::parCodaSamples(cl = cl, model = "res", 
                                          variable.names=c("beta","phi","alpha", "lambda", "gamma", "betgma"), 
                                          n.iter=nIter, 
                                          thin=thinSteps)
  }
  else{
    model = rjags::jags.model(file = model.file, data = dataList, n.chains = nChains, quiet=FALSE,
                              inits = list(inits1, inits2, inits3), n.adapt=adaptSteps)
    cat( "Burning in the MCMC chain...\n")
    update(model, n.iter=burnInSteps)  
    cat("Sampling final MCMC chain...\n" )
    codaSamples <- rjags::coda.samples(model, 
                                       variable.names=c("beta", "phi", "alpha", "lambda", "gamma", "betgma"), 
                                       n.iter=nIter, 
                                       thin=thinSteps)
  }
  return(codaSamples)
}


# Write the model file
Model = "model{
  for(i in 1:N){
    mu[i] <- inprod(X[,i], betgma[])
    for(j in 1:3){
      pi[i,j] <- exp(alpha[j] + phi[j]*mu[i])
    }
    Y[i] ~ dcat(pi[i, 1:3])
  }
  for(b in 1:2){
    alpha[b] ~dnorm(0.0, 0.04)
  }
  alpha[3] <- 0
  for(k in 1:P){
    beta[k] ~ ddexp(0, lambda)
  }
  lambda ~ dgamma(0.1, 0.1)
  for (j in 1:P){
  	betgma[j] <- beta[j]*gamma[j]
  	gamma[j] ~ dbern(0.3)
  }
  phi[1] <- 1
  phi[2] ~ dunif(0.0, 1.0)
  phi[3] <- 0
}"

JAGSFILE="BLassoVIFixedUnifGammaPoint3.bug"
cat(Model, file=JAGSFILE)

##Set Model Parameters
N <- dim(small.data.scale)[2]
P <- dim(small.data.scale)[1]
#Data List
dataList <- list("Y" = pData(small.data.scale)$RiskGroup, "X" = exprs(small.data.scale), "N" = N, "P" = P)
# Parameters to be monitored
parameters <- c("beta", "phi", "alpha", "lambda", "gamma", "betgma")

# JAGS Set-up
adaptSteps <- 1000              #number of steps to "tune" the samplers
burnInSteps <-3000             #number of steps to "burn-in" the samplers
nChains <- 3                    #number of chains to run
numSavedSteps <- 9999        #total number of steps in chains to save
thinSteps <- 3                     #number of steps to "thin" (1=keep every step)
nIter <- ceiling((numSavedSteps*thinSteps )/nChains)    #steps per chain
# Set initial values
fit.st <- vglm(dataList$Y ~ 1, multinomial(refLevel="Adverse"))
alpha.est.vec <- coef(summary(fit.st))[,1]
names(alpha.est.vec) <- NULL
set.seed(16)
inits1 <- list("beta" = rep(0.0, P),"alpha" = c(alpha.est.vec, NA), "phi" = c(NA, runif(1, 0.0,1.0), NA),
               "lambda" = rgamma(1, shape=0.1, rate=0.1))
set.seed(23)
inits2 <- list("beta" = rep(0.0, P),"alpha" = c(alpha.est.vec + 0.025, NA), "phi" = c(NA, runif(1, 0.0, 1.0), NA),
               "lambda" = rgamma(1, shape=0.1, rate=0.1))
set.seed(37)
inits3 <- list("beta" = rep(0.0, P),"alpha" = c(alpha.est.vec + 0.05, NA), "phi" = c(NA, runif(1, 0.0, 1.0), NA),
               "lambda" = rgamma(1, shape=0.1, rate=0.1))
# Find Posterior Samples
set.seed(25) # with parallelization may not reproduce exactly
codaSamples <- coda_samples(dataList, JAGSFILE, parallel=T)

############################
# Explore Posterior Samples
############################

# Trace plots for select parameters
plot(codaSamples[,1], main="alpha1")
plot(codaSamples[,2], main="alpha2")
plot(codaSamples[,7], main="beta7")
plot(codaSamples[,16], main="beta3 x gamma3")
plot(codaSamples[,30], main="gamma7")
plot(codaSamples[,34], main="lambda")
plot(codaSamples[,36], main="phi2")


# Mixing looks pretty good for most of the parameters
gelman.list <- gelman.diag(codaSamples, multivariate=FALSE) 
length(which(gelman.list$psrf[,1]>1.1))
length(which(gelman.list$psrf[,2]>1.1)) # 0 parameters fail to converge
gelman.list$psrf

mcmcChain <- as.matrix(codaSamples)

#alpha_1
boxplot(mcmcChain[,1])
summary(mcmcChain[,1]) 

# alpha_2
boxplot(mcmcChain[,2])
summary(mcmcChain[,2]) 


# phi_2 
boxplot(mcmcChain[,36])
summary(mcmcChain[,36]) 
    

# Lambda
boxplot(mcmcChain[,34])
summary(mcmcChain[,34]) 


# Find significant transcripts

# Use credible intervals
mcmcChainGamma <- mcmcChain[, 24:33]
quantile(mcmcChainBetaGamma[,1], probs=c(0.025, 0.975)) 

# Look for credible intervals that do not contain zero
betagma.ciL <- c()
betagma.ciU <- c()
for(i in 1:10){
  betagma.ciL[i] <- quantile(mcmcChainBetaGamma[,i], probs=0.025)
  betagma.ciU[i] <- quantile(mcmcChainBetaGamma[,i], probs=0.975)
}
est.sig.feats <- c()
for(i in 1:10){
  est.sig.feats[i] <- ifelse(sign(betagma.ciL[i])==sign(betagma.ciU[i]), 1, 0)
}
plot(est.sig.feats)
length(which(est.sig.feats==1)) # identifies 3 significant transcripts
which(est.sig.feats==1) 

# Try highest posterior density intervals
library("HDInterval")
test.hdi <- hdi(mcmcChainBetaGamma[,1], credMass=0.95) 
test.hdi

betagma.hdi <- list()
for(i in 1:10){
  betagma.hdi[[i]] <- hdi(mcmcChainBetaGamma[,i], credMass=0.95)
}
est.sig.feats.HDI <- c()
for(i in 1:10){
  est.sig.feats.HDI[i] <- ifelse(sign(betagma.hdi[[i]][1])==sign(betagma.hdi[[i]][2]), 1, 0)
}
plot(est.sig.feats.HDI)
length(which(est.sig.feats.HDI==1)) # This identifies 3 significant transcripts 
which(est.sig.feats.HDI==1)

# Try using Bayes Factors for VI
# Use epsilon=0.1
a <- 0.1
b <- 0.1
epsilon <- 0.1
pgamma1 <- 0.3
pgamma0 <- 1 - pgamma1
prior.odds <- (pgamma1*(b^a)*gamma(a+1))/(pgamma1*((b+epsilon)^a - b^a)*gamma(a+1) + pgamma0*a*((b+epsilon)^a)*gamma(a))
post.odds <- c()
BF <- c()
for(i in 1:dim(mcmcChainBetaGamma)[2]){
  betagma.test <- mcmcChainBetaGamma[,i]
  P.Beta.G.Eps <- sum(ifelse(abs(betagma.test) > epsilon, 1, 0))/dim(mcmcChainBetaGamma)[1]
  P.Beta.LE.Eps <- sum(ifelse(abs(betagma.test)<= epsilon, 1, 0))/dim(mcmcChainBetaGamma)[1]
  post.odds[i] <- P.Beta.G.Eps/P.Beta.LE.Eps
  BF[i] <- post.odds[i]/prior.odds
}
length(which(BF > 1)) # 3 transcripts have BF > 1
length(which(BF>5)) # 3 transcripts
which(BF>5) 
length(which(BF>10)) # 3 transcripts
which(BF>10) 

## Try variable selection focused on the gammas
post.mean.gamma <- colMeans(mcmcChainGamma)
length(which(post.mean.gamma > 0.5)) # 3 transcripts with posterior mean > 0.5
which(post.mean.gamma > 0.5)
sort(post.mean.gamma, decreasing=TRUE)[1:10] #454, 226, 171, 779, 144, 497, 417, 480, 100, 2

# Try Bayes factor for gamma 
# Testing H_0 gamma = 0 vs. Ha gamma=1

prior.odds <- pgamma1/pgamma0
post.odds <- c()
BF <- c()
for(i in 1:dim(mcmcChainGamma)[2]){
  gamma.test <- mcmcChainGamma[,i]
  P.Gamma.1 <- mean(gamma.test)
  P.Gamma.0 <- 1-mean(gamma.test)
  post.odds[i] <- P.Gamma.1/P.Gamma.0
  BF[i] <- post.odds[i]/prior.odds
}
length(which(BF > 1)) # 3 transcripts have BF > 1
length(which(BF>5)) # 3 transcripts
length(which(BF>10)) # 3 transcripts
which(BF>10) 

# Same 3 transcripts selected by  each method, extract significant transcript names
featureNames(small.data.scale)[which(BF>10)]

#######################
# Session Info
#######################

#R version 4.1.2 (2021-11-01)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Big Sur 10.16

#Matrix products: default
#LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

#locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#  [1] splines   stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] HDInterval_0.2.2    caret_6.0-90        lattice_0.20-45     ggplot2_3.3.5       VGAM_1.1-5          dclone_2.3-0       
#  [7] Matrix_1.4-0        rjags_4-12          coda_0.19-4         lubridate_1.8.0     affyio_1.64.0       stringr_1.4.0      
#  [13] purrr_0.3.4         digest_0.6.29       affy_1.72.0         GEOquery_2.62.2     Biobase_2.54.0      BiocGenerics_0.40.0
