# Include necessary packages
# library(tidyverse)

# library(mvtnorm)
library(condMVNorm)
library(MomTrunc)
# library(truncnorm)
# library(tmvtnorm)
# library(TTmoment)


library(doParallel)
# library(parallel)

#######################################################################

path=getwd()
source(paste0(substring(path, 
                        1, 
                        tail(unlist(gregexpr(pattern ='/',path)),n=1)),
              "MixCenMVReg_EM.R"))

source(paste0(substring(path, 
                        1, 
                        tail(unlist(gregexpr(pattern ='/',path)),n=1)),
              "Util_Func.R"))


load('true_pars.RData')
load('data_cens_case1_r51n1000.RData')
B=true_pars$BETA
for(i in 1:length(B)){
    B[[i]]=t(matrix(B[[i]][1,]))
}
#######################################################################

start=Sys.time()

no_cores=detectCores(logical = TRUE)  
# returns the number of available hardware threads, 
# and if it is FALSE, returns the number of physical cores
print(no_cores)

cl=makeCluster(no_cores)  
registerDoParallel(cl)  

replicate=length(data_cens_case1)


case1_gmm_bias=foreach(i=1:replicate, .errorhandling="pass") %dopar% {
  
    library(condMVNorm)
    library(MomTrunc)
    MixCenMVReg_EM(Y=data_cens_case1[[i]]$Y,
                   X=as.matrix(data_cens_case1[[i]]$X[,1]),
                   C=data_cens_case1[[i]]$censorID, 
                   G=3, pie_hat=true_pars$PIE, beta_hat=B, sigma_hat=true_pars$SIGMA,
                   Max.iter=1000, diff.tol=1e-3, print=FALSE, calc_cov=FALSE)
}
stopCluster(cl)
# save(model1.1,file = "model1.1.RData")
print(Sys.time()-start)

save(case1_gmm_bias,file='case1_gmm_bias_r51n1000.RData')

#######################################################################

