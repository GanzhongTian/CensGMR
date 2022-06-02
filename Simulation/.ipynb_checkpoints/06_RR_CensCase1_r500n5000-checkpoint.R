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
load('data_cens_case1_r500n5000.RData')
#######################################################################


start=Sys.time()

no_cores=detectCores(logical = TRUE)  
# returns the number of available hardware threads, 
# and if it is FALSE, returns the number of physical cores
# print(no_cores)

cl=makeCluster(no_cores)  
registerDoParallel(cl)  

replicate=length(data_cens_case1)
# replicate=500

case1_cens_RR=foreach(i=1:replicate, .errorhandling="remove") %dopar% {
  
    library(condMVNorm)
    library(MomTrunc)
    MixCenMVReg_EM(Y=data_cens_case1[[i]]$Y,
                   X=data_cens_case1[[i]]$X,
                   C=data_cens_case1[[i]]$censorID, 
                   G=3, pie_hat=true_pars$PIE, beta_hat=true_pars$BETA, sigma_hat=true_pars$SIGMA,
                   Max.iter=5000, diff.tol=1e-4, print=FALSE)
}
stopCluster(cl)
# save(model1.1,file = "model1.1.RData")
print(Sys.time()-start)

save(case1_cens_RR,file=paste0("case1_cens_RR_r",replicate,"n5000.RData"))