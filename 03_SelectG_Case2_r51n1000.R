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
load('data_cens_case2_r51n1000.RData')

#######################################################################

Ran_trails=32
replicate=length(data_cens_case2)

start=Sys.time()

no_cores=detectCores(logical = TRUE)  
# returns the number of available hardware threads, 
# and if it is FALSE, returns the number of physical cores


SelectG_case2=list()
for(g in 1:6){
    Solutions=list()
    for(i in 1:replicate){
          
        cl=makeCluster(no_cores)  
        registerDoParallel(cl)  

        solutions=foreach(k=1:Ran_trails, .errorhandling="pass") %dopar% {
            library(condMVNorm)
            library(MomTrunc)
            MixCenMVReg_EM(Y=data_cens_case2[[i]]$Y,
                           X=data_cens_case2[[i]]$X,
                           C=data_cens_case2[[i]]$censorID, 
                           G=g,
                           Max.iter=1000, diff.tol=1e-3, print=FALSE, calc_cov=FALSE)
        }
        stopCluster(cl)
        end=Sys.time()
        message(paste0("G = ",g,", i = ",i,", Time spent: ",round(end-start,2),units(end-start)))
        Solutions[[i]]=solutions
    }
    SelectG_case2[[g]]=Solutions
    
    save(SelectG_case2,file='SelectG_case2_r51n1000.RData')
}
    

# save(model1.1,file = "model1.1.RData")
print(Sys.time()-start)

#######################################################################
