MixCenUVReg_EM=function(Y, C, X, G=2, Max.iter=1000, 
                        pie_hat=NA, beta_hat=NA, sigma_hat=NA, diff.tol=1e-3,
                        print=TRUE, init_class=NA,calc_cov=FALSE){
    #data must be 1-dimensional vector
    #C must be 1-dimensional vector of -1,0,1 (left, no, right)
    ##############################################
    # input:    -Y    :a cont 1-dimensional vector of length n, 
    #           -X    :a design matrix containing the covariates of dimension nxp 
    #           -C    :a 1-dimensional vector of length n, taking -1,0,1 (left, no, right)
    #           
    # output:   a list of parameters and posteriors
    ##############################################

    #N=dim(X)[1]
    N=dim(Y)[1]
    P=dim(Y)[2]
    D=dim(X)[2]
    
    if(is.na(sum(as.integer(init_class)))==FALSE){
        G=length(levels(init_class))
    }
    
    colnames_Y=colnames(Y)
    
    if(is.null(colnames(Y))){
        colnames(Y)=paste(rep('Y',P), 1:P, sep="_")
    }
    
    if(all(dim(Y)==dim(C)) & P==1){
        Y=Y[,1]
        C=C[,1]
        
        interval=cbind(rep(-Inf, N),rep(Inf, N))
        for(i in 1:N){
            if(C[i]==-1){
            interval[i,]=c(-Inf, Y[i])
            }
            if(C[i]==1){
                interval[i,]=c(Y[i],Inf)
            }

        }
    }else{
        stop("Y and C must be nx1 matrice of the same size")
    }

    
    initial=FALSE
    
    if(all(is.na(pie_hat)==FALSE) & all(is.na(beta_hat)==FALSE) & all(is.na(sigma_hat)==FALSE)){

        initial=TRUE

        if(length(beta_hat)==length(sigma_hat) & length(beta_hat)==length(pie_hat)){
           G=length(beta_hat)
        }else {
           stop("Intial pie, beta, sigma lengths are different")
        }

        
        mu_hat=list()
        for(g in 1:G){
          mu_hat[[g]]=X%*%beta_hat[[g]]
        }
      }


    all_obs.LogLik=c(-Inf)
    

    # set initial iteration number:
    iter=0
    diff=Inf
    # set initial conditional probabilities:

    log.ind_density=matrix(NA,nrow=N,ncol=G) #individual contribution to likelihood (NxK)

#     while (Max.iter>=iter & abs(all_obs.LogLik[length(all_obs.LogLik)]-all_obs.LogLik[length(all_obs.LogLik)-1])>=(N*diff.tol) ){
    while(iter<Max.iter & diff>diff.tol){
        
        if(iter==0 & initial==FALSE){
      
      ## E-step: computing the conditional posterior probabilities tau_hat
      
      ### First assign every one to be 1/G or a very small value
      #             tau_hat=matrix(rep(1/G,N*G), nrow=N, ncol=G)
      # tau_hat=matrix(rep(0,N*G), nrow=N, ncol=G)
      tau_hat=matrix(rep(1e-3,N*G), nrow=N, ncol=G)
                
      if(is.na(sum(as.integer(init_class)))){

          ### Then randomly assign every one to be 1/G
          subsample=sample(1:N,round(N/10),replace=FALSE)
          tau_hat[subsample,]=t(rmultinom(round(N/10), 1, rep(1,G)/G))
          
      }else{
          ### Assign the initial class label if the label is provided
          for(g in 1:G){
              tau_hat[init_class==(levels(init_class)[g]),g]=1
          }
          
      }
      pie_hat=apply(tau_hat,2,mean)
        
      mu_hat=list()
      beta_hat=list()  
      sigma_hat=list()  
      for(g in 1:G){
        # mu_hat[[g]]=tau_hat[,g]*Y
          
        beta_hat[[g]]=solve(t(X)%*%diag(tau_hat[,g])%*%X)%*%t(X)%*%diag(tau_hat[,g])%*%Y
        mu_hat[[g]]=X%*%beta_hat[[g]]  
        # sigma_hat[[g]]=cov(tau_hat[,g]*Y)
        sigma_hat[[g]]=cov.wt(as.matrix(Y-mu_hat[[g]]),tau_hat[,g])$cov[1,1]
      }     
      
    }else{
        
        # print(mu_hat)
        # print(sigma_hat)
        # E step: computing the conditional probabilities
        for(g in 1:G){
            log.ind_density[,g]=log(pie_hat[g])+(dnorm(Y,mu_hat[[g]],sigma_hat[[g]],log=T)*(C==0)
                                                 +pnorm(Y,mu_hat[[g]],sigma_hat[[g]],log.p=T)*(C==-1)
                                                 +pnorm(Y,mu_hat[[g]],sigma_hat[[g]],log.p=T,lower.tail=F)*(C==1))
        }

        log.ind_density[is.infinite(log.ind_density)]=-999                
            
        new_tau_hat=exp(sweep(log.ind_density, 1, apply(log.ind_density,1,logSumExp))) #(NxG)
                
        # new_tau_hat=ind_density/apply(ind_density,1,sum) #(NxG)

        # new_tau_hat[is.na(new_tau_hat)]=1/G # adjust the nans 
        obs.LogLik=sum(apply(log.ind_density,1,logSumExp))
                
        all_obs.LogLik=append(all_obs.LogLik,obs.LogLik)

        # if(iter>=1){
        #     diff.posterior=max(abs(tau_hat-new_tau_hat))
        # }
#         print(diff.posterior)
        tau_hat=new_tau_hat


        ## M-step: Update pie_hat
        new_pie_hat=apply(tau_hat,2,mean)
        diff=max(abs(pie_hat-new_pie_hat))
        pie_hat=new_pie_hat
        
        # Update beta_hat
        Y_star=list()
            for(g in 1:G){
        # comput truncated normal 1st moment
                Y_star[[g]]=((C==0)*Y
                             +(C==-1)*etruncnorm(a=-Inf, b=Y, mean=mu_hat[[g]], sd=sigma_hat[[g]])
                             +(C==1)*etruncnorm(a=Y, b=Inf, mean=mu_hat[[g]], sd=sigma_hat[[g]]))


    #             Y_star=mapply(eval_ystar, Y, C, mu_hat[[g]], sigma_hat[[g]])
    #         mu_hat[[i]]=apply(tau_hat[,i] * as.matrix(Y_star),2,sum) / sum(tau_hat[,i])
                new_beta_hat=solve(t(X)%*%diag(tau_hat[,g])%*%X)%*%t(X)%*%diag(tau_hat[,g])%*%Y_star[[g]]
                diff=max(diff,max(abs(beta_hat[[g]]-new_beta_hat)))
                beta_hat[[g]]=new_beta_hat

                mu_hat[[g]]=X%*%beta_hat[[g]]
            }

        # Update sigma_hat
        
            eval_s=function(y,c,mu,sigma){
                if(c==0){
                    return((y-mu)^2)
                }
                if(c==1){
                    return((etruncnorm(a=y, b=Inf, mean=mu, sd=sigma)-mu)^2+v_censnorm(c(y,Inf), mu, sigma))
                }          
                if(c==-1){
                    return((etruncnorm(a=-Inf, b=y, mean=mu, sd=sigma)-mu)^2+v_censnorm(c(-Inf,y), mu, sigma))
                }
            }
            S_star=list()
            for(g in 1:G){

        #         s=((C==0)*(Y-mu_hat[[g]])^2
        #            +(C==-1)*((etruncnorm(a=-Inf, b=Y, mean=mu_hat[[g]], sd=sigma_hat[[g]])-mu_hat[[g]])^2
        #                              +V_CENSNORM(interval, mean=mu_hat[[g]], sd=sigma_hat[[g]]))
        #            +(C==1)*((etruncnorm(a=Y, b=Inf, mean=mu_hat[[g]], sd=sigma_hat[[g]])-mu_hat[[g]])^2
        #                              +V_CENSNORM(interval, mean=mu_hat[[g]], sd=sigma_hat[[g]])))
#             eval_s=function(y,c,mu,sigma){
#                 if(c==0){
#                     return((y-mu)^2)
#                 }
#                 if(c==1){
#                     return((etruncnorm(a=y, b=Inf, mean=mu, sd=sigma)-mu)^2+v_censnorm(c(y,Inf), mu, sigma))
#             }          
#                 if(c==-1){
#                     return((etruncnorm(a=-Inf, b=y, mean=mu, sd=sigma)-mu)^2+v_censnorm(c(-Inf,y), mu, sigma))
#                 }
#             }

                S_star[[g]]=mapply(eval_s,Y,C,as.vector(mu_hat[[g]]),sigma_hat[[g]])

                new_sigma_hat=(apply(tau_hat[,g] * as.matrix(S_star[[g]]),2,sum)/ sum(tau_hat[,g]))^0.5
                diff=max(diff,max(abs(sigma_hat[[g]]-new_sigma_hat)))
                sigma_hat[[g]]=new_sigma_hat
            }
        }
        iter = iter + 1

        if(min(unlist(sigma_hat))<1e-10){
            message("EM stopped becaused of degenerating solution!")
            break
        }
        
    }    


    
    # fisher INFO
#     score_pie=list()
#     for(g in 1:(G-1)){
#         score_pie[[g]]=tau_hat[,g]/pie_hat[g]-tau_hat[,G]/pie_hat[G]
#     }                
#     score_pie=do.call('cbind',score_pie)
#     colnames(score_pie)=combine(c("PIE"),1:(G-1))
    if(calc_cov==TRUE){
        score_beta=list()
        for(g in 1:G){
            score_beta[[g]]=diag(as.vector(tau_hat[,g]*(Y_star[[g]]-mu_hat[[g]])))%*%X/(sigma_hat[[g]])^2
            colnames(score_beta[[g]])=combine(colnames_Y,colnames(X),c(g))
        }
        score_beta=do.call('cbind',score_beta)    
        
        all_score=cbind(score_beta)
    
        obs_Info=t(all_score)%*%all_score
        
        if(any(eigen(obs_Info)$values<=1e-6)){
            cov_notPD=TRUE
            print('The Cov matrix for coefficients is not positive definite')
        }else{
            cov_notPD=FALSE
            cov=solve(obs_Info)
      }
    
    }else{
        obs_Info=NA
        cov_notPD=TRUE
        print('The Info matrix for coefficients is not computed')
    }

    
#     score_sigma=list()
#     for(g in 1:G){
#         score_sigma[[g]]=tau_hat[,g]*(S_star[[g]]-(sigma_hat[[g]])^2)/(2*(sigma_hat[[g]])^4)
#     }
#     score_sigma=do.call('cbind',score_sigma)
    
#     all_score=cbind(score_pie,score_beta,score_sigma)

    
    
    if(iter<Max.iter & diff<=diff.tol){
        converge_=TRUE
    }else{
        converge_=FALSE
    }
    
    names(pie_hat)=paste0("pie",1:G)
    names(beta_hat)=paste0("beta",1:G)
    names(sigma_hat)=paste0("sigma",1:G)
    if(print){
        print(paste0("Total Iteration = ",toString(iter)))
        print(paste0("Convergence = ",toString(converge_)))
        print(paste0("LogLik = ",toString(round(tail(all_obs.LogLik, n=1),3))))
        print(pie_hat)
        print(beta_hat)
        print(sigma_hat)
        #         print(obs_Info)
        plot(all_obs.LogLik[-1],xlab="Iterations", ylab="Log-Likelihood",)    
    }

    LogLik=tail(all_obs.LogLik, n=1)
    nparm=length(as.vector(unlist(beta_hat)))+length(unique(as.vector(unlist(sigma_hat))))+length(pie_hat)-1

    POST_PROB=tau_hat
    # POST_PROB[POST_PROB < 1e-8] = 1e-8
    entropy=-sum(POST_PROB*log(POST_PROB))
    entropy.rsquare=1 - entropy/(N*log(G))    

    AIC=-2*LogLik+nparm*2
    BIC=-2*LogLik+nparm*log(N)
    ICL=BIC+2*entropy

    OUTPUT=list()

    OUTPUT[[1]]=iter
    OUTPUT[[2]]=converge_
    OUTPUT[[3]]=nparm
    OUTPUT[[4]]=LogLik
    OUTPUT[[5]]=AIC
    OUTPUT[[6]]=BIC
    OUTPUT[[7]]=ICL

    OUTPUT[[8]]=pie_hat
    OUTPUT[[9]]=beta_hat
    OUTPUT[[10]]=sigma_hat
    OUTPUT[[11]]=tau_hat
    OUTPUT[[12]]=apply(tau_hat,1,eval_class)
    
    if(cov_notPD){
        OUTPUT[[13]]=obs_Info
        names(OUTPUT)<-c("Iterations", "Converged", "df", "LogLik","AIC", "BIC", "ICL","Pie", "Beta", "Sigma", "Posterior","Class","obs_Info")

    }else{
        OUTPUT[[13]]=cov
        names(OUTPUT)<-c("Iterations", "Converged", "df", "LogLik","AIC", "BIC", "ICL","Pie", "Beta", "Sigma", "Posterior","Class","Cov")
    }
    
    return(OUTPUT)
}
