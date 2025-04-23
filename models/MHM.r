library(nimble)
MHM <- nimbleCode({
    mu_person[1:2] <- rep(0,2)
    # Likelihood
    for (i in 1:I){
        for(j in 1:J){
            Y[i,j] ~ dbern(correct_resp_prob[i,j])
            C[i,j] ~ dbern(class_p[i,j])
            T_mat[i,j] ~ dnorm(mu_l_ij[i,j], tau_time[i,j])
            
            correct_resp_prob[i,j] <- class_p[i,j]*(irp[i,j]) + (1-class_p[i,j])*g
            irp[i,j] <- ilogit(1.7*(theta[i] - b[j]))
            
            class_p[i,j] <- ilogit(1.7*(eta[i] - kappa[j]))
            
            # Mean structure is different from LIC-RI-IRT
            mu_l_ij[i,j] <- nu_til + C[i,j]*(delta[j] + lambda)
            sigma2_time[i,j] <- C[i,j]*sigma2_epsilon[j] + (1-C[i,j])*sigma2_epsilon_til
            tau_time[i,j] <- pow(sigma2_time[i,j],-1)
            
        }
        # Prior for person parameter
        person_par[i,1:2] ~ dmnorm(mu_person[1:2], cholesky = Sigma_person[1:2,1:2], prec_param=0)
        theta[i] <- person_par[i,1]
        eta[i]   <- person_par[i,2]
        #xi[i]    <- person_par[i,3]
    }
      #Prior 
    for(j in 1:J){
        #g[j] ~ dbeta(1,1)T(,0.3)
        a[j] ~ T(dnorm(0, 1.0),0,)
        b[j] ~ dnorm(0, 1.0)
        kappa[j] ~ dnorm(0, 1.0)
        
        delta[j] ~ T(dnorm(0, 1.0),0,)
        
        sigma2_epsilon[j] ~ dgamma(1/2,1/2)
    }
    nu_til ~ T(dnorm(0, 1.0),,0)
    sigma2_epsilon_til ~ dgamma(1/2,1/2)
    g ~ T(dbeta(1,1),,0.3)
    lambda ~ T(dnorm(0,1.0),0,)
    Sigma_person[1:2,1:2] ~ dlkj_corr_cholesky(1,2)
})