library(nimble)
DLCTL <- nimbleCode({
    # Y[i,j]        :is the j-th response of the i-th student
    # C[i,j]        :is the i-th student's latent engagement state [0/1]
    # zeta[i]       : individual specific time threshold
    # class_p[i,j]  : probability of the i-th student being engaged for j-th item
    # irp[i,j]      : Probability of i-th student getting j-th item correct 
    # g             : random guessing probability (constant)
    # tau[j]           : common time threshold common to all students  
    
    # mu_person is the mean of the multivariate normal distribution from which the prior for theta and zeta is derived
    # mu_person is set to be 0 for identifiability reasons (refer Yamaguchi & Fujita, 2023 p 78 for more detials)

    mu_person[1:2] <- rep(0,2)
   for (i in 1:I){
        for(j in 1:J){
            Y[i,j] ~ dbern(correct_resp_prob[i,j])
            C[i,j] ~ dbern(class_p[i,j])
    
            #correct_resp_prob[i,j] <- C[i,j]*(irp[i,j]) + (1-C[i,j])*g[j]
            correct_resp_prob[i,j] <- class_p[i,j]*(irp[i,j]) + (1-class_p[i,j])*g
            irp[i,j] <- ilogit(1.7*(theta[i] - b[j]))
            class_p[i,j] <- ilogit(gamma*(T_mat[i,j] - (tau[j] + zeta[i])))
        }
        # Prior for person parameter
        person_par[i,1:2] ~ dmnorm(mu_person[1:2], cholesky =Sigma_person[1:2,1:2], prec_param=0)
        theta[i] <- person_par[i,1]
        zeta[i] <- person_par[i,2]
    } 
      #Prior 
    for(j in 1:J){
        #g[j] ~ dbeta(1,1)T(,0.3)
        a[j] ~ T(dnorm(0, 1.0),0,)
        b[j] ~ dnorm(0, 1.0)
        tau[j] ~ dnorm(0, 1.0)
    }
    gamma ~ T(dnorm(0,1.0),0,)
    g ~ T(dbeta(1,1),,0.3)
    Sigma_person[1:2,1:2] ~ dlkj_corr_cholesky(1,2)
})