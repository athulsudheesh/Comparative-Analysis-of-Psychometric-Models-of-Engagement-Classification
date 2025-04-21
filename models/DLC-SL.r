
library(nimble) 
DLCSL <- nimbleCode({

# Likelihood 
    # Y[i,j]        :is the j-th response of the i-th student
    # C[i,j]        :is the i-th student's latent engagement state [0/1]
    #                while responding to the j-th item
    # class_p[i,j]  : probability of the i-th student being engaged for j-th item
    # irp[i,j]      : Probability of i-th student getting j-th item correct 
    # g             : random guessing probability (constant)
    # tau[j]           : common time threshold common to all students  
    for(i in 1:I){
        for(j in 1:J){
        Y[i,j] ~ dbern(correct_resp_prob[i,j])
        C[i,j] ~ dbern(class_p[i,j])

        correct_resp_prob[i,j] <- class_p[i,j]*(irp[i,j]) + (1-class_p[i,j])*g
        irp[i,j] <- ilogit(1.7*(theta[i] - b[j]))
        class_p[i,j] <- ilogit(gamma*(T_mat[i,j] - tau[j]))
        }
    # Prior for theta 
    theta[i] ~ dnorm(mean=0,sd=1)
    }




      #Prior for other paramters
    for(j in 1:J){
        #g[j] ~ dbeta(1,4)T(,0.3)
        #a[j] ~ T(dnorm(mean=0, sd=1.0),0,)
        b[j] ~ dnorm(mean=0, sd=1.0)
        tau[j] ~ dnorm(mean=0, sd=1.0)
    }
    gamma ~ T(dnorm(mean=0,sd=1.0),0,)
    g ~ T(dbeta(1,4),,0.3)
})