library(nimble)

 EMIRT <- nimbleCode({

# Likelihood 
    # Y[i,j]        :is the j-th response of the i-th student
    # class_p is a data object that is passed and is based on a given threshold rule (common-k or CUMP)
    for(i in 1:I){
        for(j in 1:J){
        Y[i,j] ~ dbern(correct_resp_prob[i,j])

        correct_resp_prob[i,j] <- class_p[i,j]*(irp[i,j]) + (1-class_p[i,j])*g
        irp[i,j] <- ilogit(1.7*(theta[i] - b[j]))
        }
    # Prior for theta 
    theta[i] ~ dnorm(mean=0,sd=1)
    }

      #Prior for other paramters
    for(j in 1:J){

        b[j] ~ dnorm(mean=0, sd=1.0)
    }
    g ~ T(dbeta(1,4),,0.3)
})