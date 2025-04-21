files.sources = list.files("models", full.names=TRUE)
sapply(files.sources, source)

#===================================================================================================
#                                          Latent Variable Approaches
#===================================================================================================
#===========================================================================
#                            DLC-SL
#===========================================================================
DLCSL.fit <- function(dat, 
                    params.to.monitor=NULL,
                    mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 3
                                    )){
    X <- dat$X
    RT <- dat$RT

    RT[is.na(RT)] <- 0
    RT <- log(RT)
    I <- nrow(X)
    J <- ncol(X)
    const <- list(I=I,J=J)
    data <- list(Y=X,T_mat = RT)
    inits <- function(){
    list(theta=rnorm(I,0,1),
        a = abs(rnorm(J,0,1)),
        b = rnorm(J,0,1),
        #g = runif(J,0,0.3),
        g = runif(1,0,0.3),
        tau = abs(rnorm(J,0,1)),
        gamma = abs(rnorm(1,0,1)))
    }

    model <- nimbleModel(code=DLCSL, constants = const,
                    data = data, inits = list())
    
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE) 
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    nodes = list(dataNodes = dataNodes, parentNodes = parentNodes, simNodes=simNodes)
    
    mcmcfit <- nimbleMCMC(model = model,
                        monitors = c(parentNodes,params.to.monitor),
                        niter = mcmc.config$no.of.iterations,
                        nburnin = mcmc.config$no.of.burn.ins,
                        nchains = mcmc.config$no.of.chains,
                        samplesAsCodaMCMC = TRUE)
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes))
}

#============================= END OF DLC-SL ==============================


#===========================================================================
#                            DLC-TL
#===========================================================================
DLCTL.fit <- function(dat, 
                    params.to.monitor=NULL,
                    mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 3
                                    )){
    
    X <- dat$X
    RT <- dat$RT
    RT[is.na(RT)] <- 0
    RT <- log(RT)
    #RT <- scale(RT)
    I <- nrow(X)
    J <- ncol(X)
    const <- list(I=I,J=J)
    data <- list(Y=X,T_mat = RT)
    inits <- function(){
    list(person_par=matrix(rnorm(I*2,0,1),ncol=2),
        a = abs(rnorm(J,0,1)),
        b = rnorm(J,0,1),
        #g = runif(J,0,0.3),
        g = runif(1,0,0.3),
        tau = rnorm(J,0,1),
        gamma = abs(rnorm(1,0,1))
    )
    }
    model <- nimbleModel(code=DLCTL, constants = const,
                        data = data, inits = list())
    
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE) 
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    nodes = list(dataNodes = dataNodes, parentNodes = parentNodes, simNodes=simNodes)
    
    mcmcfit <- nimbleMCMC(model = model,
                            monitors = c(parentNodes,params.to.monitor),
                            niter = mcmc.config$no.of.iterations,
                            nburnin = mcmc.config$no.of.burn.ins,
                            nchains = mcmc.config$no.of.chains,
                            samplesAsCodaMCMC = TRUE)
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes))
}

#============================= END OF DLC-TL ==============================



#===========================================================================
#                            ILC-RE
#===========================================================================

ILCRE.fit <- function(dat, 
                    params.to.monitor=NULL,
                    mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 3
                                    )){
    X <- dat$X
    RT <- dat$RT
    RT[is.na(RT)] <- 0

    RT <- scale(RT)
    I <- nrow(X)
    J <- ncol(X)
    const <- list(I=I,J=J)
    data <- list(Y=X,T_mat = RT)
    Sigma <- matrix(c(1, 0.4,0,
                    0.4, 1,0.,
                    0,0.4,1),byrow = T,ncol=3)
    inits <- function(){
    list(person_par=MASS::mvrnorm(I,mu=rep(0,3),Sigma = Sigma,empirical = T),
        a = abs(rnorm(J,0,1)),
        b = rnorm(J,0,1),
        #g = runif(J,0,0.3),
        g = runif(1,0,0.3),
        kappa = rnorm(J,0,1),
        nu_til = -abs(rnorm(J,0,1)),
        delta = abs(rnorm(J,1,1)),
        lambda = abs(rnorm(J,0.5,0.5)),
        sigma2_epsilon = abs(rnorm(J,0,1)),
        sigma2_epsilon_til = abs(rnorm(1,0,1)),
        rho_theta_eta = runif(1,0.5,1),
        rho_theta_xi = runif(1,-0.5,0.5),
        rho_eta_xi = runif(1,0.5,1)
    )
    }
    model <- nimbleModel(code=ILCRE, constants = const,
                        data = data, inits = inits())
    
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE) 
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    nodes = list(dataNodes = dataNodes, parentNodes = parentNodes, simNodes=simNodes)

    mcmcfit <- nimbleMCMC(model = model,
                            monitors = c(parentNodes,params.to.monitor),
                            niter = mcmc.config$no.of.iterations,
                            nburnin = mcmc.config$no.of.burn.ins,
                            nchains = mcmc.config$no.of.chains,
                            samplesAsCodaMCMC = TRUE)
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes))                            
}
#============================= END OF ILC-RE ==============================



#===========================================================================
#                            ILC-RI
#===========================================================================

ILCRI.fit <- function(dat, 
                    params.to.monitor=NULL,
                    mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 3
                                    )){
    X <- dat$X
    RT <- dat$RT
    RT[is.na(RT)] <- 0
    RT <- scale(RT)
    I <- nrow(X)
    J <- ncol(X)
    const <- list(I=I,J=J)
    data <- list(Y=X,T_mat = RT)
    Sigma <- matrix(c(1, 0.4,0,
                    0.4, 1,0.,
                    0,0.4,1),byrow = T,ncol=3)
    inits <- function(){
    list(person_par=MASS::mvrnorm(I,mu=rep(0,3),Sigma = Sigma,empirical = T),
        a = abs(rnorm(J,0,1)),
        b = rnorm(J,0,1),
        #g = runif(J,0,0.3),
        g = runif(1,0,0.3),
        kappa = rnorm(J,0,1),
        nu_til = -abs(rnorm(J,0,1)),
        delta = abs(rnorm(J,1,1)),
        lambda = abs(rnorm(J,0.5,0.5)),
        sigma2_epsilon = abs(rnorm(J,0,1)),
        sigma2_epsilon_til = abs(rnorm(1,0,1)),
        rho_theta_eta = runif(1,0.5,1),
        rho_theta_xi = runif(1,-0.5,0.5),
        rho_eta_xi = runif(1,0.5,1)
    )
        }
    
    model <- nimbleModel(code=ILCRI, constants = const,
                        data = data, inits = list())
    
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE) 
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    nodes = list(dataNodes = dataNodes, parentNodes = parentNodes, simNodes=simNodes)

    mcmcfit <- nimbleMCMC(model = model,
                           monitors = c(parentNodes,params.to.monitor),
                            niter = mcmc.config$no.of.iterations,
                            nburnin = mcmc.config$no.of.burn.ins,
                            nchains = mcmc.config$no.of.chains,
                            samplesAsCodaMCMC = TRUE)

    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes)) 
}

#============================= END OF ILC-RI ==============================

#===================================================================================================
#                                          Threshold based Approaches
#===================================================================================================

#===========================================================================
#                            Common-K
#===========================================================================
CommonK.fit <- function(dat, k=5, 
                    params.to.monitor=NULL,
                    mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 3
                                    )){
    X <- dat$X
    RT <- dat$RT
    RT[is.na(RT)] <- 0
    RT <- RT > k
    I <- nrow(X)
    J <- ncol(X)
    const <- list(I=I,J=J)
    data <- list(Y=X,class_p = RT)
    inits <- function(){
    list(theta=rnorm(I,0,1),
        a = abs(rnorm(J,0,1)),
        b = rnorm(J,0,1),
        #g = runif(J,0,0.3),
        g = runif(1,0,0.3),
        tau = abs(rnorm(J,0,1)),
        gamma = abs(rnorm(1,0,1)))
    }

    model <- nimbleModel(code=EMIRT , constants = const,
                    data = data, inits = list())
    
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE) 
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    nodes = list(dataNodes = dataNodes, parentNodes = parentNodes, simNodes=simNodes)
    
    mcmcfit <- nimbleMCMC(model = model,
                        monitors = c(parentNodes,params.to.monitor),
                        niter = mcmc.config$no.of.iterations,
                        nburnin = mcmc.config$no.of.burn.ins,
                        nchains = mcmc.config$no.of.chains,
                        samplesAsCodaMCMC = TRUE)
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes))
}

#===========================================================================
#                            Normative Method
#===========================================================================
Normative.fit <- function(dat, k=0.1, 
                    params.to.monitor=NULL,
                    mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 3
                                    )){
    X <- dat$X
    RT <- dat$RT
    RT[is.na(RT)] <- 0
    RT <- sapply(as.data.frame(RT), function(col) {
  col > k*mean(col, na.rm = TRUE)
})
    I <- nrow(X)
    J <- ncol(X)
    const <- list(I=I,J=J)
    data <- list(Y=X,class_p = RT)
    inits <- function(){
    list(theta=rnorm(I,0,1),
        a = abs(rnorm(J,0,1)),
        b = rnorm(J,0,1),
        #g = runif(J,0,0.3),
        g = runif(1,0,0.3),
        tau = abs(rnorm(J,0,1)),
        gamma = abs(rnorm(1,0,1)))
    }

    model <- nimbleModel(code=EMIRT , constants = const,
                    data = data, inits = list())
    
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE) 
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    nodes = list(dataNodes = dataNodes, parentNodes = parentNodes, simNodes=simNodes)
    
    mcmcfit <- nimbleMCMC(model = model,
                        monitors = c(parentNodes,params.to.monitor),
                        niter = mcmc.config$no.of.iterations,
                        nburnin = mcmc.config$no.of.burn.ins,
                        nchains = mcmc.config$no.of.chains,
                        samplesAsCodaMCMC = TRUE)
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes))
}
#===========================================================================
#                            scripts
#===========================================================================
source("scripts/data_handling_utils.r")
data <- load_and_clean_data(747119)
temp <- ILCRI.fit(data, 
        mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 1
                                    ))
MCMCsummary(temp$posterior_samples)


source("scripts/data_handling_utils.r")
data <- load_and_clean_data(747119)
temp <- Normative.fit(data, 
        mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 1
                                    ))


discretized_RT <- 
