files.sources = list.files("models", full.names=TRUE)
sapply(files.sources, source)
source("scripts/threshold_functions.r")
#===================================================================================================
#                                          Latent Variable Approaches
#===================================================================================================
#===========================================================================
#                            GoMRT
#===========================================================================
GoMRT.fit <- function(dat, 
                    params.to.monitor=NULL, cross_validate = TRUE,
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

    model <- nimbleModel(code=GoMRT, constants = const,
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
                        samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, cvvals=cvvals))
}

#============================= END OF GoMRT ==============================

#===========================================================================
#                            DLC-SL
#===========================================================================
DLCSL.fit <- function(dat, 
                    params.to.monitor=NULL, cross_validate = TRUE,
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
                        samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, cvvals=cvvals))
}

#============================= END OF DLC-SL ==============================


#===========================================================================
#                            DLC-TL
#===========================================================================
DLCTL.fit <- function(dat, 
                    params.to.monitor=NULL, cross_validate = TRUE,
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
                            samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, cvvals=cvvals))
}

#============================= END OF DLC-TL ==============================


#===========================================================================
#                           MHM
#===========================================================================

MHM.fit <- function(dat, 
                    params.to.monitor=NULL, cross_validate = TRUE, 
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
    Sigma <- matrix(c(1, 0.4, 
                  0.4, 1), 
               byrow = TRUE, ncol = 2)
    inits <- function(){
    list(person_par=MASS::mvrnorm(I,mu=rep(0,2),Sigma = Sigma,empirical = T),
        a = abs(rnorm(J,0,1)),
        b = rnorm(J,0,1),
        #g = runif(J,0,0.3),
        g = runif(1,0,0.3),
        kappa = rnorm(J,0,1),
        nu_til = -abs(rnorm(1,0,1)),
        delta = abs(rnorm(J,1,1)),
        lambda = abs(rnorm(1,0.5,0.5)),
        sigma2_epsilon = abs(rnorm(J,0,1)),
        sigma2_epsilon_til = abs(rnorm(1,0,1)),
        rho_theta_eta = runif(1,0.5,1),
        rho_theta_xi = runif(1,-0.5,0.5),
        rho_eta_xi = runif(1,0.5,1),
        Sigma_person = diag(2)
    )
    }
    model <- nimbleModel(code=MHM, constants = const,
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
                            samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, cvvals=cvvals))                            
}
#============================= END OF MHM ==============================
#===========================================================================
#                            ILC-RE
#===========================================================================

ILCRE.fit <- function(dat, 
                    params.to.monitor=NULL, cross_validate = TRUE, 
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
                            samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, cvvals=cvvals))                            
}
#============================= END OF ILC-RE ==============================



#===========================================================================
#                            ILC-RI
#===========================================================================

ILCRI.fit <- function(dat, 
                    params.to.monitor=NULL, cross_validate = TRUE, 
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
                            samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, cvvals=cvvals)) 
}

#============================= END OF ILC-RI ==============================

#===================================================================================================
#                                          Threshold based Approaches
#===================================================================================================

#===========================================================================
#                            Common-K
#===========================================================================
CommonK.fit <- function(dat, k=5, 
                    params.to.monitor= NULL, cross_validate = TRUE, 
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
                        samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, cvvals=cvvals))
}

#===========================================================================
#                            Normative Method
#===========================================================================
Normative.fit <- function(dat, k=0.1, 
                    params.to.monitor=NULL, cross_validate = TRUE, 
                    mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 3
                                    )){
    X <- dat$X
    RT <- data$RT
    RT[is.na(RT)] <- 0
    k = 0.1
    thresholds = k* colMeans(RT <- RT, na.rm = TRUE)
    RT<- t(t(RT) > thresholds)
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
                        samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, 
              thresholds = thresholds, cvvals=cvvals))
}

#===========================================================================
#                            VICS
#===========================================================================
VICS.fit <- function(dat,
                    params.to.monitor=NULL, cross_validate = TRUE, 
                    mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 3
                                    )){
    X <- dat$X
    RT <- dat$RT
    RT[is.na(RT)] <- 0
    X[is.na(X)] <- 0
    thresholds <- vics_method(as.matrix(X),as.matrix(RT))
    RT<- t(t(RT) > thresholds)
    
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
                        samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, thresholds=thresholds, cvvals=cvvals))
}
#===========================================================================
#                            VII
#===========================================================================
VII.fit <- function(dat,
                    params.to.monitor=NULL, cross_validate = TRUE, 
                    mcmc.config = list(
                                    no.of.iterations = 12000,
                                    no.of.burn.ins = 2000,
                                    no.of.chains = 3
                                    )){
    X <- dat$X
    RT <- dat$RT
    RT[is.na(RT)] <- 0
    X[is.na(X)] <- 0
    thresholds <- vii_thresholds(as.matrix(X),as.matrix(RT))
    RT<- t(t(RT) > thresholds)
    
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
                        samplesAsCodaMCMC = TRUE,WAIC = TRUE)
    cvvals = NULL
    if(cross_validate){
    print("K-fold Cross Validation started...")
    mcmcconfig <- configureMCMC(model = model)
    cvvals<- runCrossValidate(MCMCconfiguration = mcmcconfig, 
                              MCMCcontrol = list(niter = mcmc.config$no.of.iterations),
                                   k = 3)
    }
    return(list(posterior_samples = mcmcfit, model = model, nodes = nodes, thresholds=thresholds, cvvals=cvvals))
}

#===========================================================================
#                            scripts
#===========================================================================
