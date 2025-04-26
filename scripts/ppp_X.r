library(progress)
files.sources = list.files("models", full.names=TRUE)
sapply(files.sources, source)
source("scripts/threshold_functions.r")
library(progress)
ppp_commonk <- function(dat,samples, Nsim =200){

    k <- 0.5
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

    res_matrix <- as.matrix(samples)
    nSamples <- nrow(res_matrix)
    postNames <- colnames(res_matrix)
    step_width <- as.integer(nSamples/Nsim)
    res_sub <- res_matrix[seq(from =1, to=nSamples,by=step_width),]
    ppScoreDist <- matrix(0,nrow(res_sub),nrow(X))
    ppDDist <- matrix(0,nrow(res_sub),ncol(X))
    ppDrep <- matrix(0,nrow(res_sub),1)
    ppDtheor <- matrix(0, nrow(res_sub),1)
    set.seed(100)
    X2 = X
    X2[is.na(X2)] <- 0
    pb <- progress_bar$new(total =100)
    for(i in 1:nrow(res_sub)){
    pb$tick()
    X_pred <- matrix(0,nrow(X),ncol(X))
   
   values(model, postNames) <- res_sub[i,]
   model$simulate(simNodes, includeData = TRUE)
   
   X_pred <- model[["Y"]]
   pi <- model[["correct_resp_prob"]]
   D_rep <- mean((X_pred - pi)^2/(pi*(1-pi)))
   D_theoretical <- mean((as.matrix(X2) - pi)^2/(pi*(1-pi)))

   ppScoreDist[i,] <- apply(X_pred,1,sum)
   ppDDist[i,] <- apply(X_pred,2,mean)
   ppDrep[i,] <- D_rep
   ppDtheor[i,] <- D_theoretical
}
    mean(ppDrep >= ppDtheor)
}


ppp_normative <- function(dat,samples, Nsim =200){

    k <- 0.5
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

    res_matrix <- as.matrix(samples)
    nSamples <- nrow(res_matrix)
    postNames <- colnames(res_matrix)
    step_width <- as.integer(nSamples/Nsim)
    res_sub <- res_matrix[seq(from =1, to=nSamples,by=step_width),]
    ppScoreDist <- matrix(0,nrow(res_sub),nrow(X))
    ppDDist <- matrix(0,nrow(res_sub),ncol(X))
    ppDrep <- matrix(0,nrow(res_sub),1)
    ppDtheor <- matrix(0, nrow(res_sub),1)
    set.seed(100)
    X2 = X
    X2[is.na(X2)] <- 0
    pb <- progress_bar$new(total =200)
    for(i in 1:nrow(res_sub)){
    pb$tick()
    X_pred <- matrix(0,nrow(X),ncol(X))
   
   values(model, postNames) <- res_sub[i,]
   model$simulate(simNodes, includeData = TRUE)
   
   X_pred <- model[["Y"]]
   pi <- model[["correct_resp_prob"]]
   D_rep <- mean((X_pred - pi)^2/(pi*(1-pi)))
   D_theoretical <- mean((as.matrix(X2) - pi)^2/(pi*(1-pi)))

   ppScoreDist[i,] <- apply(X_pred,1,sum)
   ppDDist[i,] <- apply(X_pred,2,mean)
   ppDrep[i,] <- D_rep
   ppDtheor[i,] <- D_theoretical
}
    mean(ppDrep >= ppDtheor)
}



ppp_vics <- function(dat,samples, Nsim =200){

    k <- 0.5
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

    res_matrix <- as.matrix(samples)
    nSamples <- nrow(res_matrix)
    postNames <- colnames(res_matrix)
    step_width <- as.integer(nSamples/Nsim)
    res_sub <- res_matrix[seq(from =1, to=nSamples,by=step_width),]
    ppScoreDist <- matrix(0,nrow(res_sub),nrow(X))
    ppDDist <- matrix(0,nrow(res_sub),ncol(X))
    ppDrep <- matrix(0,nrow(res_sub),1)
    ppDtheor <- matrix(0, nrow(res_sub),1)
    set.seed(100)
    X2 = X
    X2[is.na(X2)] <- 0
    pb <- progress_bar$new(total =200)
    for(i in 1:nrow(res_sub)){
    pb$tick()
    X_pred <- matrix(0,nrow(X),ncol(X))
   
   values(model, postNames) <- res_sub[i,]
   model$simulate(simNodes, includeData = TRUE)
   
   X_pred <- model[["Y"]]
   pi <- model[["correct_resp_prob"]]
   D_rep <- mean((X_pred - pi)^2/(pi*(1-pi)))
   D_theoretical <- mean((as.matrix(X2) - pi)^2/(pi*(1-pi)))

   ppScoreDist[i,] <- apply(X_pred,1,sum)
   ppDDist[i,] <- apply(X_pred,2,mean)
   ppDrep[i,] <- D_rep
   ppDtheor[i,] <- D_theoretical
}
    mean(ppDrep >= ppDtheor)
}

ppp_vii <- function(dat,samples, Nsim =200){

    k <- 0.5
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

    res_matrix <- as.matrix(samples)
    nSamples <- nrow(res_matrix)
    postNames <- colnames(res_matrix)
    step_width <- as.integer(nSamples/Nsim)
    res_sub <- res_matrix[seq(from =1, to=nSamples,by=step_width),]
    ppScoreDist <- matrix(0,nrow(res_sub),nrow(X))
    ppDDist <- matrix(0,nrow(res_sub),ncol(X))
    ppDrep <- matrix(0,nrow(res_sub),1)
    ppDtheor <- matrix(0, nrow(res_sub),1)
    set.seed(100)
    X2 = X
    X2[is.na(X2)] <- 0
    pb <- progress_bar$new(total =200)
    for(i in 1:nrow(res_sub)){
    pb$tick()
    X_pred <- matrix(0,nrow(X),ncol(X))
   
   values(model, postNames) <- res_sub[i,]
   model$simulate(simNodes, includeData = TRUE)
   
   X_pred <- model[["Y"]]
   pi <- model[["correct_resp_prob"]]
   D_rep <- mean((X_pred - pi)^2/(pi*(1-pi)))
   D_theoretical <- mean((as.matrix(X2) - pi)^2/(pi*(1-pi)))

   ppScoreDist[i,] <- apply(X_pred,1,sum)
   ppDDist[i,] <- apply(X_pred,2,mean)
   ppDrep[i,] <- D_rep
   ppDtheor[i,] <- D_theoretical
}
    mean(ppDrep >= ppDtheor)
}

ppp_gomrt <- function(dat,samples, Nsim =200){

    k <- 0.5
    X <- dat$X
    RT <- dat$RT
    RT[is.na(RT)] <- 0
    
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

    model <- nimbleModel(code=GoMRT, constants = const,
                    data = data, inits = list())
    
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE) 
    simNodes <- model$getDependencies(parentNodes, self = FALSE)

    res_matrix <- as.matrix(samples)
    nSamples <- nrow(res_matrix)
    postNames <- colnames(res_matrix)
    step_width <- as.integer(nSamples/Nsim)
    res_sub <- res_matrix[seq(from =1, to=nSamples,by=step_width),]
    ppScoreDist <- matrix(0,nrow(res_sub),nrow(X))
    ppDDist <- matrix(0,nrow(res_sub),ncol(X))
    ppDrep <- matrix(0,nrow(res_sub),1)
    ppDtheor <- matrix(0, nrow(res_sub),1)
    set.seed(100)
    X2 = X
    X2[is.na(X2)] <- 0
    pb <- progress_bar$new(total =200)
    for(i in 1:nrow(res_sub)){
    pb$tick()
    X_pred <- matrix(0,nrow(X),ncol(X))
   
   values(model, postNames) <- res_sub[i,]
   model$simulate(simNodes, includeData = TRUE)
   
   X_pred <- model[["Y"]]
   pi <- model[["correct_resp_prob"]]
   D_rep <- mean((X_pred - pi)^2/(pi*(1-pi)))
   D_theoretical <- mean((as.matrix(X2) - pi)^2/(pi*(1-pi)))

   ppScoreDist[i,] <- apply(X_pred,1,sum)
   ppDDist[i,] <- apply(X_pred,2,mean)
   ppDrep[i,] <- D_rep
   ppDtheor[i,] <- D_theoretical
}
    mean(ppDrep >= ppDtheor)
}
