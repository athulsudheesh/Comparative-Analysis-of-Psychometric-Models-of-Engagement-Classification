# Function to simulate data from the ILCRI model
simulate_ILCRI <- function(I = 100, J = 20, mean_C=0.1, seed = 123) {
  set.seed(seed)
  
  # Create correlation matrix for person parameters
  rho_theta_eta <- 0.6  # Correlation between theta and eta
  rho_theta_xi <- 0.3   # Correlation between theta and xi
  rho_eta_xi <- 0.2     # Correlation between eta and xi
  
  Sigma_person <- matrix(c(
    1, rho_theta_eta, rho_theta_xi,
    rho_theta_eta, 1, rho_eta_xi,
    rho_theta_xi, rho_eta_xi, 1
  ), nrow = 3)
  
  # Generate item parameters
  b <- rnorm(J, 0, 1)        # Item difficulties
  kappa <- rnorm(J, 0, 1)    # Item parameters for class membership
  nu_til <- abs(rnorm(J, 0, 1))  # Baseline response time parameter
  delta <- abs(rnorm(J, 0, 1))   # Time difference parameter
  lambda <- abs(rnorm(J, 0, 1))   # Time discrimination parameter
  sigma2_epsilon <- rgamma(J, 1/2, 1/2)  # Error variances
  sigma2_epsilon_til <- rgamma(1, 1/2, 1/2)  # Error variance for non-class
  g <- runif(1, 0, 0.3)  # Guessing parameter
  
  # Calculate class membership probability
    mean = mean_C
    sd = 0.05
        
    shape1 = mean^2 * (1-mean) / sd^2 - mean
    shape2 = mean * (1 - mean)^2 / sd^2 + mean - 1
    class_p <- rbeta(I*J, shape1 = shape1, shape2 = shape2)
    class_p <- matrix(class_p,I,J)
  
    # 1. Calculate average class probability for each person
    person_avg_class_p <- rowMeans(class_p)
    # 2. Convert probabilities to eta values on a logit scale
    #    We need to invert the logistic function: eta = kappa + logit(class_p)/1.7
    #    Since we're working with average probabilities across items, we'll set kappa to zero for this calculation
    eta <- qlogis(person_avg_class_p) / 1.7
# 3. Generate correlated theta and xi values based on eta
  #    Using conditional distributions from multivariate normal
  
  # Create correlation matrix
  corr_matrix <- matrix(c(
    1, rho_theta_eta, rho_theta_xi,
    rho_theta_eta, 1, rho_eta_xi,
    rho_theta_xi, rho_eta_xi, 1
  ), nrow = 3)
  
  # Standard deviations for each parameter (all set to 1 in this case)
  sds <- c(1, 1, 1)
  
  # Convert correlation matrix to covariance matrix
  cov_matrix <- diag(sds) %*% corr_matrix %*% diag(sds)
  
  # Extract conditional means and variances for theta and xi given eta
  cond_mean_theta <- eta * cov_matrix[1, 2] / cov_matrix[2, 2]
  cond_var_theta <- cov_matrix[1, 1] - (cov_matrix[1, 2]^2 / cov_matrix[2, 2])
  
  cond_mean_xi <- eta * cov_matrix[3, 2] / cov_matrix[2, 2]
  cond_var_xi <- cov_matrix[3, 3] - (cov_matrix[3, 2]^2 / cov_matrix[2, 2])
  
  # Generate theta and xi conditionally on eta
  theta <- rnorm(I, cond_mean_theta, sqrt(cond_var_theta))
  xi <- rnorm(I, cond_mean_xi, sqrt(cond_var_xi))
  # Generate person parameters
  #person_par <- MASS::mvrnorm(I, mu = rep(0, 3), Sigma = Sigma_person)
  #theta <- person_par[, 1]  # Ability
  #eta <- person_par[, 2]    # Class membership trait
  #xi <- person_par[, 3]     # Speed
  
  # Initialize output matrices
  Y <- matrix(NA, I, J)     # Binary responses
  C <- matrix(NA, I, J)     # Class memberships
  T_mat <- matrix(NA, I, J) # Response times
  
  # Generate data
  for (i in 1:I) {
    for (j in 1:J) {

      # Generate class membership indicator
      C[i, j] <- rbinom(1, 1, class_p)
      
      # Calculate item response probability
      irp <- plogis(1.7 * (theta[i] - b[j]))
      correct_resp_prob <- C[i, j] * irp + (1 - C[i, j]) * g
      
      # Generate binary response
      Y[i, j] <- rbinom(1, 1, correct_resp_prob)
      
      # Calculate response time parameters
      mu_l_ij <- nu_til[j] + lambda[j] * xi[i] + C[i, j] * delta[j]
      sigma2_time <- C[i, j] * sigma2_epsilon[j] + (1 - C[i, j]) * sigma2_epsilon_til
      
      # Generate response time
      T_mat[i, j] <- rnorm(1, mu_l_ij, sqrt(sigma2_time))
    }
  }
  
  # Return simulated data and parameters
  return(list(
    Y = Y,               # Binary responses
    C = C,               # Class memberships
    T_mat = T_mat,       # Response times
    theta = theta,       # Person abilities
    eta = eta,           # Person class membership traits 
    xi = xi,             # Person speed parameters
    b = b,               # Item difficulties
    kappa = kappa,       # Item class parameters
    nu_til = nu_til,     # Baseline response time
    delta = delta,       # Time difference parameter
    lambda = lambda,     # Time discrimination
    sigma2_epsilon = sigma2_epsilon,  # Time variances
    sigma2_epsilon_til = sigma2_epsilon_til,  # Time variance non-class
    g = g,               # Guessing parameter
    Sigma_person = Sigma_person  # Person parameter correlation matrix
  ))
}

library(ggplot2)
library(ggpubr)
plot(density(log(sim_data$T_mat)))

# Usage example
sim_data <- simulate_ILCRI(I = 60, J = 15, mean_C = 0.9)
saveRDS(sim_data, "analyses/simulations/C_90/sim_data.rds")

sim_data <- simulate_ILCRI(I = 60, J = 15, mean_C = 0.6)
saveRDS(sim_data, "analyses/simulations/C_60/sim_data.rds")

sim_data <- simulate_ILCRI(I = 60, J = 15, mean_C = 0.3)
saveRDS(sim_data, "analyses/simulations/C_30/sim_data.rds")


sim_data <- simulate_ILCRI(I = 60, J = 15, mean_C = 0.1)
saveRDS(sim_data, "analyses/simulations/C_10/sim_data.rds")
# Example: Examining the first few rows of the simulated response matrix
plot(density(exp(sim_data$T_mat)),xlim=c(0,50))




# Function to simulate data from the ILCRI model with manual class probability control
simulate_ILCRI <- function(I = 100, J = 20, mean_C=0.1, seed = 145) {
  set.seed(seed)
  
  # Create correlation matrix for person parameters
  rho_theta_eta <- 0.6  # Correlation between theta and eta
  rho_theta_xi <- 0.3   # Correlation between theta and xi
  rho_eta_xi <- 0.2     # Correlation between eta and xi
  
  Sigma_person <- matrix(c(
    1, rho_theta_eta, rho_theta_xi,
    rho_theta_eta, 1, rho_eta_xi,
    rho_theta_xi, rho_eta_xi, 1
  ), nrow = 3)
  
  # Generate item parameters
  b <- rnorm(J, 0, 1)        # Item difficulties
  kappa <- rnorm(J, 0, 1)    # Item parameters for class membership
  nu_til <- abs(rnorm(J, 0, 1))  # Baseline response time parameter
  delta <- abs(rnorm(J, 0, 1))   # Time difference parameter
  lambda <- abs(rnorm(J, 0, 1))   # Time discrimination parameter
  sigma2_epsilon <- rgamma(J, 1/2, 1/2)  # Error variances
  sigma2_epsilon_til <- rgamma(1, 1/2, 1/2)  # Error variance for non-class
  g <- runif(1, 0, 0.3)  # Guessing parameter
  
  # Generate person parameters
  person_par <- MASS::mvrnorm(I, mu = rep(0, 3), Sigma = Sigma_person)
  theta <- person_par[, 1]  # Ability
  eta <- person_par[, 2]    # Class membership trait
  xi <- person_par[, 3]     # Speed
  
  # Initialize output matrices
  Y <- matrix(NA, I, J)     # Binary responses
  C <- matrix(NA, I, J)     # Class memberships
  T_mat <- matrix(NA, I, J) # Response times
      mean = mean_C
    sd = 0.05
        
    shape1 = mean^2 * (1-mean) / sd^2 - mean
    shape2 = mean * (1 - mean)^2 / sd^2 + mean - 1
  class_p <- rbeta(I*J, shape1 = shape1, shape2 = shape2)
  class_p_matrix <- matrix(class_p,I,J) # Store class probabilities
  
  
  # Generate data
  for (i in 1:I) {
    for (j in 1:J) {
      # Get class membership probability
      class_p <- class_p_matrix[i, j]
      
      # Generate class membership indicator
      C[i, j] <- rbinom(1, 1, class_p)
      
      # Calculate item response probability
      irp <- plogis(1.7 * (theta[i] - b[j]))
      correct_resp_prob <- C[i, j] * irp + (1 - C[i, j]) * g
      
      # Generate binary response
      Y[i, j] <- rbinom(1, 1, correct_resp_prob)
      
      # Calculate response time parameters
      mu_l_ij <- nu_til[j] + lambda[j] * xi[i] + C[i, j] * delta[j]
      sigma2_time <- C[i, j] * sigma2_epsilon[j] + (1 - C[i, j]) * sigma2_epsilon_til
      
      # Generate response time
      T_mat[i, j] <- rnorm(1, mu_l_ij, sqrt(sigma2_time))
    }
  }
  
  # Return simulated data and parameters
  return(list(
    Y = Y,               # Binary responses
    C = C,               # Class memberships
    T_mat = T_mat,       # Response times
    class_p = class_p_matrix, # Class probabilities used
    theta = theta,       # Person abilities
    eta = eta,           # Person class membership traits 
    xi = xi,             # Person speed parameters
    b = b,               # Item difficulties
    kappa = kappa,       # Item class parameters
    nu_til = nu_til,     # Baseline response time
    delta = delta,       # Time difference parameter
    lambda = lambda,     # Time discrimination
    sigma2_epsilon = sigma2_epsilon,  # Time variances
    sigma2_epsilon_til = sigma2_epsilon_til,  # Time variance non-class
    g = g,               # Guessing parameter
    Sigma_person = Sigma_person  # Person parameter correlation matrix
  ))
}

# Usage example
sim_data <- simulate_ILCRI(I = 60, J = 15, mean_C = 0.9)
saveRDS(sim_data, "analyses/simulations/C_90/sim_data.rds")

sim_data <- simulate_ILCRI(I = 60, J = 15, mean_C = 0.6)
saveRDS(sim_data, "analyses/simulations/C_60/sim_data.rds")

sim_data <- simulate_ILCRI(I = 60, J = 15, mean_C = 0.3)
saveRDS(sim_data, "analyses/simulations/C_30/sim_data.rds")


sim_data <- simulate_ILCRI(I = 60, J = 15, mean_C = 0.1)
saveRDS(sim_data, "analyses/simulations/C_10/sim_data.rds")
# Example: Examining the first few rows of the simulated response matrix