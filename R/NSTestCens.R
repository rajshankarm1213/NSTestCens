#
# This script contains the code to check for the goodness-of-fit of a non-parametric hazard function estimate against any 
# given parametric model or formula of the hazard function using the Neyman Smooth Test.
#

# Load the necessary libraries
library(survival); library(pracma) # Used for Legendre Polynomials


# Simulate survival data
sim.survdata <- function(n, lambda, censor_rate){
    set.seed(123)
    T <- rexp(n, rate= lambda)
    C <- rexp(n, rate = censor_rate)
    Y <- pmin(T, C)
    delta <- as.numeric(T <= C)
    return(data.frame(Y = time, delta = status))
}

Compute.GeneralizedResiduals <- function(data){
    # Function to compute the generalized residuals
    # Args: 
    # data: A data frame containing the survival data, in the form of $time, $status
    # Returns: A data frame consisting the of generalized residuals as 'U' and the censoring indicator as 'status'

    if(!("time" %in% colnames(data)) || !("status" %in% colnames(data))){
        stop("The data frame must contain columns named 'time' and 'status'")
    }

    # Fit a Kaplan-Meier estimator to the given data, to calculate a non-parametric estimate of the cumulative hazard function
    km_fit <- survfit(Surv(time, status) ~ 1, data = data)
    surv_probs <- km_fit$surv
    surv_times <- km_fit$time


    H0 <- -log(surv_probs)

    h_obs <- approx(x = surv_times, y = H0, xout = data$time, method='constant')$y
    
    U <- 1 - exp(-h_obs) # Generalized residuals

    return(data.frame(U = U, status = data$status))
}

LegendrePoly <- function(U, degree){
    # Function to compute the Legendre Polynomials
    # Args:
    # U: A vector of generalized residuals
    # degree: The degree of the Legendre Polynomial
    # Returns: A matrix of Legendre Polynomials of the given degree
    polynomials <- sapply(1:degree, function(k) legendre(k, 2 * U - 1))
    colnames(polynomials) <- paste0("phi_", 1:degree)
    return(polynomials)
}

Adj.LegendrePoly <- function(U, delta, degree){
    # Function to adjust the Legendre Polynomial Score for censoring
    # Args:
    # U: A vector of generalized residuals
    # delta: A vector of censoring indicators
    # degree: The degree of the Legendre Polynomial
    polynomials <- LegendrePoly(U, degree)
    adj_scores <- polynomials
    # For censored observations, we integrate the basis functions over the interval [U, 1]
    for(i in 1:length(U)){
        if(delta[i] == 0){ # Censored
            phi_k <- function(x) legendre(k, 2 * x - 1)
            adj_scores[i, k] <- integrate(phi_k, lower = U[i], upper = 1)$value / (1 - U[i])
        }
    }
    return(adj_scores)
}

Compute.QSS <- function(adj_scores, n){
    # Function to compute the Quadratic Score Statistic
    # Args:
    # adj_score: A matrix of adjusted Legendre Polynomial Scores
    # n: The number of observations

    theta_k <- colMeans(adj_scores)

    # Compute the test statistic
    Q_m <- n * sum(theta_k^2)
    return(Q_m)
}

NSTest.Cens <- function(data, lambda_null, degree, basis_type = "m"){
    # Function to perform the Neyman Smooth Test to check for goodness-of-fit 
    # Args:
    # data: A data frame containing the survival data, in the form of $time, $status
    # lambda_null: A function representing the null hypothesis of the hazard function
    # degree: The degree of the Legendre Polynomial
    # basis_type: The type of basis function to use. Accepted values are 'l' for Legendre Polynomials and 'm' for monomial basis. Default is 'm'.

    # Sanity checks
    if(!("time" %in% colnames(data)) || !("status" %in% colnames(data))){
        stop("The data frame must contain columns named 'time' and 'status'")
    }
    if(degree < 1){
        stop("The degree of the Legendre Polynomial must be greater than 0")
    }
    if(basis_type != "l" && basis_type != "m"){
        stop("The basis type must be either 'l' for Legendre Polynomials or 'm' for monomial basis")
    }

    residuals <- Compute.GeneralizedResiduals(data, lambda_null)

    adj_scores <- Adj.LegendrePoly(residuals$U, residuals$status, degree)

    n <- nrow(data)
    Q_m <- Compute.QSS(adj_scores, n)

    # Compute the p-value
    p_val <- 1 - pchisq(Q_m, df = degree)

    # Print results
    list(Q_m = Q_m, p_val = p_val, basis_type = basis_type, residuals = residuals, adj_scores = adj_scores)
}