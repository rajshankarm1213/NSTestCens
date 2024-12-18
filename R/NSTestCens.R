#
# This script contains the code to check for the goodness-of-fit of a non-parametric hazard function estimate against any 
# given parametric model or formula of the hazard function using the Neyman Smooth Test.
#

# Load the necessary libraries
library(survival); library(pracma); # Used for Legendre Polynomials


# Simulate survival data
sim.survdata <- function(n, lambda, censor_rate){
    set.seed(123)
    T <- rexp(n, rate= lambda)
    C <- rexp(n, rate = censor_rate)
    Y <- pmin(T, C)
    delta <- as.numeric(T <= C)
    return(data.frame(time = Y, status = delta))
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

Monomial.Basis <- function(U, degree){
    # Function to calculate the basis function as a monomial
    # Args: 
    # U: A vector of generalized residuals
    # degree: The degree of the monomial basis
    # Returns: A matrix of monomial basis functions

    return(sapply(1:degree, function(l) U^(l)))
}

LegendrePoly <- function(U, degree){
    # Function to compute the Legendre Polynomials
    # Args:
    # U: A vector of generalized residuals
    # degree: The degree of the Legendre Polynomial
    # Returns: A matrix of Legendre Polynomials of the given degree

    # Transform U to the [-1, 1] interval
    V <- 2 * U - 1

    # Initialize a matrix to store Legendre polynomials
    legendre_matrix <- matrix(NA, nrow = length(U), ncol = degree)

    # Compute Legendre polynomials for each degree
    for (d in 1:degree) {
        L <- legendre(d, V)  # Returns (d+1) x length(U) matrix
        legendre_matrix[, d] <- L[1, ]  # Extract the first row (Legendre polynomial)
    }

    return(legendre_matrix)
}

Adj.BasisFunction <- function(U, delta, degree, basis_fn){
    # Function to adjust the Legendre Polynomial Score for censoring
    # Args:
    # U: A vector of generalized residuals
    # delta: A vector of censoring indicators
    # degree: The degree of the Legendre Polynomial
    # basis_fn: The basis function to use. 
    # Returns: A matrix of adjusted Legendre Polynomial Scores

    basis_matrix <- basis_fn(U, degree)
    adj_scores <- basis_matrix
    # For censored observations, we integrate the basis functions over the interval [U, 1]
    for(i in 1:length(U)){
        if(delta[i] == 0){ # Censored
            for(l in 1:degree){
                psi_l <- function(x) basis_fn(x, degree)[,l]
                adj_scores[i, l] <- integrate(psi_l, lower = U[i], upper = 1)$value / (1 - U[i])
            }
        }
    }
    return(adj_scores)
}

Compute.QSS <- function(adj_scores, n){
    # Function to compute the Quadratic Score Statistic
    # Args:
    # adj_score: A matrix of adjusted Legendre Polynomial Scores
    # n: The number of observations
    # Returns: The Quadratic Score Statistic

    theta_k <- colMeans(adj_scores)

    # Compute the test statistic
    Q_m <- n * sum(theta_k^2)
    return(Q_m)
}

NSTest.Cens <- function(data, degree=NA, basis_type = "m"){
    # Function to perform the Neyman Smooth Test to check for goodness-of-fit 
    # Args:
    # data: A data frame containing the survival data, in the form of $time, $status
    # degree: The degree of the Legendre Polynomial
    # basis_type: The type of basis function to use. Accepted values are 'l' for Legendre Polynomials and 'm' for monomial basis. Default is 'm'.

    
    if(!("time" %in% colnames(data)) || !("status" %in% colnames(data))){
        stop("The data frame must contain columns named 'time' and 'status'")
    }
    # If degree is NA and basis_type is 'l', then stop and ask for the degree
    if(is.na(degree)){
        stop("The degree of the Legendre Polynomial must be specified")
    }
    if(degree < 1){
        stop("The degree of the Orthogonal Polynomial must be greater than 0")
    }
    if(basis_type != "l" && basis_type != "m"){
        stop("The basis type must be either 'l' for Legendre Polynomials or 'm' for monomial basis")
    }

    basis_fn <- switch(basis_type, "l" = LegendrePoly, "m" = Monomial.Basis)
    residuals <- Compute.GeneralizedResiduals(data)

    adj_scores <- Adj.BasisFunction(residuals$U, residuals$status, degree, basis_fn)

    n <- nrow(data)
    Q_m <- Compute.QSS(adj_scores, n)

    # Compute the p-value
    p_val <- pchisq(Q_m, df = degree)

    # Return results
    return(list(Q_m = Q_m, p_val = p_val, basis_type = basis_type, residuals = residuals, adj_scores = adj_scores))
}
