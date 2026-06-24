get_ci_independent <- function(theta_hat,
                               S,
                               alpha){
  K <- length(theta_hat)
  gamma = 1-(1-alpha)^(1/K)
  z = qnorm(1-gamma/2)
  ci_lower <- theta_hat - z*S
  ci_upper <- theta_hat + z*S

  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}

get_ci_bonferroni <- function(theta_hat,
                              S,
                              alpha){
  K <- length(theta_hat)
  z = qnorm(1-(alpha/K)/2)
  ci_lower <- theta_hat - z*S
  ci_upper <- theta_hat + z*S

  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}

get_ci_rankbased_asymptotic <- function(B,
                                        theta_hat,
                                        varcovar_matrix,
                                        alpha){

  K <- length(theta_hat)
  
  thetahat_star <- MASS::mvrnorm(n = B,
                                 mu = theta_hat,
                                 Sigma = varcovar_matrix)
  
  sorted_thetahat_star <- t(apply(thetahat_star, 1, sort))

  variance_vector <- diag(varcovar_matrix)
  
  minuend <- thetahat_star^2 + matrix(variance_vector, B, K, byrow = TRUE)
  
  radicand <- t(apply(minuend, 1, sort)) - sorted_thetahat_star^2
  
  sigma_hat_star <- sqrt(
    t(apply(minuend, 1, sort)) - sorted_thetahat_star^2
    )

  sorted_theta_hat <- sort(theta_hat)

  compute_max <- function(b) {
    t_b <- max(abs(
      (sorted_thetahat_star[b, ] - sorted_theta_hat) /
        sigma_hat_star[b,]
    ))
    return(t_b)
  }

  t_star <- sapply(1:B, compute_max)

  t_hat <- quantile(t_star, probs = 1 - alpha)

  sigma_hat <- sqrt(
    sort(theta_hat^2 + variance_vector) - sorted_theta_hat^2
    )
  ci_lower <- sorted_theta_hat - t_hat*sigma_hat
  ci_upper <- sorted_theta_hat + t_hat*sigma_hat

  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}

get_ci_rankbased_level2bs <- function(B,
                                      C,
                                      theta_hat,
                                      varcovar_matrix,
                                      alpha) {

  K <- length(theta_hat)
  sorted_theta_hat <- sort(theta_hat)

  thetahat_star <- MASS::mvrnorm(n = B,
                                 mu = theta_hat,
                                 Sigma = varcovar_matrix) # B x K

  sorted_thetahat_star <- t(apply(thetahat_star, 1, sort)) # B x K

  generate_level2_data <- function(mu){MASS::mvrnorm(n = 1,
                                                     mu,
                                                     Sigma = varcovar_matrix)}

  thetahat_double_star <- # list of length B, each a K x C matrix
    apply(thetahat_star,
          1,
          function(thetahat_b) {replicate(C, # expected is a K x C matrix
                                          generate_level2_data(mu=thetahat_b))},
          simplify = FALSE
          )

  sorted_thetahat_double_star <- lapply(thetahat_double_star, 
                                        function(x)apply(x, 2, sort))

  # for each matrix b and for each row k
  # K x C matrix arg
  compute_sigma_hat <- function(mat, n){
    apply(mat, 
          1,
          function(row){
            sqrt(
              sum((row - mean(row))^2)/(n-1))
            })}
  
  # output must be B lists, each a vector of length K
  sigma_hat_star <- lapply(sorted_thetahat_double_star, # K by C 
                           function(x) compute_sigma_hat(x, C) 
                           )

  compute_max <- function(b) {
    t_b <- abs(
      (sorted_thetahat_star[b, ] - sorted_theta_hat) /
        sigma_hat_star[[b]]
    )
    max(t_b)
  }

  t_star <- sapply(1:B, compute_max)

  t_hat <- quantile(t_star, probs = 1 - alpha)

  sigma_hat <- compute_sigma_hat(t(sorted_thetahat_star),B)

  ci_lower <- sorted_theta_hat - t_hat*sigma_hat
  ci_upper <- sorted_theta_hat + t_hat*sigma_hat

  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}

get_ci_nonrankbased <- function(B, 
                                theta_hat,
                                alpha, 
                                varcovar_matrix) {

  K <- length(theta_hat)

  thetahat_star <- MASS::mvrnorm(n = B,
                                 mu = theta_hat,
                                 Sigma = varcovar_matrix) # B x K
  
  t_star <- apply(thetahat_star, 
                  1, 
                  function(x) max(
                    abs(
                    (x - theta_hat) / sqrt(diag(varcovar_matrix))
                    )
                    ))  
  
  t_hat <- quantile(t_star, probs = 1 - alpha)

  ci_lower <- theta_hat - t_hat*sqrt(diag(varcovar_matrix))
  ci_upper <- theta_hat + t_hat*sqrt(diag(varcovar_matrix))
  return(list(
    ci_lower = ci_lower,
    ci_upper = ci_upper
  ))
}