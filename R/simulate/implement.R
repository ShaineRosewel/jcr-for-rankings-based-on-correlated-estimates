source("R/simulate/compute_metrics.R")
library(doParallel)

block_corr <- function(block_sizes, rho_within, rho_between) {
  B <- length(block_sizes)
  N <- sum(block_sizes)
  
  if (length(rho_within) == 1) {
    rho_within <- rep(rho_within, B)
  }
  
  Sigma <- matrix(rho_between, nrow = N, ncol = N)
  diag(Sigma) <- 1
  
  idx <- 1
  for (b in seq_len(B)) {
    block_idx <- idx:(idx + block_sizes[b] - 1)
    Sigma[block_idx, block_idx] <- rho_within[b]
    diag(Sigma[block_idx, block_idx]) <- 1
    idx <- idx + block_sizes[b]
  }
  
  return(Sigma)
}


cl=parallel::makeCluster(15)
registerDoParallel(cl)
case_start <- Sys.time()


df <- readRDS("mean_travel_time_ranking_2011.rds")
mean <- 23.8
sds <- c(2) #c(3.6, 6.0) # disparity
Ks <- c(10) #c(20, 30, 40, 50)

filename <- "balanced-2-testing-reprod-from-mod"

for (sd in sds) {
  for (K in Ks) {
    cat("RESULTS FOR SD =", sd, "; K =", K, "\n")
    
    corr_matrix <- block_corr(block_sizes=c(.5*K,.5*K), #2B CASE
                              rho_within=c(0.9, 0.1), 
                              rho_between=0.0)
    
    set.seed(123974)
    true_theta <- rnorm(K, mean, sd)
    se <- df$S[1:K]
    variance_vector <- se^2
    delta <- diag(variance_vector)
    varcovar_matrix <- delta^(1/2) %*% corr_matrix %*% delta^(1/2)
    out <- implement_algorithm(
      true_theta = true_theta,
      K = K, 
      reps=5000,
      B = 500,
      alpha = 0.05,
      C = 300,
      varcovar_matrix = varcovar_matrix)
    case <- paste(paste0("K",gsub("\\.", "p", as.character(K))),
                  paste0("sd",gsub("\\.", "p", as.character(sd))),
                  paste0("r",gsub("\\.", "p", filename)),
                  sep = "_")
    
    case_end <- Sys.time()
    case_runtime <- as.numeric(difftime(case_end, case_start, units = "mins"))
    cat("Finished at:", format(case_end, "%I:%M %p"), "\n")
    cat("Runtime for", case, round(case_runtime, 2), "minutes\n")
    print(colMeans(out))
    
    filename <- paste0(
      case,
      ".rds")
    
    saveRDS(out, file = filename)
  }
}
stopCluster(cl)
