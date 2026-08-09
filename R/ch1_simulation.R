# 1. CLEAN UP
try(future::plan(sequential), silent = TRUE)
closeAllConnections()
gc() 

# 2. Load core packages
library(scoringRules)
library(MASS)
library(future.apply)
library(ggplot2)

# 3. Parameters
n_sims <- 800; n <- 2000; alpha <- 0.1
epsilon <- delta <- seq(-0.6, 0.6, length.out = 7)
params  <- expand.grid(epsilon = epsilon, delta = delta)
m        <- nrow(params)
eta_base <- 1.0

# 4. Pre-calculate c-values
c_mat <- matrix(NA, m, m); high <- 1e10
for (i in 1:m){
  for (j in (1:m)[-i]){
    tmp <- if (params$delta[i] != params$delta[j]) {
      ((params$epsilon[i])*sqrt(1+params$delta[j])-(params$epsilon[j])*sqrt(1+params$delta[i]))/(sqrt(1+params$delta[j])-sqrt(1+params$delta[i]))
    } else { numeric(0) }
    c_mat[i,j] <- max(abs(crps_norm(c(high,-high,tmp),params$epsilon[i],sqrt(1+params$delta[i])) - 
                            crps_norm(c(high,-high,tmp),params$epsilon[j],sqrt(1+params$delta[j]))))
  }
}
lam_fix <- 1 / (2 * c_mat)

# 5. Parallel Setup
plan(multisession, workers = parallel::detectCores() - 1)

# 6. Worker Function
run_sim_fast <- function(sim_id) {
  set.seed(sim_id)
  y <- rep(NA, n); y[1] <- rnorm(1)
  for (i in 2:n){ y[i] <- rnorm(1, mean = y[i - 1]) }
  means <- c(0, y[1 : (n - 1)])
  L <- matrix(NA, n, m)
  for (i in 1:m){ L[, i] <- crps_norm(y, mean = means + params$epsilon[i], sd = sqrt(1 + params$delta[i])) }
  d <- array(0, c(m, m, n))
  for (i in 1:m) { for (j in (1:m)[-i]) { d[i, j, ] <- L[, i] - L[, j] } }
  
  E_std <- array(1, c(m, m, n)); E_kel <- array(1, c(m, m, n))
  cd <- matrix(0, m, m); cv <- matrix(0, m, m)
  for (t in 1:n) {
    for (i in 1:m) {
      for (j in (1:m)[-i]) {
        E_std[i,j,t] <- (if(t==1) 1 else E_std[i,j,t-1]) * (1 + lam_fix[i,j] * d[i,j,t])
        l_t <- if(t==1) lam_fix[i,j] else max(0, min(cd[i,j]/(cv[i,j]+1e-9), 0.9/c_mat[i,j]))
        E_kel[i,j,t] <- (if(t==1) 1 else E_kel[i,j,t-1]) * (1 + l_t * d[i,j,t])
        cd[i,j] <- cd[i,j] + d[i,j,t]; cv[i,j] <- cv[i,j] + d[i,j,t]^2
      }
    }
  }
  
  EE <- list(Sym=matrix(0,m,n), FF=matrix(0,m,n), AF=matrix(0,m,n), FK=matrix(0,m,n), Full=matrix(0,m,n))
  cda <- matrix(0, m, m); cva <- matrix(0, m, m)
  for(t in 1:n){
    for(i in 1:m){
      others <- (1:m)[-i]
      EE$Sym[i,t] <- mean(E_kel[i,others,t])
      if(t>1){
        t_prev <- t - 1
        # Adaptive Eta: 1 / sum_{j != i} Var(d_i,j, 0:t-1)
        emp_vars = (cva[i, others] / t_prev) - (cda[i, others] / t_prev)^2
        eta_t = eta_base / (0.01 + sum(pmax(emp_vars, 0)))
        
        # Fixed Eta weights (e.g., eta=5)
        w_f <- exp(5*cda[i,others]-max(5*cda[i,others])); w_f <- w_f/sum(w_f)
        # Adaptive Eta weights
        w_a <- exp(eta_t*cda[i,others]-max(eta_t*cda[i,others])); w_a <- w_a/sum(w_a)
        
        EE$FF[i,t]=sum(w_f*E_std[i,others,t]); EE$AF[i,t]=sum(w_a*E_std[i,others,t])
        EE$FK[i,t]=sum(w_f*E_kel[i,others,t]); EE$Full[i,t]=sum(w_a*E_kel[i,others,t])
      } else { EE$FF[i,1]=EE$AF[i,1]=EE$FK[i,1]=EE$Full[i,1]=EE$Sym[i,1] }
    }
    cda=cda+d[,,t]; cva=cva+d[,,t]^2
  }
  
  res <- matrix(0, 5, n)
  for(k in 1:5){
    curr <- 1:m; dat <- EE[[k]]
    for(t in 1:n){
      elig <- which(dat[,t] < 1/alpha); curr <- elig[elig %in% curr]; res[k,t] <- length(curr)
    }
  }
  return(res)
}

# 7. Execute
cat("Running sims \n")
results <- future_lapply(1:n_sims, run_sim_fast, future.seed = TRUE)

# 8. Aggregate and Plot
avg_size <- Reduce("+", results) / n_sims
plot_df <- data.frame(t=rep(1:n, 5), Size=as.vector(t(avg_size)), 
                      Strategy=rep(c("Symmetric", "FixEta/FixLam", "AdaEta/FixLam", "FixEta/aGRAPALam", "AdaEta/aGRAPALam"), each=n))

ggplot(plot_df, aes(x=t, y=Size, color=Strategy)) + 
  geom_line(linewidth=1) + theme_bw() +
  scale_color_manual(values=c("gray", "red", "royalblue", "darkgreen", "purple"))

