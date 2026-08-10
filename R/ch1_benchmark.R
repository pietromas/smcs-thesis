################################################################################
## ch1_benchmark.R -- data for Figures 1.2 and 1.3
##
## Benchmark against discrete argmin inference, five methods:
##
##   SMCS-w     : adaptively weighted SMCS, known penalty variance
##   SMCS-w-est : the same, but with the penalty variance ESTIMATED from the
##                data and corrected by the one-dimensional epoch construction
##                of Chapter 2. Per pair (i,k), the penalty variance is the
##                sample variance frozen on the doubling epoch grid
##                t_j = 8, 16, ..., inflated by the exact lower quantile of
##                chisq(nu)/nu at eps = delta/(K*(d-1)) per event, with
##                delta = 0.05, so each model's anytime guarantee is
##                alpha + delta. No oracle ingredient at all.
##   DA-plug, DA-adj : the two dimension-agnostic argmin tests
##   Bonferroni : one-sided t-test with Bonferroni correction
##
## Usage:   Rscript R/ch1_benchmark.R <scen: a|b|c> <rho> <null: 0|1> <nsims>
## Output:  results/bench_<scen><null>_rho<rho>.rds and a summary line
################################################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4)
  stop("Usage: Rscript R/ch1_benchmark.R <scen a|b|c> <rho> <null 0|1> <nsims>")
scen  <- args[1]; rho <- as.numeric(args[2])
isnull <- as.integer(args[3]) == 1L; ns <- as.integer(args[4])
if (!scen %in% c("a","b","c")) stop("scen must be a, b or c")

d      <- 100
nTot   <- 1000
alpha  <- 0.05
eta0   <- 1
t0     <- 10
grid   <- seq(100, nTot, by = 20)
kappa  <- 1e-8
delta  <- 0.05                                  # failure budget of the correction

## epoch grid for the estimated-variance correction
tks  <- 2^(3:9); tks <- tks[tks <= nTot]        # 8,16,...,512  (K = 7 epochs)
K    <- length(tks)
eps1 <- delta / (K * (d - 1))                   # budget per epoch and per pair
q1   <- qchisq(eps1, df = tks - 2) / (tks - 2)  # exact 1-D quantiles (no MC)

build_mu <- function(scen, d, isnull){
  if (scen == "a") mu <- c(0.1, 0, rep(0.1, d - 2))
  if (scen == "b") mu <- c(0.2, 0.1 + (seq(2, d) - 2) / (d - 2) * 0.9)
  if (scen == "c") mu <- c(0.05, 0, 0, 0, rep(10, d - 4))
  if (isnull) mu[1] <- min(mu)
  mu
}
mu    <- build_mu(scen, d, isnull)
Sigma <- rho^abs(outer(1:d, 1:d, "-")); L <- chol(Sigma)
Vd    <- outer(diag(Sigma), diag(Sigma), "+") - 2 * Sigma   # true variances
diag(Vd) <- 1

rowLSE <- function(A){ m <- apply(A, 1, max); m + log(rowSums(exp(A - m))) }

## Both SMCS variants in one pass (they share data, weights and mean estimates).
run_smcs2 <- function(X){
  Csum  <- matrix(0, d, d); Csum2 <- matrix(0, d, d)
  logE  <- numeric(d); logE2 <- numeric(d)
  excl  <- rep(FALSE, d); excl2 <- rep(FALSE, d)
  t1a <- Inf; t1b <- Inf
  path <- path2 <- numeric(length(grid)); gi <- 1L
  thr  <- log(1/alpha); Vhat <- NULL; ke <- 0L
  for (t in 1:nTot){
    M <- outer(X[t, ], X[t, ], "-")
    ## epoch refit of the estimated penalty variances (frozen within epochs)
    kk <- findInterval(t, tks)
    if (kk >= 1L && kk != ke){
      ke <- kk
      s2 <- pmax((Csum2 - Csum^2 / (t - 1)) / (t - 2), 1e-10)
      Vhat <- s2 / q1[ke]; diag(Vhat) <- 1
    }
    if (t > t0){
      MUh <- Csum / (t - 1)
      lw  <- eta0 * Csum / sqrt(t - 1); diag(lw) <- -Inf
      lse_w <- rowLSE(lw)
      ## known-variance version
      LAM <- pmax(MUh, 0) / Vd
      A   <- LAM * M - 0.5 * LAM^2 * Vd
      logE <- logE + (rowLSE(lw + A) - lse_w)
      nx <- !excl & (logE >= thr)
      if (nx[1] && !is.finite(t1a)) t1a <- t
      excl <- excl | nx
      ## estimated-variance version (epoch-corrected penalty)
      if (!is.null(Vhat)){
        LAM2 <- pmax(MUh, 0) / Vhat
        A2   <- LAM2 * M - 0.5 * LAM2^2 * Vhat
        logE2 <- logE2 + (rowLSE(lw + A2) - lse_w)
        nx2 <- !excl2 & (logE2 >= thr)
        if (nx2[1] && !is.finite(t1b)) t1b <- t
        excl2 <- excl2 | nx2
      }
    }
    Csum <- Csum + M; Csum2 <- Csum2 + M^2
    if (gi <= length(grid) && t == grid[gi]){
      path[gi] <- sum(!excl); path2[gi] <- sum(!excl2); gi <- gi + 1L }
  }
  list(size = sum(!excl),  rej1 = is.finite(t1a), path = path,
       size2 = sum(!excl2), rej12 = is.finite(t1b), path2 = path2)
}

da_test <- function(X, t, r, adj){
  m  <- floor(t / 2); X1 <- X[1:m, , drop = FALSE]; X2 <- X[(m+1):(2*m), , drop = FALSE]
  m2 <- colMeans(X2)
  if (!adj){ s <- setdiff(order(m2), r)[1] }
  else {
    Dif <- X2 - X2[, r]
    v2  <- pmax(apply(Dif, 2, var), kappa); v2[r] <- Inf
    sc  <- (m2 - m2[r]) / sqrt(v2); sc[r] <- Inf
    s   <- which.min(sc)
  }
  w  <- X1[, r] - X1[, s]
  sqrt(m) * mean(w) > qnorm(1 - alpha) * sd(w)
}

bonf_test <- function(X, t, r){
  Xt <- X[1:t, , drop = FALSE]
  Dif <- Xt[, r] - Xt
  st  <- sqrt(t) * colMeans(Dif) / pmax(apply(Dif, 2, sd), 1e-12); st[r] <- -Inf
  max(st) > qnorm(1 - alpha / (d - 1))
}

set.seed(1000 + round(100 * rho) + 7 * isnull + match(scen, c("a","b","c")))
res <- list(
  smcs_size = numeric(ns), smcs_rej1_any = logical(ns), smcs_path = matrix(0, ns, length(grid)),
  smcs2_size = numeric(ns), smcs2_rej1_any = logical(ns), smcs2_path = matrix(0, ns, length(grid)),
  dap_size = numeric(ns), dap_rej1_fix = logical(ns), dap_rej1_any = logical(ns),
  daa_size = numeric(ns), daa_rej1_fix = logical(ns), daa_rej1_any = logical(ns),
  bon_size = numeric(ns), bon_rej1_fix = logical(ns), bon_rej1_any = logical(ns))

for (s in 1:ns){
  X <- matrix(rnorm(nTot * d), nTot, d) %*% L + matrix(mu, nTot, d, byrow = TRUE)
  sm <- run_smcs2(X)
  res$smcs_size[s]  <- sm$size;  res$smcs_rej1_any[s]  <- sm$rej1;  res$smcs_path[s, ]  <- sm$path
  res$smcs2_size[s] <- sm$size2; res$smcs2_rej1_any[s] <- sm$rej12; res$smcs2_path[s, ] <- sm$path2
  dap <- daa <- bon <- logical(d)
  for (r in 1:d){
    dap[r] <- da_test(X, nTot, r, adj = FALSE)
    daa[r] <- da_test(X, nTot, r, adj = TRUE)
    bon[r] <- bonf_test(X, nTot, r)
  }
  res$dap_size[s] <- sum(!dap); res$daa_size[s] <- sum(!daa); res$bon_size[s] <- sum(!bon)
  res$dap_rej1_fix[s] <- dap[1]; res$daa_rej1_fix[s] <- daa[1]; res$bon_rej1_fix[s] <- bon[1]
  ra <- rp <- rb <- FALSE
  for (t in grid){
    if (!rp) rp <- da_test(X, t, 1, adj = FALSE)
    if (!ra) ra <- da_test(X, t, 1, adj = TRUE)
    if (!rb) rb <- bonf_test(X, t, 1)
    if (rp && ra && rb) break
  }
  res$dap_rej1_any[s] <- rp; res$daa_rej1_any[s] <- ra; res$bon_rej1_any[s] <- rb
}

meta <- list(scen = scen, rho = rho, isnull = isnull, ns = ns, d = d, nTot = nTot,
             alpha = alpha, grid = grid, eta0 = eta0, delta = delta, tks = tks)
saveRDS(list(meta = meta, res = res),
        sprintf("results/bench_%s%d_rho%s.rds", scen, isnull,
                gsub("-", "m", sprintf("%.1f", rho))))
cat(sprintf(
 "%s%s rho=%+.1f ns=%d | setsize SMCS/SMCSest/DAp/DAa/Bon = %.1f/%.1f/%.1f/%.1f/%.1f | rej1 ANY SMCS/SMCSest = %.3f/%.3f | rej1 fix DAp/DAa/Bon = %.3f/%.3f/%.3f | rej1 ANY DAp/DAa/Bon = %.3f/%.3f/%.3f\n",
 scen, ifelse(isnull, "0(null)", "(alt)"), rho, ns,
 mean(res$smcs_size), mean(res$smcs2_size), mean(res$dap_size), mean(res$daa_size), mean(res$bon_size),
 mean(res$smcs_rej1_any), mean(res$smcs2_rej1_any),
 mean(res$dap_rej1_fix), mean(res$daa_rej1_fix), mean(res$bon_rej1_fix),
 mean(res$dap_rej1_any), mean(res$daa_rej1_any), mean(res$bon_rej1_any)))
