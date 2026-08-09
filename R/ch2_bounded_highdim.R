################################################################################
## ch2_bounded_highdim.R -- Tables 2.2 and 2.3
##
## Conditionally bounded loss differences, m = 10 models (p = 9 comparisons),
## mixed dependence from a two-factor loading model.
##
##   d_{i,t} = mu_i + a U_t,  a = 0.8,  U_t uniform(-1,1)^p, Gaussian copula
##   comparisons 1-3 block P : positive loading on factor 1
##   comparisons 4-6 block N : negative loading on factor 1  (P-N negative)
##   comparisons 7-8 block S : load on factor 2 only
##   comparison  9   ISOLATED: |rho| <= 0.03 against everything
##
## Usage: Rscript R/ch2_bounded_highdim.R [nreps] [n] [alpha] [refit] [seed]
##        defaults: 3000 500 0.05 1 20260723
################################################################################

args   <- commandArgs(trailingOnly = TRUE)
nreps  <- if (length(args) >= 1) as.integer(args[1]) else 3000L
n      <- if (length(args) >= 2) as.integer(args[2]) else 500L
alpha  <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
refit  <- if (length(args) >= 4) as.integer(args[4]) else 1L
seed0  <- if (length(args) >= 5) as.integer(args[5]) else 20260723L

suppressMessages(library(quadprog))

a      <- 0.8
p      <- 9
cvec   <- rep(1, p)
warmup <- 12L
eps    <- 1e-12
thr    <- log(1 / alpha)

cat(sprintf("bounded_highdim: m = %d models (p = %d comparisons), n = %d, alpha = %.3f\n",
            p + 1, p, n, alpha))
cat(sprintf("                 nreps = %d, refit every %d steps, warmup = %d, seed = %d\n\n",
            nreps, refit, warmup, seed0))

## ---------------------------------------------------------------------------
## 1. Mixed correlation structure (two factors + one isolated comparison)
## ---------------------------------------------------------------------------
B <- rbind(
  c( 0.85,  0.15), c( 0.75,  0.25), c( 0.65,  0.10),
  c(-0.85,  0.20), c(-0.75, -0.15), c(-0.65,  0.25),
  c( 0.05,  0.60), c( 0.00, -0.55), c( 0.02,  0.04)
)
stopifnot(nrow(B) == p, all(rowSums(B^2) < 1))
Rz <- B %*% t(B) + diag(1 - rowSums(B^2))
ev <- eigen(Rz, symmetric = TRUE, only.values = TRUE)$values
stopifnot(min(ev) > 1e-8)

Ru      <- (6 / pi) * asin(Rz / 2); diag(Ru) <- 1
Sigma_d <- (a^2 / 3) * Ru
Lchol   <- chol(Rz)

offd <- Ru[upper.tri(Ru)]
cat("correlation structure of the loss differences:\n")
cat(sprintf("  min = %+.3f | max = %+.3f | mean = %+.3f | %d negative, %d |r|<0.1\n",
            min(offd), max(offd), mean(offd), sum(offd < 0), sum(abs(offd) < 0.1)))
cat(sprintf("  min eigenvalue of Rz = %.4f\n", min(ev)))
cat(sprintf("  block means: PP = %+.2f | NN = %+.2f | PN = %+.2f | S(7,8) = %+.2f\n",
            mean(Ru[1:3, 1:3][upper.tri(diag(3))]),
            mean(Ru[4:6, 4:6][upper.tri(diag(3))]),
            mean(Ru[1:3, 4:6]), Ru[7, 8]))
cat(sprintf("  CONTROLLED PAIR: r(1,4) = %+.3f | r(1,2) = %+.3f | ISOLATED comp 9 max|r| = %.3f\n",
            Ru[1, 4], Ru[1, 2], max(abs(Ru[9, -9]))))

set.seed(seed0)
Zc <- matrix(rnorm(2e5 * p), 2e5, p) %*% Lchol
cat(sprintf("  MC check: max |Sigma_hat - Sigma_exact| = %.5f\n\n",
            max(abs(cov(a * (2 * pnorm(Zc) - 1)) - Sigma_d))))
rm(Zc)

## ---------------------------------------------------------------------------
## 2. Growth-optimal bet on Gamma_i = {lambda >= 0, c'lambda <= 1}
## ---------------------------------------------------------------------------
solve_bet <- function(muh, Mh, ridge = 1e-7){
  if (all(muh <= 0)) return(rep(0, p))
  out <- try(solve.QP(Dmat = Mh + diag(ridge + 1e-10, p), dvec = muh,
                      Amat = cbind(diag(p), -cvec), bvec = c(rep(0, p), -1),
                      meq = 0), silent = TRUE)
  if (inherits(out, "try-error")) return(rep(0, p))
  lam <- pmax(out$solution, 0)
  s <- sum(cvec * lam); if (s > 1) lam <- lam / s
  lam
}

lse <- function(v){ mx <- max(v); mx + log(sum(exp(v - mx))) }

## ---------------------------------------------------------------------------
## 3. One repetition
## ---------------------------------------------------------------------------
one_rep <- function(mu, lam_mv_or, lam_av_or){
  Z <- matrix(rnorm(n * p), n, p) %*% Lchol
  D <- sweep(a * (2 * pnorm(Z) - 1), 2, mu, "+")

  logE_j <- numeric(p); logE_j_or <- numeric(p)
  logE_mv <- 0;         logE_mv_or <- 0
  S1 <- numeric(p);     S2 <- matrix(0, p, p)
  lam_mv <- numeric(p)
  tau <- rep(NA_integer_, 4)

  for (t in 1:n){
    d <- D[t, ]

    if (t > warmup){
      n0  <- t - 1
      muh <- S1 / n0
      if (t == warmup + 1L || ((t - warmup - 1L) %% refit == 0L)){
        lam_mv <- solve_bet(muh, S2 / n0)
      }
      lam_av <- pmin(pmax(muh / (diag(S2) / n0 + eps), 0), 1 / cvec)
    } else {
      lam_av <- numeric(p); lam_mv <- numeric(p)
    }

    if (is.na(tau[1])){
      logE_j <- logE_j + log(pmax(1 + lam_av * d, eps))
      if (lse(logE_j) - log(p) >= thr) tau[1] <- t
    }
    if (is.na(tau[2])){
      logE_mv <- logE_mv + log(max(1 + sum(lam_mv * d), eps))
      if (logE_mv >= thr) tau[2] <- t
    }
    if (is.na(tau[3])){
      logE_j_or <- logE_j_or + log(pmax(1 + lam_av_or * d, eps))
      if (lse(logE_j_or) - log(p) >= thr) tau[3] <- t
    }
    if (is.na(tau[4])){
      logE_mv_or <- logE_mv_or + log(max(1 + sum(lam_mv_or * d), eps))
      if (logE_mv_or >= thr) tau[4] <- t
    }

    if (!anyNA(tau)) break
    S1 <- S1 + d; S2 <- S2 + tcrossprod(d)
  }
  tau
}

## ---------------------------------------------------------------------------
## 4. Scenarios
## ---------------------------------------------------------------------------
mk <- function(v){ stopifnot(length(v) == p, all(abs(v) <= 1 - a + 1e-12)); v }

scen <- list(
  pair_compl   = mk(c(.15,0,0, .15,0,0, 0,0,0)),
  pair_redund  = mk(c(.15,.15,0, 0,0,0, 0,0,0)),
  six_compl    = mk(c(.12,.12,.12, .12,.12,.12, 0,0,0)),
  three_redund = mk(c(.12,.12,.12, 0,0,0, 0,0,0)),
  nine_graded  = mk(c(.14,.11,.08, .14,.11,.08, .14,.11,.08)),
  four_mixed   = mk(c(.13,0,0, .13,0,0, .13,.13,0)),
  single_corr  = mk(c(.15,0,0, 0,0,0, 0,0,0)),
  single_iso   = mk(c(0,0,0, 0,0,0, 0,0,.15)),
  null_tied    = mk(numeric(p)),
  null_partial = mk(c(0,0,0, -.10,-.10,-.10, 0,0,0))
)
is_null <- grepl("^null", names(scen))

## ---------------------------------------------------------------------------
## 5. Run
## ---------------------------------------------------------------------------
meth <- c("AVG", "MV", "AVGor", "MVor")
res  <- list()

for (s in seq_along(scen)){
  nm <- names(scen)[s]; mu <- scen[[s]]
  M_true     <- Sigma_d + tcrossprod(mu)
  lam_mv_or  <- solve_bet(mu, M_true)
  lam_av_or  <- pmin(pmax(mu / (diag(M_true) + eps), 0), 1 / cvec)
  nact       <- sum(lam_mv_or > 1e-8)

  set.seed(seed0 + 1000L * s)
  TAU <- matrix(NA_integer_, nreps, 4, dimnames = list(NULL, meth))
  t0  <- Sys.time()
  for (r in 1:nreps) TAU[r, ] <- one_rep(mu, lam_mv_or, lam_av_or)
  el  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  rej  <- colMeans(!is.na(TAU))
  all4 <- complete.cases(TAU)
  mt   <- if (any(all4)) colMeans(TAU[all4, , drop = FALSE]) else rep(NA, 4)
  md   <- if (any(all4)) apply(TAU[all4, , drop = FALSE], 2, median) else rep(NA, 4)

  res[[nm]] <- list(mu = mu, tau = TAU, rej = rej, n = n, alpha = alpha,
                    lam_mv_or = lam_mv_or, lam_av_or = lam_av_or,
                    nact = nact, Sigma_d = Sigma_d, nreps = nreps)
  saveRDS(res[[nm]], sprintf("results/mvhd_%s.rds", nm))

  if (is_null[s]){
    cat(sprintf("RES scen=%-13s NULL | anytime Type-I AVG/MV/AVGor/MVor = %.3f/%.3f/%.3f/%.3f | %.1f min\n",
                nm, rej[1], rej[2], rej[3], rej[4], el))
  } else {
    cat(sprintf(paste0("RES scen=%-13s nviol=%d act=%d | rej AVG/MV/AVGor/MVor = %.3f/%.3f/%.3f/%.3f",
                       " | meanTau = %.1f/%.1f/%.1f/%.1f | medTau = %.0f/%.0f/%.0f/%.0f",
                       " | ratio AVG/MV = %.2f | %.1f min\n"),
                nm, sum(mu > 0), nact, rej[1], rej[2], rej[3], rej[4],
                mt[1], mt[2], mt[3], mt[4], md[1], md[2], md[3], md[4],
                mt[1] / mt[2], el))
  }
  flush.console()
}

saveRDS(list(res = res, Rz = Rz, Ru = Ru, Sigma_d = Sigma_d, B = B,
             n = n, alpha = alpha, nreps = nreps, a = a, p = p),
        "results/mvhd_all.rds")
cat("\ndone. wrote results/mvhd_<scenario>.rds and results/mvhd_all.rds\n")
