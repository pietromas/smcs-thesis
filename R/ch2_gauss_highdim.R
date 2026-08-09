################################################################################
## ch2_gauss_highdim.R -- Tables 2.6, 2.7, 2.8
##
## Gaussian loss differences, m = 10 models (p = 9 comparisons), mixed
## dependence. Compares the per-comparison averaging construction with the
## multivariate bet, each with its oracle and its uncorrected plug-in, all
## under the epoch correction of Section 2.5.
##
##   d_{i,t} ~ N(mu_i, Sigma),  Sigma = R  (unit variances)
##   e-variable: exp(lambda' d - 0.5 lambda' SigmaHat lambda),  lambda >= 0
##
##   AVG    per-comparison, epoch-corrected variances, symmetric average
##   MV     multivariate, epoch-corrected Wishart-inflated covariance
##   AVGor  oracle AVG (true sigma^2, no inflation)
##   MVor   oracle MV  (true Sigma,   no inflation)
##   AVGnv  naive AVG  (plain sample variances)   [validity check]
##   MVnv   naive MV   (plain sample covariance)  [validity check]
##
## Usage: Rscript R/ch2_gauss_highdim.R [nreps] [n] [alpha] [refit] [delta] [seed]
##        defaults: 2000 1000 0.05 1 0.05 20260724
################################################################################

args   <- commandArgs(trailingOnly = TRUE)
nreps  <- if (length(args) >= 1) as.integer(args[1]) else 2000L
n      <- if (length(args) >= 2) as.integer(args[2]) else 1000L
alpha  <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
refit  <- if (length(args) >= 4) as.integer(args[4]) else 1L
delta  <- if (length(args) >= 5) as.numeric(args[5]) else 0.05
seed0  <- if (length(args) >= 6) as.integer(args[6]) else 20260724L

suppressMessages(library(quadprog))

p        <- 9
eps      <- 1e-12
thr      <- log(1 / alpha)
delta_mv <- delta / 2
delta_un <- delta / 2

cat(sprintf("gauss_highdim: m = %d models (p = %d comparisons), n = %d, alpha = %.3f\n",
            p + 1, p, n, alpha))
cat(sprintf("               nreps = %d, bet refit every %d steps, delta = %.3f",
            nreps, refit, delta))
cat(sprintf(" (MV %.4f / UNI %.4f), seed = %d\n\n", delta_mv, delta_un, seed0))

## ---------------------------------------------------------------------------
## 1. Correlation structure: identical to the bounded study
## ---------------------------------------------------------------------------
B <- rbind(
  c( 0.85,  0.15), c( 0.75,  0.25), c( 0.65,  0.10),
  c(-0.85,  0.20), c(-0.75, -0.15), c(-0.65,  0.25),
  c( 0.05,  0.60), c( 0.00, -0.55), c( 0.02,  0.04)
)
stopifnot(nrow(B) == p, all(rowSums(B^2) < 1))
Rz <- B %*% t(B) + diag(1 - rowSums(B^2))
R  <- (6 / pi) * asin(Rz / 2); diag(R) <- 1
Sigma <- R
Lch   <- chol(Sigma)
stopifnot(min(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values) > 1e-8)

offd <- R[upper.tri(R)]
cat("correlation structure of the loss differences:\n")
cat(sprintf("  min = %+.3f | max = %+.3f | mean = %+.3f | %d negative, %d |r|<0.1\n",
            min(offd), max(offd), mean(offd), sum(offd < 0), sum(abs(offd) < 0.1)))
cat(sprintf("  block means: PP = %+.2f | NN = %+.2f | PN = %+.2f | S(7,8) = %+.2f\n",
            mean(R[1:3,1:3][upper.tri(diag(3))]), mean(R[4:6,4:6][upper.tri(diag(3))]),
            mean(R[1:3,4:6]), R[7,8]))
cat(sprintf("  CONTROLLED PAIR: rho(1,4) = %+.3f | rho(1,2) = %+.3f | ISOLATED comp 9 max|rho| = %.3f\n\n",
            R[1,4], R[1,2], max(abs(R[9,-9]))))

## ---------------------------------------------------------------------------
## 2. Epoch schedule and inflation quantiles
## ---------------------------------------------------------------------------
tk_all <- 2^((1:20) + 2); tk_all <- tk_all[tk_all <= n]
tk     <- tk_all[tk_all >= p + 2]     # need a non-singular Wishart
K      <- length(tk)
stopifnot(K >= 2)
nu_k   <- tk - 2
eps_mv <- delta_mv / K
eps_un <- delta_un / (K * p)

wishart_qmin <- function(nu, p, e, nsim = 2e5, chunk = 2e4, sd = 7L){
  set.seed(sd + nu)
  out <- numeric(nsim); done <- 0L
  while (done < nsim){
    b <- min(chunk, nsim - done)
    W <- stats::rWishart(b, nu, diag(p))
    out[(done + 1L):(done + b)] <-
      apply(W, 3, function(M) min(eigen(M / nu, symmetric = TRUE,
                                        only.values = TRUE)$values))
    done <- done + b
  }
  unname(quantile(out, e, type = 1))
}

cat(sprintf("epoch schedule (K = %d), eps_MV = %.5f, eps_UNI = %.6f\n", K, eps_mv, eps_un))
q_mv <- numeric(K); q_un <- numeric(K)
for (k in 1:K){
  q_mv[k] <- wishart_qmin(nu_k[k], p, eps_mv)
  q_un[k] <- qchisq(eps_un, nu_k[k]) / nu_k[k]
  cat(sprintf("  t_%d = %4d (nu = %4d) | q_MV = %.4f -> inflation %6.2f | q_UNI = %.4f -> inflation %5.2f\n",
              k, tk[k], nu_k[k], q_mv[k], 1/q_mv[k], q_un[k], 1/q_un[k]))
}
cat("\n")

epoch_of <- integer(n)
for (t in 1:n) epoch_of[t] <- sum(tk <= t)

## ---------------------------------------------------------------------------
## 3. Growth-optimal bet on the nonnegative orthant
## ---------------------------------------------------------------------------
solve_bet <- function(muh, Sh, ridge = 1e-8){
  if (all(muh <= 0)) return(rep(0, p))
  out <- try(solve.QP(Dmat = Sh + diag(ridge, p), dvec = muh,
                      Amat = diag(p), bvec = rep(0, p), meq = 0), silent = TRUE)
  if (inherits(out, "try-error")) return(rep(0, p))
  pmax(out$solution, 0)
}

lse <- function(v){ mx <- max(v); mx + log(sum(exp(v - mx))) }

## ---------------------------------------------------------------------------
## 4. One repetition: six strategies on the same stream
## ---------------------------------------------------------------------------
meth <- c("AVG", "MV", "AVGor", "MVor", "AVGnv", "MVnv")

one_rep <- function(mu, lam_mv_or, lam_av_or){
  D <- matrix(rnorm(n * p), n, p) %*% Lch
  D <- sweep(D, 2, mu, "+")

  lg <- numeric(p); lg_or <- numeric(p); lg_nv <- numeric(p)
  lmv <- 0; lmv_or <- 0; lmv_nv <- 0
  tau <- rep(NA_integer_, 6)

  Sh <- NULL; sh <- NULL; Sh_nv <- NULL; sh_nv <- NULL
  lam_mv <- numeric(p); lam_av <- numeric(p)
  lam_mv_nv <- numeric(p); lam_av_nv <- numeric(p)
  cur_ep <- 0L

  for (t in 1:n){
    d <- D[t, ]

    ## epoch refit of the frozen penalty, using data before t
    ke <- epoch_of[t]
    if (ke > cur_ep){
      cur_ep <- ke
      Sfull  <- cov(D[1:(tk[ke] - 1), , drop = FALSE])
      Sh_nv  <- Sfull;                sh_nv <- diag(Sfull)
      Sh     <- Sfull / q_mv[ke];     sh    <- diag(Sfull) / q_un[ke]
    }

    ## bet directions
    if (ke > 0 && ((t - tk[1]) %% refit == 0L || t == tk[ke])){
      muh       <- colMeans(D[1:(t - 1), , drop = FALSE])
      lam_mv    <- solve_bet(muh, Sh)
      lam_av    <- pmax(muh / sh, 0)
      lam_mv_nv <- solve_bet(muh, Sh_nv)
      lam_av_nv <- pmax(muh / sh_nv, 0)
    }

    if (ke > 0){
      if (is.na(tau[1])){
        lg <- lg + (lam_av * d - 0.5 * lam_av^2 * sh)
        if (lse(lg) - log(p) >= thr) tau[1] <- t
      }
      if (is.na(tau[2])){
        lmv <- lmv + sum(lam_mv * d) - 0.5 * drop(lam_mv %*% Sh %*% lam_mv)
        if (lmv >= thr) tau[2] <- t
      }
      if (is.na(tau[5])){
        lg_nv <- lg_nv + (lam_av_nv * d - 0.5 * lam_av_nv^2 * sh_nv)
        if (lse(lg_nv) - log(p) >= thr) tau[5] <- t
      }
      if (is.na(tau[6])){
        lmv_nv <- lmv_nv + sum(lam_mv_nv * d) -
                  0.5 * drop(lam_mv_nv %*% Sh_nv %*% lam_mv_nv)
        if (lmv_nv >= thr) tau[6] <- t
      }
    }
    ## oracles bet from t = 1, true moments, no inflation
    if (is.na(tau[3])){
      lg_or <- lg_or + (lam_av_or * d - 0.5 * lam_av_or^2 * diag(Sigma))
      if (lse(lg_or) - log(p) >= thr) tau[3] <- t
    }
    if (is.na(tau[4])){
      lmv_or <- lmv_or + sum(lam_mv_or * d) -
                0.5 * drop(lam_mv_or %*% Sigma %*% lam_mv_or)
      if (lmv_or >= thr) tau[4] <- t
    }

    if (!anyNA(tau)) break
  }
  tau
}

## ---------------------------------------------------------------------------
## 5. Scenarios (unit noise scale)
## ---------------------------------------------------------------------------
mk <- function(v){ stopifnot(length(v) == p); v }

scen <- list(
  pair_compl   = mk(c(.22,0,0, .22,0,0, 0,0,0)),
  pair_redund  = mk(c(.22,.22,0, 0,0,0, 0,0,0)),
  six_compl    = mk(c(.18,.18,.18, .18,.18,.18, 0,0,0)),
  three_redund = mk(c(.18,.18,.18, 0,0,0, 0,0,0)),
  nine_graded  = mk(c(.22,.17,.12, .22,.17,.12, .22,.17,.12)),
  four_mixed   = mk(c(.20,0,0, .20,0,0, .20,.20,0)),
  single_corr  = mk(c(.24,0,0, 0,0,0, 0,0,0)),
  single_iso   = mk(c(0,0,0, 0,0,0, 0,0,.24)),
  null_tied    = mk(numeric(p)),
  null_partial = mk(c(0,0,0, -.15,-.15,-.15, 0,0,0))
)
is_null <- grepl("^null", names(scen))

## ---------------------------------------------------------------------------
## 6. Run
## ---------------------------------------------------------------------------
res <- list()
for (s in seq_along(scen)){
  nm <- names(scen)[s]; mu <- scen[[s]]
  lam_mv_or <- solve_bet(mu, Sigma)
  lam_av_or <- pmax(mu / diag(Sigma), 0)
  nact      <- sum(lam_mv_or > 1e-8)

  set.seed(seed0 + 1000L * s)
  TAU <- matrix(NA_integer_, nreps, 6, dimnames = list(NULL, meth))
  t0  <- Sys.time()
  for (r in 1:nreps) TAU[r, ] <- one_rep(mu, lam_mv_or, lam_av_or)
  el  <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  rej  <- colMeans(!is.na(TAU))
  ok4  <- complete.cases(TAU[, 1:4, drop = FALSE])
  mt   <- if (any(ok4)) colMeans(TAU[ok4, 1:4, drop = FALSE]) else rep(NA, 4)
  md   <- if (any(ok4)) apply(TAU[ok4, 1:4, drop = FALSE], 2, median) else rep(NA, 4)

  res[[nm]] <- list(mu = mu, tau = TAU, rej = rej, nact = nact, n = n,
                    alpha = alpha, delta = delta, lam_mv_or = lam_mv_or,
                    Sigma = Sigma, tk = tk, q_mv = q_mv, q_un = q_un,
                    nreps = nreps)
  saveRDS(res[[nm]], sprintf("results/gshd_%s.rds", nm))

  if (is_null[s]){
    cat(sprintf(paste0("RES scen=%-13s NULL | anytime Type-I AVG/MV = %.3f/%.3f",
                       " | oracle AVGor/MVor = %.3f/%.3f | NAIVE AVGnv/MVnv = %.3f/%.3f | %.1f min\n"),
                nm, rej[1], rej[2], rej[3], rej[4], rej[5], rej[6], el))
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

saveRDS(list(res = res, R = R, Sigma = Sigma, B = B, tk = tk, q_mv = q_mv,
             q_un = q_un, n = n, alpha = alpha, delta = delta, nreps = nreps,
             p = p, refit = refit),
        "results/gshd_all.rds")
cat("\ndone. wrote results/gshd_<scenario>.rds and results/gshd_all.rds\n")
