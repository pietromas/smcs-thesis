################################################################################
## ch2_bounded_p3.R -- Table 2.1
##
## Conditionally bounded loss differences, m = 4 models (p = m-1 = 3
## comparisons), equicorrelated Gaussian copula.
##
##   d_{i,t} = mu_i + a U_t,  a = 0.8,  U_t uniform(-1,1)^p
##
##   AVG    per-comparison aGRAPA, symmetric average
##   MV     multivariate plug-in bet on Gamma_i
##   AVGor  oracle AVG (true moments)
##   MVor   oracle MV  (true moments)
##
## Usage: Rscript R/ch2_bounded_p3.R [nreps] [n] [alpha] [seed]
##        defaults: 1000 500 0.05 20260726
################################################################################

args  <- commandArgs(trailingOnly = TRUE)
nreps <- if (length(args) >= 1) as.integer(args[1]) else 1000L
n     <- if (length(args) >= 2) as.integer(args[2]) else 500L
alpha <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
seed0 <- if (length(args) >= 4) as.integer(args[4]) else 20260726L

suppressMessages(library(quadprog))

a      <- 0.8
p      <- 3
cvec   <- rep(1, p)
warmup <- 6L
eps    <- 1e-12
thr    <- log(1 / alpha)

cat(sprintf("bounded_p3: m = %d models (p = %d comparisons), n = %d, alpha = %.3f, nreps = %d\n\n",
            p + 1, p, n, alpha, nreps))

## exact correlation of d after the copula transform
rho_d <- function(rz) (6 / pi) * asin(rz / 2)

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

one_rep <- function(mu, Lch, lam_mv_or, lam_av_or){
  Z <- matrix(rnorm(n * p), n, p) %*% Lch
  D <- sweep(a * (2 * pnorm(Z) - 1), 2, mu, "+")

  lg <- numeric(p); lg_or <- numeric(p); lmv <- 0; lmv_or <- 0
  S1 <- numeric(p); S2 <- matrix(0, p, p)
  lam_mv <- numeric(p); tau <- rep(NA_integer_, 4)

  for (t in 1:n){
    d <- D[t, ]
    if (t > warmup){
      n0 <- t - 1
      muh <- S1 / n0
      lam_mv <- solve_bet(muh, S2 / n0)
      lam_av <- pmin(pmax(muh / (diag(S2) / n0 + eps), 0), 1 / cvec)
    } else {
      lam_av <- numeric(p); lam_mv <- numeric(p)
    }
    if (is.na(tau[1])){
      lg <- lg + log(pmax(1 + lam_av * d, eps))
      if (lse(lg) - log(p) >= thr) tau[1] <- t
    }
    if (is.na(tau[2])){
      lmv <- lmv + log(max(1 + sum(lam_mv * d), eps))
      if (lmv >= thr) tau[2] <- t
    }
    if (is.na(tau[3])){
      lg_or <- lg_or + log(pmax(1 + lam_av_or * d, eps))
      if (lse(lg_or) - log(p) >= thr) tau[3] <- t
    }
    if (is.na(tau[4])){
      lmv_or <- lmv_or + log(max(1 + sum(lam_mv_or * d), eps))
      if (lmv_or >= thr) tau[4] <- t
    }
    if (!anyNA(tau)) break
    S1 <- S1 + d; S2 <- S2 + tcrossprod(d)
  }
  tau
}

## scenarios: mean vector and copula equicorrelation
mk <- function(v){ stopifnot(length(v) == p, all(abs(v) <= 1 - a + 1e-12)); v }
scen <- list(
  list(nm = "three equal",  mu = mk(c(.10,.10,.10)), rho = -0.4),
  list(nm = "three equal",  mu = mk(c(.10,.10,.10)), rho =  0.0),
  list(nm = "three equal",  mu = mk(c(.10,.10,.10)), rho =  0.8),
  list(nm = "two equal",    mu = mk(c(.12,.12,0)),   rho =  0.0),
  list(nm = "three graded", mu = mk(c(.15,.10,.05)), rho = -0.4),
  list(nm = "single",       mu = mk(c(.15,0,0)),     rho =  0.0),
  list(nm = "null tied",    mu = mk(c(0,0,0)),       rho =  0.0),
  list(nm = "null tied",    mu = mk(c(0,0,0)),       rho =  0.8),
  list(nm = "null partial", mu = mk(c(0,0,-.15)),    rho =  0.0),
  list(nm = "null partial", mu = mk(c(0,0,-.15)),    rho =  0.8)
)

meth <- c("AVG","MV","AVGor","MVor")
res  <- list()

for (s in seq_along(scen)){
  z <- scen[[s]]; mu <- z$mu; rz <- z$rho
  Rz <- matrix(rz, p, p); diag(Rz) <- 1
  stopifnot(min(eigen(Rz, symmetric = TRUE, only.values = TRUE)$values) > 1e-8)
  Lch   <- chol(Rz)
  Sig_d <- (a^2 / 3) * rho_d(Rz); diag(Sig_d) <- a^2 / 3
  M_true <- Sig_d + tcrossprod(mu)

  lam_mv_or <- solve_bet(mu, M_true)
  lam_av_or <- pmin(pmax(mu / (diag(M_true) + eps), 0), 1 / cvec)

  set.seed(seed0 + 100L * s)
  TAU <- matrix(NA_integer_, nreps, 4, dimnames = list(NULL, meth))
  for (r in 1:nreps) TAU[r, ] <- one_rep(mu, Lch, lam_mv_or, lam_av_or)

  rej <- colMeans(!is.na(TAU))
  ok  <- complete.cases(TAU)
  mt  <- if (any(ok)) colMeans(TAU[ok, , drop = FALSE]) else rep(NA, 4)

  res[[s]] <- list(nm = z$nm, mu = mu, rho = rz, tau = TAU, rej = rej,
                   nact = sum(lam_mv_or > 1e-8))

  if (grepl("^null", z$nm)){
    cat(sprintf("RES %-13s rho=%+.1f NULL | Type-I AVG/MV/AVGor/MVor = %.3f/%.3f/%.3f/%.3f\n",
                z$nm, rz, rej[1], rej[2], rej[3], rej[4]))
  } else {
    cat(sprintf(paste0("RES %-13s rho=%+.1f mu=(%s) act=%d | rej = %.3f/%.3f/%.3f/%.3f",
                       " | meanTau = %.1f/%.1f/%.1f/%.1f\n"),
                z$nm, rz, paste(sprintf("%.2f", mu), collapse = ","),
                sum(lam_mv_or > 1e-8), rej[1], rej[2], rej[3], rej[4],
                mt[1], mt[2], mt[3], mt[4]))
  }
  flush.console()
}

saveRDS(res, "results/bounded_p3_all.rds")
cat("\ndone. wrote results/bounded_p3_all.rds\n")
