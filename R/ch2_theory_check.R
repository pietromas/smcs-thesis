################################################################################
## ch2_theory_check.R -- Table 2.4
## Theoretical rejection-time bounds against observed oracle times.
## Run after ch2_bounded_highdim.R.
################################################################################
suppressMessages(library(quadprog))

A  <- readRDS("results/mvhd_all.rds")
R  <- A$Ru; Rz <- A$Rz; a <- A$a; p <- A$p
s2 <- a^2 / 3; alpha <- 0.05; L <- log(1 / alpha); Sig <- s2 * R
cvec <- rep(1, p)

set.seed(1); NMC <- 1e6
U <- 2 * pnorm(matrix(rnorm(NMC * p), NMC, p) %*% chol(Rz)) - 1

rate_exact <- function(lam, mu) mean(log1p(drop(U %*% (a * lam)) + sum(lam * mu)))

rate_single <- function(mu){
  lam <- min(max(mu / (s2 + mu^2), 0), 1)
  integrate(function(u) log1p(lam * (mu + a * u)) / 2, -1, 1)$value
}

bet_mv <- function(mu){
  M <- Sig + outer(mu, mu)
  o <- solve.QP(M + diag(1e-9, p), mu, cbind(diag(p), -cvec), c(rep(0, p), -1))
  lam <- pmax(o$solution, 0); s <- sum(lam); if (s > 1) lam <- lam / s
  lam
}

cat(sprintf("%-15s %2s %8s %8s %8s %8s\n",
            "scenario", "k", "tAVGth", "AVGor", "tMVth", "MVor"))
for (nm in names(A$res)){
  z <- A$res[[nm]]; mu <- z$mu
  if (all(mu <= 0)) next
  k    <- sum(abs(mu - max(mu)) < 1e-12)
  lmax <- max(sapply(mu[mu > 0], rate_single))
  b    <- log(p / alpha - p + 1)
  tA   <- if (k == 1) b / lmax else NA
  tM   <- L / rate_exact(bet_mv(mu), mu)
  o    <- colMeans(z$tau[complete.cases(z$tau), , drop = FALSE])
  cat(sprintf("%-15s %2d %8s %8.1f %8.1f %8.1f\n", nm, k,
              ifelse(is.na(tA), "---", sprintf("%.1f", tA)),
              o["AVGor"], tM, o["MVor"]))
}
