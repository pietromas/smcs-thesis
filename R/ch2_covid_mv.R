################################################################################
## ch2_covid_mv.R -- multivariate bet on the Covid-19 Forecast Hub data
##
## Runs the joint bet of Chapter 2 alongside the averaging constructions of
## Chapter 1 on the same stream of Hub forecasts. The loss differences are
## conditionally bounded with the pairwise, time-varying bound
##   c_{ij,t} = max(q, 1-q) |log f_{i,t} - log f_{j,t}|,
## so the polytope Gamma_{i,t} changes every week and the bet is re-solved
## against it at each round.
##
## Note on scope: the data are strongly serially dependent, so the growth-rate
## theory of Chapter 2, which assumes iid loss differences, does not apply here.
## Validity does, since the bounded-case argument uses only the conditional
## mean. We therefore report set sizes and exclusion times, not comparisons
## against theoretical rejection-time predictions.
##
## Usage:  Rscript R/ch2_covid_mv.R [quantile] [alpha] [mmin]
##         defaults: quantile = 0.5, alpha = 0.1, mmin = 15
## Output: results/covid_mv.rds, figures/covid_mv_setsize.pdf
################################################################################
suppressMessages(library(quadprog))

args  <- commandArgs(trailingOnly = TRUE)
qv    <- if (length(args) >= 1) as.numeric(args[1]) else 0.5
alpha <- if (length(args) >= 2) as.numeric(args[2]) else 0.1
mmin  <- if (length(args) >= 3) as.integer(args[3]) else 15
eta0  <- 1
eps   <- 1e-6

load("data/forecasts_death.rda"); load("data/truth.rda")

## ---- data preparation, identical to ch1_covid.R -----------------------------
fd <- as.data.frame(forecasts_death)
fd$quantile <- as.numeric(as.character(fd$quantile))
fd <- fd[!is.na(fd$quantile) & abs(fd$quantile - qv) < 1e-6,
         c("model", "target_end_date", "value")]
fd$target_end_date <- as.Date(fd$target_end_date)
fd$value <- as.numeric(fd$value)
fd <- fd[!is.na(fd$value), ]
fd <- fd[!duplicated(fd[, c("model", "target_end_date")], fromLast = TRUE), ]

tr <- as.data.frame(truth)
tr$target_end_date <- as.Date(tr$target_end_date)
tr$value <- as.numeric(tr$value)
tr <- tr[order(tr$target_end_date), ]
tr <- tr[!duplicated(tr$target_end_date), c("target_end_date", "value")]
fd <- fd[fd$target_end_date %in% tr$target_end_date, ]

dates  <- sort(unique(fd$target_end_date))
models <- sort(unique(fd$model))
Fm <- matrix(NA_real_, length(models), length(dates),
             dimnames = list(models, as.character(dates)))
Fm[cbind(match(fd$model, models), match(fd$target_end_date, dates))] <- fd$value

counts <- colSums(!is.na(Fm))
ok <- counts >= mmin
if (!any(ok)) stop("no week has >= mmin models; lower mmin")
r <- rle(ok); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
runs <- which(r$values); best <- runs[which.max(r$lengths[runs])]
win  <- starts[best]:ends[best]

covw <- rowMeans(!is.na(Fm[, win, drop = FALSE]))
Fm   <- Fm[covw >= 0.9, win, drop = FALSE]
for (i in seq_len(nrow(Fm))){
  v <- Fm[i, ]
  if (anyNA(v)){
    if (is.na(v[1])) v[1] <- v[which(!is.na(v))[1]]
    for (t in 2:length(v)) if (is.na(v[t])) v[t] <- v[t - 1]
    Fm[i, ] <- v
  }
}
dates <- as.Date(colnames(Fm))
y     <- tr$value[match(dates, tr$target_end_date)]
m <- nrow(Fm); Tn <- ncol(Fm); p <- m - 1
cat(sprintf("window %s to %s, T = %d weeks, m = %d models, p = %d comparisons\n",
            min(dates), max(dates), Tn, m, p))

pinball <- function(y, f, q) ifelse(y >= f, q * (y - f), (1 - q) * (f - y))
logF <- log(eps + Fm); logy <- log(eps + y)
Lmat <- matrix(0, m, Tn)
for (t in 1:Tn) Lmat[, t] <- pinball(logy[t], logF[, t], qv)

## ---- multivariate bet on the time-varying polytope --------------------------
warmup <- p + 4     # need more observations than coordinates before betting

solve_bet <- function(muh, Mh, cvec, ridge = 1e-7){
  if (all(muh <= 0)) return(rep(0, length(muh)))
  k <- length(muh)
  out <- try(solve.QP(Dmat = Mh + diag(ridge, k), dvec = muh,
                      Amat = cbind(diag(k), -cvec), bvec = c(rep(0, k), -1),
                      meq = 0), silent = TRUE)
  if (inherits(out, "try-error")) return(rep(0, k))
  lam <- pmax(out$solution, 0)
  s <- sum(cvec * lam); if (s > 1) lam <- lam / s
  lam
}

run_mv <- function(){
  ## one e-process per tested model, each betting jointly on its own vector
  S1 <- vector("list", m); S2 <- vector("list", m)
  for (i in 1:m){ S1[[i]] <- numeric(p); S2[[i]] <- matrix(0, p, p) }
  logE <- numeric(m); excl <- rep(FALSE, m)
  size <- integer(Tn); size[1] <- m
  excl_time <- rep(NA_integer_, m); nact <- matrix(0, m, Tn)

  for (t in 2:Tn){
    for (i in 1:m){
      oth  <- (1:m)[-i]
      dvec <- Lmat[i, t] - Lmat[oth, t]
      cvec <- pmax(max(qv, 1 - qv) * abs(logF[i, t] - logF[oth, t]), eps)
      if (!excl[i]){
        if (t > warmup){
          n0  <- t - 2
          muh <- S1[[i]] / n0
          Mh  <- S2[[i]] / n0
          lam <- solve_bet(muh, Mh, cvec)
        } else lam <- numeric(p)
        nact[i, t] <- sum(lam > 1e-8)
        logE[i] <- logE[i] + log(max(1 + sum(lam * dvec), eps))
        if (logE[i] >= log(1 / alpha)){ excl[i] <- TRUE; excl_time[i] <- t }
      }
      S1[[i]] <- S1[[i]] + dvec
      S2[[i]] <- S2[[i]] + tcrossprod(dvec)
    }
    size[t] <- sum(!excl)
  }
  list(size = size, excl = excl, excl_time = excl_time, logE = logE, nact = nact)
}

## ---- averaging constructions, as in ch1_covid.R ------------------------------
run_avg <- function(weights = c("uniform", "softmax"), fraction = c("fixed", "agrapa")){
  weights <- match.arg(weights); fraction <- match.arg(fraction)
  Cs <- matrix(0, m, m); S1 <- matrix(0, m, m); S2 <- matrix(0, m, m)
  logE_i <- numeric(m); excl <- rep(FALSE, m)
  size <- integer(Tn); size[1] <- m; excl_time <- rep(NA_integer_, m)
  for (t in 2:Tn){
    D  <- outer(Lmat[, t], Lmat[, t], "-")
    Cb <- max(qv, 1 - qv) * abs(outer(logF[, t], logF[, t], "-"))
    if (fraction == "fixed"){
      LAM <- 1 / (2 * Cb + eps)
    } else {
      n0  <- t - 2
      muh <- if (n0 > 0) S1 / n0 else matrix(0, m, m)
      m2h <- if (n0 > 0) S2 / n0 else matrix(1, m, m)
      LAM <- pmin(pmax(muh / (m2h + eps), 0), 1 / (Cb + eps))
    }
    if (weights == "uniform"){
      W <- matrix(1 / (m - 1), m, m)
    } else {
      A <- eta0 * Cs; diag(A) <- -Inf
      A <- A - apply(A, 1, max)
      W <- exp(A); W <- W / rowSums(W)
    }
    diag(W) <- 0; W <- W / rowSums(W)
    fac <- rowSums(W * (1 + LAM * D))
    logE_i <- logE_i + log(pmax(fac, eps))
    newx <- !excl & (logE_i >= log(1 / alpha))
    excl_time[newx] <- t; excl <- excl | newx
    size[t] <- sum(!excl)
    Cs <- Cs + D; S1 <- S1 + D; S2 <- S2 + D^2
  }
  list(size = size, excl = excl, excl_time = excl_time, logE = logE_i)
}

st <- run_avg("uniform", "fixed")
sa <- run_avg("softmax", "agrapa")
mv <- run_mv()
res <- list(standard = st, `softmax+aGRAPA` = sa, multivariate = mv)

cat("\nfinal set sizes:\n")
for (nm in names(res))
  cat(sprintf("  %-16s %2d of %2d | excluded %2d\n",
              nm, res[[nm]]$size[Tn], m, sum(res[[nm]]$excl)))

cat(sprintf("\nmultivariate bet: no bet before week %d; mean number of comparisons\n", warmup))
cat(sprintf("bet on, over weeks %d to %d and surviving models: %.1f of %d\n",
            warmup + 1, Tn, mean(mv$nact[!mv$excl, (warmup+1):Tn]), p))

cat("\nmodels in the final set:\n")
for (nm in names(res))
  cat(sprintf("  %-16s %s\n", nm, paste(rownames(Fm)[!res[[nm]]$excl], collapse = ", ")))

saveRDS(list(res = res, dates = dates, models = rownames(Fm), m = m, Tn = Tn,
             qv = qv, alpha = alpha, warmup = warmup),
        "results/covid_mv.rds")

cols <- c("#7d8a2e", "#c0392b", "#2471a3")
pdf("figures/covid_mv_setsize.pdf", width = 8.5, height = 5)
plot(dates, st$size, type = "s", lwd = 2, col = cols[1], ylim = c(0, m + 1),
     xlab = "target end date", ylab = "size of the confidence set",
     main = sprintf("Covid-19 Hub: averaging against the joint bet (q = %.2f, alpha = %.2f)",
                    qv, alpha))
lines(dates, sa$size, type = "s", lwd = 2, col = cols[2])
lines(dates, mv$size, type = "s", lwd = 2, col = cols[3])
abline(v = dates[warmup], lty = 3)
legend("bottomleft", c("standard", "softmax + aGRAPA", "multivariate"),
       col = cols, lwd = 2, bty = "n")
dev.off()
cat("\nwrote figures/covid_mv_setsize.pdf and results/covid_mv.rds\n")
