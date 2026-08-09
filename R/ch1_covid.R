################################################################################
## covid_case_study_v2.R -- robust version. Same four strategies, but with a
## proper evaluation-window selection for the messy Hub data, plus diagnostics
## printed at every stage so problems are visible.
##
## Window selection: find the longest CONTIGUOUS run of weekly dates on which
## at least `mmin` models all have forecasts, then keep exactly the models that
## are complete on that window.
##
## Usage:  Rscript covid_case_study_v2.R [quantile] [alpha] [mmin]
##         defaults: quantile = 0.5, alpha = 0.1, mmin = 15
################################################################################

args  <- commandArgs(trailingOnly = TRUE)
qv    <- if (length(args) >= 1) as.numeric(args[1]) else 0.5
alpha <- if (length(args) >= 2) as.numeric(args[2]) else 0.1
mmin  <- if (length(args) >= 3) as.integer(args[3]) else 15
eta0  <- 1
eps   <- 1e-6

load("forecasts_death.rda"); load("truth.rda")

fd <- as.data.frame(forecasts_death)
cat(sprintf("raw forecasts: %d rows, %d models, quantile column class = %s\n",
            nrow(fd), length(unique(fd$model)), class(fd$quantile)[1]))
fd$quantile <- as.numeric(as.character(fd$quantile))
fd <- fd[!is.na(fd$quantile) & abs(fd$quantile - qv) < 1e-6,
         c("model", "target_end_date", "value")]
fd$target_end_date <- as.Date(fd$target_end_date)
fd$value <- as.numeric(fd$value)
fd <- fd[!is.na(fd$value), ]
fd <- fd[!duplicated(fd[, c("model", "target_end_date")], fromLast = TRUE), ]
cat(sprintf("after quantile filter: %d rows, %d models, dates %s to %s\n",
            nrow(fd), length(unique(fd$model)),
            min(fd$target_end_date), max(fd$target_end_date)))

tr <- as.data.frame(truth)
tr$target_end_date <- as.Date(tr$target_end_date)
tr$value <- as.numeric(tr$value)
tr <- tr[order(tr$target_end_date), ]
tr <- tr[!duplicated(tr$target_end_date), c("target_end_date", "value")]
fd <- fd[fd$target_end_date %in% tr$target_end_date, ]

## models x dates matrix
dates  <- sort(unique(fd$target_end_date))
models <- sort(unique(fd$model))
Fm <- matrix(NA_real_, length(models), length(dates),
             dimnames = list(models, as.character(dates)))
Fm[cbind(match(fd$model, models), match(fd$target_end_date, dates))] <- fd$value

## ---- window selection: longest contiguous run with >= mmin complete models --
counts <- colSums(!is.na(Fm))
cat("models available per week (5-number summary): ",
    paste(round(fivenum(counts)), collapse = " "), "\n")
ok <- counts >= mmin
if (!any(ok)) stop("no week has >= mmin models; lower mmin (3rd argument)")
r <- rle(ok); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
runs <- which(r$values)
best <- runs[which.max(r$lengths[runs])]
win  <- starts[best]:ends[best]
## models with >= 90% coverage on the window; small gaps filled by carrying the
## last submitted forecast forward (predictable, so validity is unaffected)
covw <- rowMeans(!is.na(Fm[, win, drop = FALSE]))
keep <- covw >= 0.9
Fm   <- Fm[keep, win, drop = FALSE]
nimp <- sum(is.na(Fm))
for (i in seq_len(nrow(Fm))){
  v <- Fm[i, ]
  if (anyNA(v)){
    if (is.na(v[1])) v[1] <- v[which(!is.na(v))[1]]
    for (t in 2:length(v)) if (is.na(v[t])) v[t] <- v[t - 1]
    Fm[i, ] <- v
  }
}
cat(sprintf("imputed %d missing forecasts (%.1f%%) by carrying the last forecast forward\n",
            nimp, 100 * nimp / length(Fm)))
dates <- as.Date(colnames(Fm))
y     <- tr$value[match(dates, tr$target_end_date)]
m <- nrow(Fm); Tn <- ncol(Fm)
cat(sprintf("evaluation window: %s to %s, T = %d weeks, m = %d models kept\n",
            min(dates), max(dates), Tn, m))
if (m < 3 || Tn < 10) stop("window too small; adjust mmin")

## ---------------- losses on the log scale ------------------------------------
pinball <- function(y, f, q) ifelse(y >= f, q * (y - f), (1 - q) * (f - y))
logF <- log(eps + Fm); logy <- log(eps + y)
Lmat <- matrix(0, m, Tn)
for (t in 1:Tn) Lmat[, t] <- pinball(logy[t], logF[, t], qv)

## ---------------- four strategies, marginal coverage --------------------------
run_strategy <- function(weights = c("uniform", "softmax", "poly"),
                         fraction = c("fixed", "agrapa")){
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
    } else if (weights == "softmax"){
      A <- eta0 * Cs; diag(A) <- -Inf
      A <- A - apply(A, 1, max)
      W <- exp(A); W <- W / rowSums(W)
    } else {
      P <- pmax(Cs, 0)^2; diag(P) <- 0
      rs <- rowSums(P)
      W <- P / pmax(rs, eps)
      W[rs == 0, ] <- 1 / (m - 1)
    }
    diag(W) <- 0; W <- W / rowSums(W)
    fac <- rowSums(W * (1 + LAM * D))
    logE_i <- logE_i + log(pmax(fac, eps))
    newx <- !excl & (logE_i >= log(1 / alpha))
    excl_time[newx] <- t
    excl <- excl | newx
    size[t] <- sum(!excl)
    Cs <- Cs + D; S1 <- S1 + D; S2 <- S2 + D^2
  }
  list(size = size, excl = excl, excl_time = excl_time, logE = logE_i)
}

st <- run_strategy("uniform", "fixed")
sf <- run_strategy("softmax", "fixed")
sa <- run_strategy("softmax", "agrapa")
pa <- run_strategy("poly",    "agrapa")
res <- list(standard = st, `softmax+fixed` = sf,
            `softmax+aGRAPA` = sa, `poly+aGRAPA` = pa)

cat("\nfinal set sizes and e-process growth diagnostics:\n")
for (nm in names(res)){
  o <- res[[nm]]
  cat(sprintf("  %-15s final size %2d of %2d | excluded %2d | max final logE = %6.2f (threshold %.2f)\n",
              nm, o$size[Tn], m, sum(o$excl), max(o$logE), log(1/alpha)))
}
cat("\nmodels in the final set (softmax + aGRAPA):\n  ",
    paste(rownames(Fm)[!sa$excl], collapse = ", "), "\n")

cols <- c("#7d8a2e", "#2471a3", "#c0392b", "#e67e22")
pdf("covid_setsize.pdf", width = 8.5, height = 5)
plot(dates, st$size, type = "s", lwd = 2, col = cols[1], ylim = c(0, m + 1),
     xlab = "target end date", ylab = "size of the confidence set",
     main = sprintf("Covid-19 Forecast Hub, 1 wk ahead inc death (q = %.2f, alpha = %.2f)",
                    qv, alpha))
lines(dates, sf$size, type = "s", lwd = 2, col = cols[2])
lines(dates, sa$size, type = "s", lwd = 2, col = cols[3])
lines(dates, pa$size, type = "s", lwd = 2, col = cols[4])
legend("bottomleft", names(res), col = cols, lwd = 2, bty = "n")
dev.off()
cat("\nwrote covid_setsize.pdf\n")
