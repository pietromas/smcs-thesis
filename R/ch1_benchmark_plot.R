################################################################################
## ch1_benchmark_plot.R -- Figures 1.2 and 1.3
##
## Figures including the fifth method SMCS-w-est (estimated, epoch-corrected
## variances). Reads the results/bench_*.rds files written by ch1_benchmark.R,
## which must be run for every combination of scenario, rho and null first.
##
## Usage:  Rscript R/ch1_benchmark_plot.R
## Output: figures/benchmark_setsize.pdf, figures/benchmark_type1.pdf
################################################################################
rhos  <- c(0, 0.4, 0.8); scens <- c("a", "b", "c")
lab   <- c(a = "(a) dense near-ties", b = "(b) graded competitors",
           c = "(c) small gap, far mass")
meths <- c("SMCS-w (known var)", "SMCS-w-est", "DA-plug", "DA-adj", "Bonferroni")
cols  <- c("#c0392b", "#e67e22", "#2471a3", "#7fb3d5", "#7d8a2e")

fname <- function(scen, nul, rho)
  sprintf("bench_%s%d_rho%s.rds", scen, nul, gsub("-", "m", sprintf("%.1f", rho)))

grab <- function(scen, nul, rho){
  f <- file.path("results", fname(scen, nul, rho))
  if (!file.exists(f)){ warning("missing ", f); return(NULL) }
  readRDS(f)
}

pdf("figures/benchmark_setsize.pdf", width = 10.5, height = 3.8)
par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.6, 0.8))
for (sc in scens){
  M <- matrix(NA, 5, length(rhos))
  for (j in seq_along(rhos)){
    o <- grab(sc, 0, rhos[j]); if (is.null(o)) next
    M[, j] <- c(mean(o$res$smcs_size), mean(o$res$smcs2_size),
                mean(o$res$dap_size), mean(o$res$daa_size), mean(o$res$bon_size))
  }
  yl <- c(0, max(100 * (sc == "a"), max(M, na.rm = TRUE) * 1.15))
  plot(NULL, xlim = range(rhos) + c(-.05, .05), ylim = yl,
       xlab = expression(rho[loss]), ylab = "average set size at t = 2n",
       main = lab[sc], xaxt = "n")
  axis(1, at = rhos)
  for (k in 1:5){ lines(rhos, M[k, ], col = cols[k], lwd = 2)
                  points(rhos, M[k, ], col = cols[k], pch = 16, cex = 1.2) }
  if (sc == "a") legend("bottomleft", meths, col = cols, lwd = 2, pch = 16,
                        bty = "n", cex = .85)
}
dev.off()

pdf("figures/benchmark_type1.pdf", width = 12.5, height = 4.0)
par(mfrow = c(1, 3), mar = c(4.6, 4.2, 2.6, 0.8))
bcols <- c("#aed6f1", "#2471a3", "#d6dbdf", "#5d6d7e", "#d5e8a8", "#7d8a2e",
           "#c0392b", "#e67e22")
bnames <- c("DA-plug fixed", "DA-plug anytime", "DA-adj fixed", "DA-adj anytime",
            "Bonf fixed", "Bonf anytime", "SMCS-w (anytime)", "SMCS-w-est (anytime)")
for (sc in scens){
  H <- matrix(NA, 8, length(rhos))
  for (j in seq_along(rhos)){
    o <- grab(sc, 1, rhos[j]); if (is.null(o)) next
    H[, j] <- c(mean(o$res$dap_rej1_fix), mean(o$res$dap_rej1_any),
                mean(o$res$daa_rej1_fix), mean(o$res$daa_rej1_any),
                mean(o$res$bon_rej1_fix), mean(o$res$bon_rej1_any),
                mean(o$res$smcs_rej1_any), mean(o$res$smcs2_rej1_any))
  }
  barplot(H, beside = TRUE, col = bcols,
          ylim = c(0, max(0.3, max(H, na.rm = TRUE) * 1.15)),
          names.arg = sprintf("rho = %.1f", rhos),
          ylab = "Type-I error (reject model 1)", main = paste0(lab[sc], ": null"))
  abline(h = 0.05, lty = 2, lwd = 1.4)
  if (sc == "a") legend("topleft", bnames, fill = bcols, bty = "n", cex = .75)
}
dev.off()
cat("wrote figures/benchmark_setsize.pdf and figures/benchmark_type1.pdf\n")
