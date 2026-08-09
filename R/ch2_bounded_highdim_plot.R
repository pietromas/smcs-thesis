################################################################################
## ch2_bounded_highdim_plot.R -- Figure 2.2
## Rejection-time CDFs, panels grouped by contrast.
## Run after ch2_bounded_highdim.R.
################################################################################
A <- readRDS("results/mvhd_all.rds"); res <- A$res; n <- A$n

panels <- c("pair_compl",  "pair_redund",
            "six_compl",   "three_redund",
            "nine_graded", "four_mixed",
            "single_corr", "single_iso")
panels <- panels[panels %in% names(res)]

titles <- c(pair_compl   = "Two violations, complementary",
            pair_redund  = "Two violations, redundant",
            six_compl    = "Six violations, complementary",
            three_redund = "Three violations, redundant",
            nine_graded  = "Nine graded violations, mixed",
            four_mixed   = "Four violations, mixed",
            single_corr  = "Single violation, correlated",
            single_iso   = "Single violation, isolated")

cols <- c(AVG = "#2471a3", MV = "#c0392b", AVGor = "#2471a3", MVor = "#c0392b")
ltys <- c(AVG = 1, MV = 1, AVGor = 3, MVor = 3)
grid <- 0:n
nr   <- ceiling(length(panels) / 2)

pdf("figures/mvhd_rejection_curves.pdf", width = 9, height = 2.7 * nr)
par(mfrow = c(nr, 2), mar = c(4, 4, 2.5, 1))
for (nm in panels){
  TAU <- res[[nm]]$tau; N <- nrow(TAU)
  plot(NA, xlim = c(0, n), ylim = c(0, 1), xlab = "t",
       ylab = "rejection rate by t", main = titles[nm])
  for (m in c("AVG", "MV", "AVGor", "MVor")){
    y <- sapply(grid, function(x) sum(TAU[, m] <= x, na.rm = TRUE)) / N
    lines(grid, y, col = cols[m], lty = ltys[m], lwd = 2)
  }
  if (nm == panels[1])
    legend("bottomright", bty = "n", lwd = 2,
           col = cols[c("AVG","MV","AVGor","MVor")],
           lty = ltys[c("AVG","MV","AVGor","MVor")],
           legend = c("Average (aGRAPA)", "Multivariate (plug-in)",
                      "Average (oracle)", "Multivariate (oracle)"))
}
dev.off()

cat("\nscenario        rejAVG rejMV   tauAVG   tauMV tauAVGor tauMVor  ratio\n")
for (nm in names(res)){
  TAU <- res[[nm]]$tau; a4 <- complete.cases(TAU)
  mt  <- if (any(a4)) colMeans(TAU[a4, , drop = FALSE]) else rep(NA, 4)
  cat(sprintf("%-15s %6.3f %6.3f %8.1f %7.1f %8.1f %7.1f %6.2f\n",
              nm, mean(!is.na(TAU[,1])), mean(!is.na(TAU[,2])),
              mt[1], mt[2], mt[3], mt[4], mt[1]/mt[2]))
}
cat("\nwrote figures/mvhd_rejection_curves.pdf\n")
