################################################################################
## ch1_benchmark_all.R -- runs ch1_benchmark.R over all 18 cells
## (3 scenarios x 3 correlations x {alternative, null}).
##
## Usage:  Rscript R/ch1_benchmark_all.R [nsims]   (default 500)
################################################################################
args <- commandArgs(trailingOnly = TRUE)
ns   <- if (length(args) >= 1) args[1] else "500"
for (scen in c("a","b","c"))
  for (rho in c("0","0.4","0.8"))
    for (nul in c("0","1"))
      system2("Rscript", c("R/ch1_benchmark.R", scen, rho, nul, ns))
