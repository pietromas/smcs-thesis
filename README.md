# Sequential Model Confidence Sets: adaptive weights and multivariate betting

Code to reproduce the simulation studies and figures of the MSc thesis
*Sequential Model Confidence Sets and New Horizons*, ETH Zurich.

## Requirements

R (>= 4.0) and the packages:

    install.packages(c("quadprog", "scoringRules", "MASS",
                       "future.apply", "ggplot2"))

The Covid-19 case study additionally needs several non-CRAN packages; see the
header of `R/load_covid_data.R`.

## Contents

| Script | Produces |
|---|---|
| `R/ch1_simulation.R` | Figure 1.1 |
| `R/ch1_benchmark.R` | one cell of the benchmark |
| `R/ch1_benchmark_all.R` | runs all 18 benchmark cells |
| `R/ch1_benchmark_plot.R` | Figures 1.2 and 1.3 |
| `R/load_covid_data.R` | downloads the Covid-19 Hub data into `data/` |
| `R/ch1_covid.R` | Figure 1.4 |
| `R/ch2_bounded_p3.R` | Table 2.1 |
| `R/ch2_bounded_highdim.R` | Tables 2.2 and 2.3 |
| `R/ch2_bounded_highdim_plot.R` | Figure 2.2 |
| `R/ch2_theory_check.R` | Table 2.4 |
| `R/ch2_gauss_highdim.R` | Tables 2.6 to 2.8 |
| `R/ch2_gauss_highdim_plot.R` | Figure 2.3 |

## Running

From the repository root. Simulation scripts write `.rds` files to `results/`
and print one summary line per scenario. Plotting scripts read those files and
write PDFs to `figures/`.

    Rscript R/ch2_bounded_highdim.R
    Rscript R/ch2_bounded_highdim_plot.R
    Rscript R/ch2_theory_check.R

    Rscript R/ch1_benchmark_all.R
    Rscript R/ch1_benchmark_plot.R

Each script takes its parameters as command-line arguments and defaults to the
values used in the thesis. Seeds are fixed, so the reported numbers reproduce.

## Data

Section 1.6 uses the public Covid-19 Forecast Hub data, which are not
redistributed here. Run

    Rscript R/load_covid_data.R

once to download them into `data/`. The script header lists the non-CRAN
packages that must be installed first.

## License

MIT, see `LICENSE`.
