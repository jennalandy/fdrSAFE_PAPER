# fdrSAFE: Selective Aggregation for Local False Discovery Rate (fdr) Estimation

For the `fdrSAFE` R software package, see the [jennalandy/fdrSAFE](https://github.com/jennalandy/fdrSAFE) repository.

This repository contains the code to replicate all results reported in the paper [*fdrSAFE: Selective Aggregation for Local False Discovery Rate Estimation*](https://arxiv.org/abs/2401.12865). See details in the Simulation Study and Experimental Application sections of our paper.

Throughout this repository, we use [quarto](https://quarto.org/) documents which can be edited and run with RStudio, Jupyter Lab, or Visual Studio Code.

### Simulation Studies

Our simulation studies are in R scripts. When each script is run, it will log progress and results in a new sub-directory. Scripts assume you are in the `simulation_studies` directory. We use [`sink`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/sink) to log progress; if you terminate a test early you will need to run `sink()` for output to show up in the console again. Note that these scripts take many hours to run.

See [Simulation Study README](simulations/README.md) for more details.

### Experimental Application

Our experimental application relies on the Platinum Spike dataset<sup>1</sup>.

See [Experimental README](experimental/README.md) for more details.

### References

[1] Q. Zhu, J.C. Miecznikowski, and M.S. Halfon. Preferred analysis methods for affymetrix genechips. II. an expanded, balanced, wholly-defined spike-in dataset. *BMC Bioinformatics*, 11:285, 2010. doi:https://doi.org/10.1186/1471-2105-11-285.
