# Conjugator R Package

Quantifying plasmid conjugation just got easier! This package is designed to make the estimation and reporting of plasmid conjugation rates from liquid mating cultures easier, more accurate, and more comparable. This functionality can also be explored in a [Shiny app](https://ibz-shiny.ethz.ch/jhuisman/conjugator/). 

The two recommended methods to estimate plasmid conjugation rates are **the Simonsen end-point formula** (Simonsen *et al*, J.Gen. Microbiol, 1990), applicable if all strains involved in conjugation grow at the same rate, and **the ASM end-point formula**, which relaxes this assumption (Huisman *et al*, 2020). 

Both of these formulae will cease to be accurate once either the contribution of transconjugants to the overall conjugation becomes substantial, or the recipient dynamics are dominated by conjugation rather than growth. Based on the relative timing of these events, we derived the **critical time** within which the conjugation rate estimates remain valid. 

More details about these methods can be found in [our paper](https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1).

### Installation 

The quickest way to install Conjugator is to use the devtools package:
```{r}
# Make sure you have devtools; this installs the package from CRAN
install.packages("devtools")
# Then use devtools to install the conjugator package from github
library(devtools)
install_github("JSHuisman/conjugator", build_vignettes = TRUE)
```
The package devtools generally requires the user to have a working development environment, which includes [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, and Xcode on Mac.

Alternatively, you can also download the Conjugator package locally (select the tar.gz file and save it to your computer) and then run the following:
```{r}
# path_to_file indicates where you saved the tar.gz file + full file name and extension
install.packages("path_to_file", repos = NULL, type = "source")
```
Note: the Conjugator package is not yet listed on CRAN, so the default install.packages in R will not work!

### Usage 

This package has two main functions:
(i) Estimate conjugation rates from experimental data
```{r}
# the DRT_example is included with the package
estimate_conj_rate(DRT_example, "ASM")
```
(ii) Estimate the critical times
```{r}
# DRT and TRT refer to the first and second conjugation experiment
# TRT_example is also included in the package
estimate_crit_time(DRT_example, TRT_example, tol_factor = 10)
```

For a walk-through of the most important package functionality, check out the [vignette](https://jshuisman.github.io/conjugator/articles/conjugator.html).

### Contact and citation 
This package was created by Jana S. Huisman, in collaboration with Fabienne Benz, Sarah J.N. Duxbury, J. Arjan G.M. de Visser, Alex R. Hall, Egil A.J. Fischer, and Sebastian Bonhoeffer.

For questions or suggestions concerning the package and the implemented methods please file an issue on this git repository.

If the conjugator app or R package helped you in your work, please cite the manuscript:<br/>
[Jana S. Huisman, Fabienne Benz, Sarah J.N. Duxbury, J. Arjan G.M. de Visser, Alex R. Hall, Egil A.J. Fischer, Sebastian Bonhoeffer. Estimating the rate of plasmid transfer in liquid mating cultures. *BioRxiv* 2020.](https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1)
