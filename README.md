# Conjugator R Package

Quantifying plasmid conjugation just got easier! This package is designed to make the estimation and reporting of plasmid conjugation rates from liquid mating cultures easier, more accurate, and more comparable. 

The two recommended methods to estimate plasmid conjugation rates are **the Simonsen end-point formula** (Simonsen *et al*, J.Gen. Microbiol, 1990), applicable if all strains involved in conjugation grow at the same rate, and **the ASM end-point formula**, which relaxes this assumption ([Huisman *et al*, 2020](https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1)). 

Both of these formulae will cease to be accurate once either the contribution of transconjugants to the overall conjugation becomes substantial, or the recipient dynamics are dominated by conjugation rather than growth. Based on the relative timing of these events, we derived the **critical time** within which the conjugation rate estimates remain valid. 

### Installation 

The quickest way to install Conjugator is to use the devtools package:
```{r}
# Make sure you have devtools; this installs the package from CRAN
install.packages("devtools")
# Then use devtools to install the conjugator package from github
devtools::install_github("JSHuisman/conjugator")
```
Alternatively, you can also download the Conjugator package locally (select the tar.gz file and save it to your computer) and then run the following:
```{r}
# path_to_file indicates where you saved the tar.gz file + full file name and extension
install.packages("path_to_file", repos = NULL, type = "source")
```

### Usage 

This package carries out two main tasks:
- estimate_conj_rate() estimates conjugation rates from experimental data
- estimate_crit_time() reports the critical times

For a walk-through of the most important package functionality, check out the [vignette](./vignettes).

### Contact and citation 
This package was created by Jana S. Huisman, in collaboration with Fabienne Benz, Sarah J.N. Duxbury, J. Arjan G.M. de Visser, Alex R. Hall, Egil A.J. Fischer, and Sebastian Bonhoeffer.

For questions or suggestions concerning the package and the implemented methods please file an issue on this git repository.

If the conjugator app or R package helped you in your work, please cite the manuscript:<br/>
[Jana S. Huisman, Fabienne Benz, Sarah J.N. Duxbury, J. Arjan G.M. de Visser, Alex R. Hall, Egil A.J. Fischer, Sebastian Bonhoeffer. Estimating the rate of plasmid transfer in liquid mating cultures. *BioRxiv* 2020.](https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1)
