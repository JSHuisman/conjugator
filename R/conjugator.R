#' conjugator: A package to compute conjugation rate estimates
#'
#' The conjugator package provides the functions for
#' conjugation rate estimation, the assessment of critical times,
#' and basic growth rate estimation.
#' Details about these methods can be found in Huisman et al., 2020
#' \url{https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1}
#'
#'This functionality can also be explored in a Shiny app,
#'available at \url{https://ibz-shiny.ethz.ch/jhuisman/conjugator/}.
#'
#' @section Conjugation rate estimation:
#' The estimate_conj_rate function allows the user to compute a variety
#' of conjugation rate estimates from experimental data.
#' Currently implemented are the Simonsen and Approximate Simonsen endpoint formulae,
#' T/D, T/DR, T/(T+R), as well as T/sqrt(DR) and log(T/sqrt(DR))
#' (called Dionisio and Gama, respectively).
#'
#' @section Critical time estimation:
#' The applicability of the formulae to assess conjugation rates breaks down
#' once transconjugants become the dominant source of conjugation event, and
#' once the recipient population is depleted by conjugation. Each of these
#' processes determines a "critical time". We offer functions to estimate these
#' critical times from experimental data.
#'
#' @section Growth rate estimation:
#' The package also contains some basic functions to estimate population
#' growth rates from OD data. This is for convenience, and only in
#' service of the conjugation rate estimation. We take no responsibility
#' for the correctness of the growth results, and generally recommend
#' using a specialized package (e.g. growthcurver) for this purpose.
#'
#' @docType package
#' @name conjugator
#'
#' @importFrom stats na.omit
#' @importFrom stats uniroot
NULL
