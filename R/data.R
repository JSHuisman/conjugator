#' Example conjugation rates (DRT)
#'
#' An example DRT dataset containing made up conjugation rates
#'
#' @format A data frame with 3 rows and 10 variables:
#' \describe{
#'   \item{ID}{sample ID}
#'   \item{psi.D}{growth rate of the donors, per hour}
#'   \item{psi.R}{growth rate of the recipients, per hour}
#'   \item{psi.T}{growth rate of the transconjugants, per hour}
#'   \item{D.0}{initial population size of the donors, in CFU/mL}
#'   \item{R.0}{initial population size of the recipients, in CFU/mL}
#'   \item{D.t}{final population size of the donors, in CFU/mL}
#'   \item{R.t}{final population size of the recipients, in CFU/mL}
#'   \item{T.t}{final population size of the transconjugants, in CFU/mL}
#'   \item{t}{Measurement timepoint, in hours}
#' }
"DRT_example"

#' Example conjugation rates (TRT)
#'
#' An example TRT dataset containing made up conjugation rates
#'
#' @format A data frame with 3 rows and 10 variables:
#' \describe{
#'   \item{ID}{sample ID}
#'   \item{psi.D}{growth rate of the donors, per hour}
#'   \item{psi.R}{growth rate of the recipients, per hour}
#'   \item{psi.T}{growth rate of the transconjugants, per hour}
#'   \item{D.0}{initial population size of the donors, in CFU/mL}
#'   \item{R.0}{initial population size of the recipients, in CFU/mL}
#'   \item{D.t}{final population size of the donors, in CFU/mL}
#'   \item{R.t}{final population size of the recipients, in CFU/mL}
#'   \item{T.t}{final population size of the transconjugants, in CFU/mL}
#'   \item{t}{Measurement timepoint, in hours}
#' }
"TRT_example"

#' Example OD measurements
#'
#' @format A data frame with 3 rows and 10 variables:
#' \describe{
#'   \item{Time_h}{time in hours}
#'   \item{Donor_OD}{OD of the donors, over 12 hours}
#'   \item{Recipient_OD}{OD of the recipients, over 12 hours}
#' }
"growth_example"
