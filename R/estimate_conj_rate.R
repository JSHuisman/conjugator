###########################################################
## FUNCTIONS TO ESTIMATE CONJUGATION RATES
## Authors: Jana S. Huisman & Sebastian Bonhoeffer
###########################################################

#' Estimate plasmid conjugation rates
#'
#' \code{estimate_conj_rate} returns the conjugation
#' rate estimates from data.
#'
#' @param data Dataframe. Required columns vary depending on
#' the method used to estimate the conjugation rates. The following
#' columns are typically required: D.t, R.t, T.t are final population sizes;
#' D.0, R.0 are initial pop sizes;
#' psi.D, psi.R, psi.T are growth rates;
#' t is the time of measurement (assumes start at 0).
#' Optional columns: T.0, psi.max (for the Simonsen method).
#' @param method String or list of strings. Describes method(s) used to
#' estimate the conjugation rate. Takes the following
#' values:
#' \itemize{
#'    \item SM (Simonsen Method): The population growth rate can be supplied with
#'  psi.N or psi.max; otherwise the method will choose
#'  the maximum among psi.D, psi.R, psi.T;
#'    \item ASM (Approximate Simonsen Method): ;
#'    \item TD: \eqn{T/D};
#'    \item T_DR: \eqn{T/DR};
#'    \item T_RT: \eqn{T/(R+T)}: The fraction
#'   of plasmid-carrying recipients;
#'    \item Dionisio: \eqn{T/\sqrt{DR}};
#'    \item Gama: \eqn{log10(T/\sqrt{DR})}.
#' }
#' @param id_cols List of strings. Column names in the data
#' that should be treated as identifiers (will be returned in output).
#' @param verbose Boolean. Should the method return warnings?
#'
#' @return Dataframe with the conjugation rate estimates.
#' Contains 3 columns: sample IDs, estimates, methods.
#' @examples
#' estimate_conj_rate(DRT_example)
#' estimate_conj_rate(DRT_example, method = c('SM', 'T_DR', 'Dionisio'))
#' @family conjugation rate functions
#' @seealso \code{\link{estimate_crit_time}} for estimation of the
#' critical time.
#' @aliases conjugation conj_rate
#' @export
estimate_conj_rate <- function(data, method = c("SM", "ASM"), id_cols = c('ID'), verbose = T){
  data <- as.data.frame(data)
  .check_conj_input(data, method)

  # Loop through methods
  all_results = data.frame()
  for (method_i in method){
    result <- .estimate_conj_per_method(data, method_i, id_cols, verbose)

    all_results <- rbind(all_results, result)
  }

  return(all_results)
}

#' Subroutine to compute individual conjugation rate estimates
#'
#' \code{.estimate_conj_per_method} returns the conjugation
#' rate estimates from data (for a single method).
#'
#' @inheritParams estimate_conj_rate
#'
#' @return Dataframe with conjugation rate estimates.
#' Contains 3 columns: sample IDs, estimates, methods.
.estimate_conj_per_method <- function(data, method = "ASM", id_cols = c('ID'), verbose = T){

  estimate <- switch(method,
                     SM = .estimate_SM(data),
                     ASM = .estimate_ASM(data, verbose),
                     TD = .estimate_TD(data),
                     T_DR = .estimate_T_DR(data),
                     T_RT = .estimate_T_RT(data),
                     Dionisio = .estimate_Dionisio(data),
                     Gama = .estimate_Gama(data)
  )

  if (is.null(estimate)) stop("Unknown method. Type ?estimate_conj_rate to check which methods are allowed.")

  # combine output
  if(length(intersect(id_cols, colnames(data))) < 1) {
    id_cols = 'ID'
    data[,'ID'] <- as.factor(1:nrow(data))
  }
  result <- data.frame(lapply(data[id_cols], factor), estimate = as.numeric(estimate),
                       method = as.character(method))

  return(result)
}

###########################################################
#' Check conjugation rate estimation input
#' @inheritParams estimate_conj_rate
.check_conj_input <- function(data, method){
  target_cols = c("D.t", "R.t", "T.t", "D.0", "R.0",
                   "psi.D", "psi.R", "psi.T", "t")
  target_cols_present = intersect(target_cols, colnames(data))

  if(length(target_cols_present)<=1){
    stop("The input data does not contain enough data columns.")
  }

  if(any(sapply(data[, target_cols_present], class) != 'numeric')){
    stop("The input data contains non-numeric columns.")
  }

}

#' Check if necessary columns are present
#' @inheritParams estimate_conj_rate
#' @param necessary_cols list of columns to check the presence of.
.check_cols_present <- function(data, necessary_cols){
  if(!all(necessary_cols %in% colnames(data))){
    missing_cols = setdiff(necessary_cols, colnames(data))
    stop(paste0("Missing strictly necessary columns in the input data.
                The following columns are missing: ", paste(missing_cols, collapse = ", ")) )
  }
}

###########################################################

#' Simonsen conjugation rate estimate
#'
#' \code{.estimate_SM} returns the Simonsen
#'  formula from data.
#'  The population growth rate can be supplied with
#'  psi.N or psi.max; otherwise the method will choose
#'  the maximum among psi.D, psi.R, psi.T.
#'
#' @inheritParams estimate_conj_rate
#' @return estimate vector
#' @aliases gamma_max
.estimate_SM <- function(data){
  necessary_cols = c("D.t", "R.t", "T.t", "D.0", "R.0")
  .check_cols_present(data, necessary_cols)

  growth_cols = c("psi.R", "psi.D", "psi.T", "psi.N", "psi.max")
  if(!any(growth_cols %in% colnames(data))){
    stop("Please supply at least one growth rate column in the input data.
         Accepted names: psi.R, psi.D, psi.T, psi.N, or psi.max")
  }

  # add starting transconjugant population size
  if (!"T.0" %in% colnames(data)){
    data$T.0 = rep(0, dim(data)[1])
  } else if (any(is.na(data[,"T.0"]))){
    data$T.0[is.na(data$T.0)] = 0
  }

  # add common population growth rate
  if ("psi.max" %in% colnames(data)){
    #don't do anything
  } else if ("psi.N" %in% colnames(data)){
    data$psi.max = data$psi.N
  } else {
    present_growth_rates = intersect(colnames(data), c("psi.R", "psi.D", "psi.T"))
    data$psi.max = as.numeric(apply(data, 1, function(x) max(x[present_growth_rates])))
  }

  estimate <- data[, 'psi.max']*log(1+data[,'T.t']*(
    data[, 'D.t']+data[, 'R.t']+data[, 'T.t'])/(data[, 'R.t']*data[, 'D.t']))/
    (data[, 'D.t']+data[, 'R.t']+data[, 'T.t']-data[, 'D.0']-data[, 'R.0']-data[, 'T.0'])
  return(estimate)
}


#' Estimate ASM conjugation rates
#'
#' \code{.estimate_ASM} returns the ASM
#'  formula from data.
#'
#' Details can be found in Huisman et al., 2020
#' \url{https://www.biorxiv.org/content/10.1101/2020.03.09.980862v1}
#'
#' @inheritParams estimate_conj_rate
#' @return estimate vector
#' @aliases gammaD
.estimate_ASM <- function(data, verbose = T){
  necessary_cols = c("psi.D", "psi.R", "psi.T", "D.t", "R.t", "T.t", "D.0", "R.0", "t")
  .check_cols_present(data, necessary_cols)

  estimate_func <- function(i){
    data_row <- data[i, ]

    if (data_row[['psi.D']] + data_row[['psi.R']] == data_row[['psi.T']]){
      result <- data_row[['T.t']]/
        (data_row[['D.0']]*data_row[['R.0']] * exp(data_row[['psi.T']] *data_row[['t']])*data_row[['t']])
    } else {
      result <- (data_row[['psi.D']] + data_row[['psi.R']] - data_row[['psi.T']])*data_row[['T.t']]/
        (data_row[['D.t']]*data_row[['R.t']] - data_row[['D.0']]*data_row[['R.0']] *
           exp(data_row[['psi.T']] *data_row[['t']]))
    }

    return(result)
  }

  estimate <- sapply(1:nrow(data), estimate_func)

  if (any(estimate < 0, na.rm = T) & verbose){
    warning('Negative conjugation rates estimated by the ASM.
            The experiment likely exceeded the critical time.')
  }
  estimate[estimate < 0] <- NA

  return(estimate)
}

#' Estimate the transconjugant fraction
#'
#' \code{.estimate_T_RT} returns the fraction
#'   of plasmid-carrying recipients.
#'
#' @inheritParams estimate_conj_rate
#' @return estimate vector
#' @aliases transconjugant_frac
.estimate_T_RT <- function(data){
  necessary_cols = c("R.t", "T.t")
  .check_cols_present(data, necessary_cols)

  estimate <- data[, 'T.t']/(data[, 'R.t'] + data[, 'T.t'])
  return(estimate)
}

#' Estimate T/D from data
#'
#' \code{.estimate_TD} returns the T/D estimate.
#' @return estimate vector
#' @inheritParams estimate_conj_rate
.estimate_TD <- function(data){
  necessary_cols = c("D.t", "T.t")
  .check_cols_present(data, necessary_cols)

  estimate <- data[, 'T.t']/data[, 'D.t']
  return(estimate)
}

#' Estimate T/DR from data
#'
#' \code{.estimate_T_DR} returns the T/DR estimate.
#' @return estimate vector
#' @inheritParams estimate_conj_rate
.estimate_T_DR <- function(data){
  necessary_cols = c("D.t", "R.t", "T.t")
  .check_cols_present(data, necessary_cols)

  estimate <- data[, 'T.t']/(data[, 'D.t'] * data[, 'R.t'])
  return(estimate)
}

#' Dionisio estimate from data
#'
#' \code{.estimate_Dionisio} returns the
#' Dionisio estimate: T/sqrt(DR).
#'
#' Dionisio F, Matic I, Radman M, Rodrigues OR, Taddei F.
#' Plasmids spread very fast in heterogeneous bacterial communities.
#' Genetics 2002 Dec;162(4):1525–32.
#'
#' @inheritParams estimate_conj_rate
#' @return estimate vector
.estimate_Dionisio <- function(data){
  necessary_cols = c("D.t", "R.t", "T.t")
  .check_cols_present(data, necessary_cols)

  estimate <- data[, 'T.t']/sqrt(data[, 'R.t']*data[, 'D.t'])
  return(estimate)
}

#' Gama estimate from data
#'
#' \code{.estimate_Gama} returns the Gama
#' estimate: log10(T/sqrt(DR)).
#' This is the log10 of the Dionisio estimate.
#'
#' João Alves Gama, Rita Zilhão, and Francisco Dionisio.
#' Multiple plasmid interference - Pledging allegiance to my enemy’s enemy.
#' Plasmid, 93(August):17–23, 2017.
#'
#' @inheritParams estimate_conj_rate
#' @return estimate vector
.estimate_Gama <- function(data){
  necessary_cols = c("D.t", "R.t", "T.t")
  .check_cols_present(data, necessary_cols)

  estimate <- log10(data[, 'T.t']/sqrt(data[, 'R.t']*data[, 'D.t']))
  return(estimate)
}



