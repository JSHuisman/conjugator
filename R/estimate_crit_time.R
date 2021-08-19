###########################################################
## FUNCTIONS TO ESTIMATE CRITICAL TIMES
## Authors: Jana S. Huisman
###########################################################

#' Estimate critical times
#'
#' \code{estimate_crit_time} returns the
#'  conjugation rates and critical times corresponding
#'  to a conjugation experiment.
#'  For the ASM to return a proper conjugation estimate
#'  the time point of measurement must be smaller than
#'  min(tcrit1, tcrit2, tcrit3).
#'
#' @param DRT Dataframe. Data corresponding to the first
#' conjugation experiment.
#' Necessary columns: D.0, R.0 are initial pop sizes;
#' psi.D, psi.R, psi.T are growth rates.
#' If gamma.D is not among the columns,
#' the function will attempt to estimate gamma.D from
#' the data given. This requires D.t, R.t, T.t
#' (final population sizes) and t (the time of measurement).
#' If gamma.T is among the columns, the input under TRT will
#' be ignored.
#' @param TRT Dataframe. Data corresponding to the second
#' conjugation experiment. Same columns needed as for DRT.
#' @param tol_factor Double. The imprecision factor: this
#' factor indicates how close the system is allowed to
#' approach a state where it violates the ASM assumptions
#' (smaller values allow greater violation).
#' @param id_cols List of strings. Specifies the column
#' names that should be treated as identifiers (will be returned in output).
#' @param verbose Boolean. Should warnings be returned?
#' @return Dataframe of critical time estimates.
#' @examples
#' estimate_crit_time(DRT_example, TRT_example)
#' estimate_crit_time(DRT_example, TRT_example, id_cols = c('ID', 't'))
#' @family critical time functions
#' @aliases critical_time crit_time
#' @export
estimate_crit_time <- function(DRT, TRT = NULL, tol_factor = 10, id_cols = c('ID'), verbose = T){
  if(!is.data.frame(DRT)) stop("The DRT input must be a dataframe.")

  ## Check ID_cols of DRT
  if(length(intersect(id_cols, colnames(DRT))) < 1) {
    DRT[,'ID'] <- as.factor(1:nrow(DRT))
    id_cols <- 'ID'
  }

  data <- .get_crit_time_input(DRT, TRT, id_cols, verbose)

  compute_tcrit <- .check_crit_time_input(data, verbose)

  t_crit_df <- data.frame(data[id_cols], tcrit1 = rep(NA, nrow(data)),
                          tcrit2 = rep(NA, nrow(data)), tcrit3 = rep(NA, nrow(data)),
                          min_tcrit = rep(NA, nrow(data)))

  if (compute_tcrit[1]){
    t_crit_df[, 'tcrit1'] = sapply(1:nrow(data), function(i){.estimate_tcrit1(data[i, ], tol_factor = 10)})
  }

  if (compute_tcrit[2]){
    t_crit_df[,'tcrit2'] = sapply(1:nrow(data), function(i){.estimate_tcrit2(data[i, ], tol_factor = 10)})
  }

  if (compute_tcrit[3]){
    t_crit_df[,'tcrit3'] = sapply(1:nrow(data), function(i){.estimate_tcrit3(data[i, ], tol_factor = 10)})
  }

  # add the minimum
  t_crit_df[,'min_tcrit'] = sapply(1:nrow(t_crit_df), function(x){min(t_crit_df[x,'tcrit1'],
                                                                t_crit_df[x,'tcrit2'],
                                                                t_crit_df[x,'tcrit3'], na.rm = T)})

  # return full estimate
  if (! 'gamma.D' %in% colnames(data)) data[, 'gamma.D'] <- NA
  if (! 'gamma.T' %in% colnames(data)) data[, 'gamma.T'] <- NA
  full_estimate <- merge(data[, c(id_cols, 'gamma.D', 'gamma.T')], as.data.frame(t_crit_df), by = id_cols)

  return(full_estimate)
}


#' Estimate tcrit1
#' @inheritParams estimate_crit_time
#' @param data Dataframe combining both DRT and TRT data.
#' @return Estimates for tcrit1
.estimate_tcrit1 <- function(data, tol_factor = 10){
  if (data[['psi.D']] + data[['psi.R']] > data[['psi.T']]){
    tcrit1 <- log((data[['psi.D']] + data[['psi.R']])/
                    (tol_factor * data[['gamma.T']] * data[['R.0']] ))/data[['psi.R']]
  } else {
    tcrit1 <- log((data[['psi.T']])/
                    (tol_factor * data[['gamma.T']] * data[['R.0']] ))/data[['psi.R']]
  }

  tcrit1 <- signif(tcrit1, digits = 3)
  return(tcrit1)
}

#' Estimate tcrit2
#' @inheritParams estimate_crit_time
#' @param data Dataframe combining both DRT and TRT data.
#' @return Estimate for tcrit2
.estimate_tcrit2 <- function(data, tol_factor = 10){
  tcrit2 <- log(data[['psi.R']]/(tol_factor * data[['gamma.D']] * data[['D.0']]))/
    data[['psi.D']]
  tcrit2 <- signif(tcrit2, digits = 3)
  return(tcrit2)
}

#' Estimate tcrit3
#' @inheritParams estimate_crit_time
#' @param data Dataframe combining both DRT and TRT data.
#' @return Estimates for tcrit3
.estimate_tcrit3 <- function(data, tol_factor = 10){
  if (data[['psi.D']] + data[['psi.R']] > data[['psi.T']]){
    tcrit3 <- log(data[['psi.R']] *
                    (data[['psi.D']] + data[['psi.R']] - data[['psi.T']])/
                    (tol_factor*data[['gamma.D']] * data[['gamma.T']] * data[['D.0']] * data[['R.0']]))/
      (data[['psi.D']] + data[['psi.R']])
  } else if (data[['psi.D']] + data[['psi.R']] == data[['psi.T']]){
    # tcrit3 <- log(data[['psi.R']] /
    #                 (tol_factor*data[['gamma.D']] * data[['gamma.T']] * data[['D.0']] * data[['R.0']]))/
    #   (data[['psi.D']] + data[['psi.R']])
    tcrit3 <- uniroot(.t3_rootfunc, c(-100, 100), data, tol_factor)$root
  } else {
    tcrit3 <- log(data[['psi.R']] *
                    (data[['psi.T']] - data[['psi.D']] - data[['psi.R']])/
                    (tol_factor*data[['gamma.D']] * data[['gamma.T']] * data[['D.0']] * data[['R.0']]))/
      data[['psi.T']]
  }

  tcrit3 <- signif(tcrit3, digits = 3)
  return(tcrit3)
}

# For use in the rootfinder for psiD + psiR = psiT
.t3_rootfunc <- function(t, data, tol_factor = 10){
  return_val <- data[['psi.R']] - tol_factor*data[['gamma.D']] * data[['gamma.T']] *
    data[['D.0']] * data[['R.0']]*exp(data[['psi.T']]*t)*t
  return(return_val)
}


###########################################################
#' Shape data to estimate critical times
#'
#' \code{.get_crit_time_input} returns the dataframe needed
#' to estimate the critical times corresponding
#' to a conjugation experiment.
#'
#' Note: the number of columns in the output can vary,
#' depending on whether gamma.T can be estimated or not.
#'
#' @inheritParams estimate_crit_time
#'
.get_crit_time_input <- function(DRT, TRT, id_cols = c('ID'), verbose = T){

  ## Compute gamma.D if it does not exist already
  if (! "gamma.D" %in% colnames(DRT)){
    first_result <- try(estimate_conj_rate(DRT, method = c("ASM")), silent = T)
    if ('try-error' %in% class(first_result)){
      if (verbose){
        warning("Unable to estimate gamma.D from the DRT input. Please add gamma.D as column directly,
           or supply the necessary input columns to estimate gamma.D .\n")
      }
    } else {

      if (any(is.na(first_result[2]))){
        if (verbose){
          warning(paste0('The ASM method generated NAs in row ', which(is.na(first_result[2])),
                         ', using SM to determine gamma.D instead.\n'))
        }
        first_result <- estimate_conj_rate(DRT, method = c("SM"))

      }

      DRT[, 'gamma.D'] <- first_result[2]
    }
  }

  ## Create minimal crit_data
  target_cols <- c(id_cols, "psi.D", "psi.R", "psi.T", "D.0", "R.0", "gamma.D", "gamma.T")
  present_targets <- intersect(colnames(DRT), target_cols)
  crit_data <- DRT[, present_targets]

  ## If gamma.T is present, we're done
  if ("gamma.T" %in% colnames(crit_data)) return(crit_data)

  ## Add gamma.T if it does not exist already
  if (! "gamma.T" %in% colnames(TRT)){
    second_result <- try(estimate_conj_rate(TRT, method = c("ASM")), silent = T)

    if ('try-error' %in% class(second_result)){
      if(verbose){
      warning("Unable to estimate gamma.T from the TRT input. Please add gamma.T as column directly,
           or supply the necessary input columns to estimate gamma.T .\n")
      }
      return(crit_data)
    }

    if (any(is.na(second_result[2]))){
      if (verbose){
      warning(paste0('The ASM method generated NAs in row ', which(is.na(second_result[2])),
                     ', using SM to determine gamma.T instead.\n'))
        }
      second_result <- estimate_conj_rate(TRT, method = c("SM"))
    }

    TRT[, 'gamma.T'] <- second_result[2]
  }

  ## Check ID cols of TRT
  if(length(intersect(id_cols, colnames(TRT))) < 1) {
    TRT[,'ID'] <- as.factor(1:nrow(TRT))
    id_cols <- 'ID'
  }
  if(!all(intersect(id_cols, colnames(TRT)) == intersect(id_cols, colnames(DRT))) ){
    stop('The id_cols of DRT and TRT do not match.')
  }
  if (nrow(merge(DRT, TRT, by = id_cols)) == 0){
    stop('No matching IDs could be found between DRT and TRT.')
  }
  if (!identical(DRT[id_cols], TRT[id_cols])&verbose){
    warning('The entries in the id_cols of DRT and TRT are not equal. Maximum matching set will be used.')
  }

  ## Combining both results
  crit_data <- merge(crit_data, TRT[, c(id_cols, 'gamma.T')], by = id_cols)

  return(crit_data)
}
###########################################################
#' Checks input to estimate critical times
#' @inheritParams estimate_crit_time
#' @param data Dataframe combining both DRT and TRT data.
.check_crit_time_input <- function(data, verbose = T){
  target_cols_t1 = c("R.0", "psi.R",
                    "psi.D", "psi.T","gamma.T")
  target_cols_t2 = c("D.0", "psi.R",
                     "psi.D","gamma.D")
  target_cols_t3 = c("D.0", "R.0", "psi.R",
                     "psi.D", "psi.T","gamma.D","gamma.T")
  compute_tcrit = c(T, T, T)

  ## Error if no results are possible
  if (!(all(target_cols_t1 %in% colnames(data))|
        all(target_cols_t2 %in% colnames(data)) ) ){
    stop("The input data does not contain the data columns needed to compute at least one critical time.")
  }

  ## Warnings for incomplete results
  if(!all(target_cols_t1 %in% colnames(data)) ){
    if (verbose){
    warning("The input data does not contain the data columns needed to compute tcrit1.")
    }
    compute_tcrit[1] = F
  }
  if (!all(target_cols_t2 %in% colnames(data)) ){
    if (verbose){
      warning("The input data does not contain the data columns needed to compute tcrit2.")
    }
    compute_tcrit[2] = F
  }
  if (!all(target_cols_t3 %in% colnames(data)) ){
    if (verbose){
      warning("The input data does not contain the data columns needed to compute tcrit3.")
    }
    compute_tcrit[3] = F
  }

  ## Require numeric columns
  target_cols_present = intersect(colnames(data), target_cols_t3)
  if(any(sapply(data[, target_cols_present], class) != 'numeric')){
    stop("The input data contains non-numeric columns")
  }

  return(compute_tcrit)
}

###########################################################
#' Scan critical times
#'
#' \code{scan_crit_time} returns the
#' conjugation rates and critical times corresponding to
#' a DRT conjugation experiment, assuming different gamma.T values.
#'  For the ASM to return a proper conjugation estimate
#'  the time point of measurement must be smaller than
#'  min(tcrit1, tcrit2, tcrit3).
#'
#' @inheritParams estimate_crit_time
#' @param mult_seq Vector. Describes the multiplication factors
#' to take into account for gamma.T.
#' @return Dataframe of estimated critical times for the input DRT
#' data and a range of gamma.T values.
#' @examples
#' scan_crit_time(DRT_example)
#' scan_crit_time(DRT_example, id_cols = c('ID', 't'), mult_seq = 10**seq(0, 3, 0.1))
#' @family critical time functions
#' @seealso \code{\link{estimate_conj_rate}} for estimation of the
#' conjugation rates.
#' @aliases scan_crit time_scan crittime_scan
#' @export
scan_crit_time <- function(DRT, tol_factor = 10, id_cols = c('ID'), mult_seq = 10**seq(-2, 5, 1)){

  if (!'gamma.D' %in% colnames(DRT)){
    ASM_estimate = try(.estimate_ASM(DRT))
    if ('try-error' %in% class(ASM_estimate)){
      stop("No gamma.D estimate was present in the input data and none could be estimated.")
    } else {
      DRT[,'gamma.D'] = ASM_estimate
    }
  }

  add_gammaT <- function(mult_factor){cbind(DRT, 'gamma.T' = mult_factor * DRT[, 'gamma.D'],
                                            'mult_factor' = mult_factor)}
  full_data <- do.call("rbind", lapply(mult_seq, add_gammaT))

  tcrit_df <- estimate_crit_time(full_data, TRT = NULL, tol_factor, id_cols = c(id_cols, 'mult_factor'))

  return(tcrit_df)
}


