###########################################################
## FUNCTIONS FOR GROWTH CURVE FITTING
## Authors: Jana S. Huisman & Sebastian Bonhoeffer
###########################################################

### Logistic Model ########################################
#' Fit input data to a logistic growth function
#'
#' \code{estimate_growth_rate} fits the input OD data to a
#'   logistic growth function.
#' It differs from growthcurver by fitting the logged data and function.
#' It can also subtract an estimated value for blank wells (taking the minimum
#' across all readings), and takes into account measurements only
#' after t_start.
#'
#' @param OD_data Dataframe of OD data. Expects a column for the time (t_col),
#' and will interpret all other columns as OD data for different samples (identified
#' by their column name).
#' @param t_start Earliest possible time point from which to start the fit.
#' @param t_col Name of the time column (defaults to 'time').
#' @param blank_tf Boolean. Should a blank value be subtracted?
#' @param verbose Boolean. Should warnings be returned?
#' @return Dataframe of estimated growth parameters. Contains columns
#' name, carrying_cap, growth_rate, N0, t_min and t_stat.
#' @examples
#' estimate_growth_rate(growth_example, t_start = 1, t_col = 'Time_h')
#' @family growth rate functions
#' @aliases growth growth_rate
#' @export
estimate_growth_rate <- function(OD_data, t_start = 0, t_col = 'time', blank_tf = T, verbose = T){
  .check_growth_input(OD_data, t_col)

  samples = setdiff(colnames(OD_data), c(t_col))
  n_samples = length(samples)

  # we subtract the minimum OD everywhere
  if (blank_tf){
    blank <- min(OD_data[, samples], na.rm = T)
  } else {
    blank <- 0
  }

  growth_fit <- data.frame(name = samples, carrying_cap = NA,
                    growth_rate = NA, N0 = NA, t_min = NA, t_stat = NA)

  # select only trustworthy data points after t_start
  time_sub_OD_data = OD_data[OD_data[, t_col] > t_start, ]

  for (i in 1:n_samples){
    colname <- samples[i]

    # some rows should be non-NA, otherwise .check_growth_input would flag it
    sub_OD_data = time_sub_OD_data[!is.na(time_sub_OD_data[, colname]), ]
    if(nrow(sub_OD_data)<1) {
      if (verbose){
      warning(paste0('Missing data for sample ', colname))
      }
      next
    }

    min_t <- min(sub_OD_data[[t_col]])
    OD_col <- sub_OD_data[[colname]] - blank
    time_col <- sub_OD_data[[t_col]] - min_t
    init_growth_rate <- (log(OD_col[2])-log(OD_col[1]))/(time_col[2]-time_col[1])
    if (is.infinite(init_growth_rate)){
      init_growth_rate <- 20
    }

    # 2 types of errors occur: (i) NAs in the .logist function, (ii) fit errors in nls
    # either results in a try-error and will be caught
    model <- try(stats::nls(log(OD_col) ~ suppressWarnings(log(.logist(time_col, carrying_cap, growth_rate, N0))),
                start = list(carrying_cap = max(OD_col),
                           growth_rate = init_growth_rate,
                           N0 = min(OD_col) ) ), silent = T )
    if('try-error' %in% class(model)) {
      if (verbose){
        warning(paste0('Unable to fit a logistic growth function for ', colname))
      }
      next
    }

    coeffs <- summary(model)$coeff[, 1]
    t_stat_0 = log(coeffs["carrying_cap"]/coeffs["N0"])/coeffs["growth_rate"]

    growth_fit[i,] <- c(samples[i], signif(coeffs, digits = 4),
                        signif(min_t, digits = 4),
                        signif(t_stat_0 + min_t, digits = 4))
  }

  growth_fit[, c('carrying_cap', 'growth_rate', 'N0', 't_min', 't_stat')] <-
    sapply(c('carrying_cap', 'growth_rate', 'N0', 't_min', 't_stat'),
           function(x){as.numeric(growth_fit[,x])})

  return(growth_fit)
}

#' Logistic growth function
#'
#' \code{.logist} computes the logistic growth function.
#'
#' @param t Time vector.
#' @param carrying_cap Carrying capacity.
#' @param growth_rate Growth rate.
#' @param N0 Initial population size.
#' @return Final population size after log growth with the given parameters.
#' @family growth rate functions
.logist <- function(t, carrying_cap, growth_rate, N0){
  result <- carrying_cap/(1 + (carrying_cap - N0)/N0*exp(-growth_rate*t))
  return(result)
}

#' Check growth input data
#' \code{.check_growth_input} performs general checks on the input data.
#'
#' @inheritParams estimate_growth_rate
#'
.check_growth_input <- function(OD_data, t_col){

  if(! t_col %in% colnames(OD_data)){
    stop("The input data does not contain the time column specified.")
  }

  sample_cols = setdiff(colnames(OD_data), t_col)
  if(any(sapply(OD_data[, sample_cols], class) != 'numeric')){
    stop("The input data contains non-numeric columns.")
  }
}

###########################################################
#' Growth fit for plotting
#' \code{get_growth_fit_for_plot} computes a vector of
#' the fit of the growth data, for comparison against the
#' OD data.
#'
#' @inheritParams estimate_growth_rate
#' @param growth_fit Output of estimate_growth_rate for the OD_data;
#' requires columns name, carrying_cap, growth_rate, N0, t_min, t_stat
#' @return Simulated log growth matching the input fits.
#' @examples
#' growth_fit <- estimate_growth_rate(growth_example, t_start = 1, t_col = 'Time_h')
#' get_growth_fit_for_plot(growth_example, t_col = 'Time_h', growth_fit)
#' @family growth rate functions
#' @aliases growth_fit
#' @export
get_growth_fit_for_plot <- function(OD_data, t_col = 'time', growth_fit, blank_tf = T, verbose = T){
  .check_growth_input(OD_data, t_col)
  growth_fit_cols = c("name","carrying_cap", "growth_rate","N0","t_min","t_stat")
  if (!all(growth_fit_cols %in% colnames(growth_fit))&verbose){
    stop('The growth fit data is corrupted: missing necessary columns.')
  }

  samples = growth_fit[['name']]
  if (blank_tf){
    blank <- min(OD_data[, samples], na.rm = T)
  } else {
    blank <- 0
  }

  growth_fit_df <- data.frame()
  for (i in seq_along(samples)){
    if (any(is.na(growth_fit[i, ]))){
      if (verbose){
        warning(paste0('Fit data missing for sample ', samples[i]))
      }
      next
    }
    t_min_ind = which(signif(OD_data[[t_col]], digits = 4) == growth_fit[i,'t_min'])
    sub_OD_data = OD_data[t_min_ind:nrow(OD_data), ]
    short_time_col <- sub_OD_data[[t_col]] - growth_fit[i,'t_min']

    sample_df <- data.frame(ID = samples[i],
                            time = OD_data[[t_col]],
                            OD = c(rep(0, t_min_ind-1),
                                   .logist(short_time_col, growth_fit[i, 'carrying_cap'],
                                        growth_fit[i, 'growth_rate'], growth_fit[i, 'N0']) ) + blank )

    growth_fit_df <- rbind(growth_fit_df, sample_df)
  }

  return(growth_fit_df)
}



