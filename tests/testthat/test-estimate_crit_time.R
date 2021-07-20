###########################################################
## Unit tests for the critical time estimation
## Authors: Jana S. Huisman
###########################################################
context("Critical Time Estimation")
library(conjugator)

###########################################################
test_that("example input data returns right result", {
  example_data <- data.frame('psi.D' = 1, 'psi.R' = 1.5, 'psi.T' = 1.5,
                             'D.0' = 1e5, 'R.0' = 1e5,
                             'gamma.D' = 1e-11, 'gamma.T' = 1e-9)

  reference_result <- data.frame('ID' = factor(1), 'gamma.D' = 1e-11, 'gamma.T' = 1e-9,
                                 'tcrit1' = 5.22, 'tcrit2' = 11.9,
                                 'tcrit3' = 8.45, 'min_tcrit' = 5.22)

  expect_equal(estimate_crit_time(example_data, tol_factor = 10), reference_result)
  expect_equal(suppressWarnings(estimate_crit_time(example_data[1:6], tol_factor = 10)),
               cbind(reference_result[1:2], 'gamma.T' = NA, 'tcrit1' = NA,
                     reference_result['tcrit2'], 'tcrit3' = NA,
                     'min_tcrit' = unname(reference_result['tcrit2'])) )
  expect_equal(suppressWarnings(estimate_crit_time(example_data[c(1:5, 7)], tol_factor = 10)),
               cbind(reference_result['ID'], 'gamma.D' = NA, reference_result[3:4], 'tcrit2' = NA,
                     'tcrit3' = NA, 'min_tcrit' = unname(reference_result['tcrit1'])) )
})

test_that("example input DRT twice returns right result", {
  reference_result <- data.frame('ID' = factor(c("A1", "A2", "A3")),
                                 'gamma.D' = c(9.004904e-10, 9.004904e-10, 8.008783e-10),
                                 'gamma.T' = c(9.004904e-10, 9.004904e-10, 8.008783e-10),
                                 'tcrit1' = c(9.9, 9.9, 10.0),
                                 'tcrit2' = c(11.6, 11.6, 11.8),
                                 'tcrit3' = c(11.6, 11.6, 11.6),
                                 'min_tcrit' = c(9.9, 9.9, 10.0))

  expect_equal(estimate_crit_time(DRT_example, tol_factor = 10, TRT = DRT_example), reference_result)
})


test_that("wrong input data format results in errors", {
  example_data <- data.frame('psi.D' = TRUE, 'psi.R' = 1.5, 'psi.T' = 1.5,
                             'D.0' = 1e5, 'R.0' = 1e5,
                             'gamma.D' = 1e-11, 'gamma.T' = 1e-9)

  expect_error(estimate_crit_time(DRT = example_data, tol_factor = 10), "non-numeric columns")
  expect_error(estimate_crit_time(DRT = NULL, TRT = NULL, tol_factor = 10), "must be a dataframe.")
})

test_that("missing input data columns results in errors", {
  expect_error(estimate_crit_time(DRT_example[,c("D.0")]), "must be a dataframe")
  expect_error(suppressWarnings(estimate_crit_time(DRT_example[,c("D.0", "psi.D")])), "compute at least one critical time")
  expect_error(suppressWarnings(estimate_crit_time(DRT_example[1:5], TRT = NULL, tol_factor = 10)), "at least one critical time.")
})

test_that("partial input data columns results in warnings", {
  partial_data <- data.frame('psi.D' = 1, 'psi.R' = 1.5, 'psi.T' = 1.5,
                             'D.0' = 1e5, 'R.0' = 1e5)

  expect_warning(estimate_crit_time(cbind(partial_data, 'gamma.D' = 1e-11)), "needed to compute tcrit1")
  expect_warning(estimate_crit_time(cbind(partial_data, 'gamma.T' = 1e-9)), "needed to compute tcrit2")
})

###########################################################
test_that("estimate output is a dataframe with the right columns", {
  result <- estimate_crit_time(DRT_example, TRT = TRT_example, tol_factor = 10)
  expect_equal(colnames(result), c("ID", "gamma.D", "gamma.T", "tcrit1", "tcrit2", "tcrit3", "min_tcrit"))
  expect_equal(dim(result), c(nrow(DRT_example), 7))
})

test_that("ID cols are present in output", {
  result <- estimate_crit_time(DRT_example, TRT = TRT_example, tol_factor = 10)
  expect_equal(result$ID, DRT_example$ID)

  test_DRT <- DRT_example
  test_DRT[, 'ID_2'] <- 'test'
  expect_error(estimate_crit_time(test_DRT, TRT = TRT_example, tol_factor = 10, id_cols = c('ID', 'ID_2')),
               "id_cols of DRT and TRT do not match.")
  test_TRT <- TRT_example
  test_TRT[, 'ID_2'] <- 'test'
  result <- estimate_crit_time(test_DRT, TRT = test_TRT, tol_factor = 10, id_cols = c('ID', 'ID_2'))
  expect_true(all(c('ID', 'ID_2') %in% colnames(result)))

})


test_that("IDs are matched properly", {
  newTRT <- rbind(TRT_example, cbind('ID' = factor(1:3), TRT_example[2:10]))
  expect_warning(estimate_crit_time(DRT_example, TRT = newTRT, tol_factor = 10), 'id_cols of DRT and TRT are not equal')
  result_newTRT <- suppressWarnings(estimate_crit_time(DRT_example, TRT = newTRT, tol_factor = 10))

  newDRT <- rbind(DRT_example, cbind('ID' = factor(1:3), DRT_example[2:10]))
  expect_warning(estimate_crit_time(newDRT, TRT = TRT_example, tol_factor = 10), 'id_cols of DRT and TRT are not equal')
  result_newDRT <- suppressWarnings(estimate_crit_time(newDRT, TRT = TRT_example, tol_factor = 10))

  expect_equal(dim(result_newTRT), c(nrow(DRT_example), 7))
  expect_equal(dim(result_newDRT), c(nrow(DRT_example), 7))
  expect_error(estimate_crit_time(DRT_example, TRT = cbind('ID' = factor(1:3), TRT_example[2:10]),
                                  tol_factor = 10), 'No matching IDs')
})

