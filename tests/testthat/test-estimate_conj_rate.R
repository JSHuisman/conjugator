###########################################################
## Unit tests for the conjugation rate estimation
## Authors: Jana S. Huisman
###########################################################
context("Conjugation Rate Estimation")
library(conjugator)

###########################################################

test_that("wrong input data format results in errors", {
  expect_error(estimate_conj_rate(DRT_example[,c("D.t")], "TD"), "does not contain enough data columns")
  expect_error(estimate_conj_rate(c(0.4, 0.7, 0.8), "ASM"), "does not contain enough data columns")
  expect_error(estimate_conj_rate(data.frame(D.t = c(100, 1000), T.t = c(NA, NA)), "TD"), "non-numeric columns")
})

test_that("missing input data columns results in errors", {
  expect_silent(estimate_conj_rate(DRT_example[,c("D.t", "T.t")], "TD"))
  expect_error(estimate_conj_rate(DRT_example[,c("D.t", "T.t")], "ASM"), "Missing strictly necessary columns")
  expect_error(estimate_conj_rate(DRT_example[,c("D.t", "T.t")], "SM"), "Missing strictly necessary columns")
  expect_error(estimate_conj_rate(DRT_example[,c("D.t", "T.t")], "T_DR"), "Missing strictly necessary columns")
  expect_error(estimate_conj_rate(DRT_example[,c("D.t", "T.t")], "T_RT"), "Missing strictly necessary columns")
  expect_error(estimate_conj_rate(DRT_example[,c("D.t", "T.t")], "Dionisio"), "Missing strictly necessary columns")
  expect_error(estimate_conj_rate(DRT_example[,c("D.t", "T.t")], "Gama"), "Missing strictly necessary columns")
})

test_that("missing method argument results in default", {
  expect_equal(estimate_conj_rate(DRT_example), estimate_conj_rate(DRT_example, c("SM", "ASM")))
})

test_that("unknown method results in error", {
  expect_error(estimate_conj_rate(DRT_example, "test"))
})

###########################################################

test_that("Zero transconjugants results in zero rates", {
  bad_example = DRT_example
  bad_example$T.t = 0

  expect_equal(estimate_conj_rate(bad_example, "SM")$estimate, c(0, 0, 0))
  expect_equal(estimate_conj_rate(bad_example, "SM")$estimate, estimate_conj_rate(bad_example, "ASM")$estimate)
})


test_that("ASM is the same as SM for small, equal growth rates", {
  create_equal_example <- function(x) {
    all_equal_example = DRT_example
    all_equal_example[, c('psi.R', 'psi.D', 'psi.T')] = x
    return(all_equal_example)}

  expect_equal(estimate_conj_rate(create_equal_example(1), "SM")$estimate,
               estimate_conj_rate(create_equal_example(1), "ASM")$estimate, tolerance=1e-10)
  expect_equal(estimate_conj_rate(create_equal_example(0.01), "SM")$estimate,
               estimate_conj_rate(create_equal_example(0.01), "ASM")$estimate, tolerance=1e-13)

  # why does this fail?? because the example is not rescaled for large final populations
  # expect_equal(estimate_conj_rate(create_equal_example(10), "SM")$estimate,
  #              estimate_conj_rate(create_equal_example(10), "ASM")$estimate, tolerance=1e-8)
})

test_that("SM uses maximum growth rate", {
  expect_silent(estimate_conj_rate(DRT_example[,setdiff(colnames(DRT_example), "psi.R")], "SM"))
  expect_equal(estimate_conj_rate(DRT_example[,setdiff(colnames(DRT_example), "psi.D")], "SM"), estimate_conj_rate(DRT_example, "SM"))
})

###########################################################

test_that("any combination of methods runs", {
  methods = c("SM", "ASM", "TD", "T_DR", "T_RT", "Dionisio", "Gama")
  expect_silent(estimate_conj_rate(DRT_example, sample(methods, 3)) )
  expect_silent(estimate_conj_rate(DRT_example, sample(methods, 7, replace = T)) )
})

test_that("estimate output is a dataframe with the right columns", {
  result <- estimate_conj_rate(DRT_example)

  test_example <- DRT_example
  colnames(test_example)[1] <- 'test'
  test_result <- estimate_conj_rate(test_example, c('SM', 'ASM'))

  expect_equal(colnames(result), c('ID', 'estimate', 'method'))
  expect_equal(colnames(test_result), c('ID', 'estimate', 'method'))
  expect_equal(class(result$ID), 'factor')
  expect_equal(class(result$estimate), 'numeric')
  expect_equal(class(result$method), 'character')
})

test_that("all ID cols are present in output", {
  test_DRT <- DRT_example
  test_DRT[, 'ID_2'] <- 'test'
  result <- estimate_conj_rate(test_DRT, id_cols = c('ID', 'ID_2'))
  expect_true(all(c('ID', 'ID_2') %in% colnames(result)))
  expect_false(all(c('ID', 'ID_2') %in% colnames(estimate_conj_rate(test_DRT) ) ))
})
