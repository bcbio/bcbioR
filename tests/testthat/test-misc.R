library(bcbioR)


test_that("rnaseq samplesheet", {

  expect_no_error(bcbio_nfcore_check(system.file("extdata", "rnaseq_good.csv", package = "bcbioR")))
  expect_error(bcbio_nfcore_check(system.file("extdata", "rnaseq_missing_col.csv", package = "bcbioR")))
  expect_error(bcbio_nfcore_check(system.file("extdata", "rnaseq_wnumber.csv", package = "bcbioR")))
  expect_warning(bcbio_nfcore_check(system.file("extdata", "rnaseq_na.csv", package = "bcbioR")))
})
