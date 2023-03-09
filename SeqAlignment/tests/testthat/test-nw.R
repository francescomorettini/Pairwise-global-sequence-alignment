test_that("returns error if not character", {
  input <- c('AGT', 1, 1, -1, -2)
  expect_error(nw_alignment(input))
})

test_that("apply the NW algorithm", {
  expected <- list('ACGGTTGC','A-GCGT-C')
  expect_identical(nw_alignment('ACGGTTGC','AGCGTC', 1, -1, -2), expected)
})

test_that("returns error if not admitted characters", {
  input <- c('AGT', 'ACR', 1, -1, -2)
  expect_error(nw_alignment(input))
})

test_that("gap value is negative", {
  expect_warning(nw_alignment('AGT', 'ACT', 1, -1, 2))
})

test_that("match score is greater than mismatch", {
  expect_warning(nw_alignment('AGT', 'ACT', 1, 2, -1))
})