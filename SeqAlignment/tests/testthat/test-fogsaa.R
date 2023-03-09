test_that("apply the FOGSAA algorithm", {
  expected <- list('ACGGTTGC', 'A-GCGT-C')
  expect_equal(fogsaa_alignment('ACGGTTGC','AGCGTC', 1, -1, -2), expected)
})

test_that("returns error if not character", {
  input <- c('AGT', 1, 1, -1, -2)
  expect_error(fogsaa_alignment(input))
})

test_that("returns error if not admitted characters", {
  input <- c('AGT', 'ACR', 1, -1, -2)
  expect_error(fogsaa_alignment(input))
})

test_that("gap value is negative", {
  expect_error(fogsaa_alignment('AGT', 'ACT', 1, -1, 2))
})

test_that("match score is greater than mismatch", {
  expect_error(fogsaa_alignment('AGT', 'ACT', 1, 2, -1))
})






