test_that("returns error if wrong score", {
  expected <- -2
  expect_equal(alignment_score('ACGGTTGC', 'A-GCGT-C', 1, -1, -2), expected)
})

test_that("gap value is negative", {
  expect_warning(alignment_score('AGT', 'ACT', 1, -1, 2))
})

test_that("match score is greater than mismatch", {
  expect_warning(alignment_score('AGT', 'ACT', 1, 2, -1))
})

