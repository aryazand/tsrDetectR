test_that("strand_coverage produces 2 object", {
  expected <- 2
  result <- strand_coverage(cmv_proseq_sample) |> length()
  expect_equal(expected, result)
})
