test_that("findtsr_peaks finds the correct starts!", {

  x <- S4Vectors::Rle(c(rep(0,10),2,0,0,0,0,3,3,6,3,3,0,0,0,0,2,rep(0,10)))
  y <- findtsr_peaks(x, bg_size = 15, height_above_bg = 2)

  expect_equal(start(y), c(10,15,24))
})

test_that("findtsr_peaks finds the correct ends!", {

  x <- S4Vectors::Rle(c(rep(0,10),2,0,0,0,0,3,3,6,3,3,0,0,0,0,2,rep(0,10)))
  y <- findtsr_peaks(x, bg_size = 15, height_above_bg = 2)

  expect_equal(end(y), c(12,21,26))
})
