test_that("runfun produces correct Rle object", {
  expect_equal(as.vector(runfun(Rle(1:10), width =2, FUN = sum)),
               c(3,5,7,9,11,13,15,17,19))
})
