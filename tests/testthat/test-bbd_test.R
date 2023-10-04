test_that("Criged-block decomposition approach works", {
  load("bbd.RData")
  bbd_res <- solver_bbd(S, Lambda)
  expect_equal(bbd_res_check[c("X_est")],
               bbd_res[c("X_est")],
               tolerance = 1e-5)
})
