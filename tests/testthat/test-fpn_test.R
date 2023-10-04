test_that("Fast Projected Newton-like Algorithm works", {
  load("fpn.RData")
  fpn_res <- solver_fpn(S, Lambda)
  expect_equal(fpn_res_check[c("X_est", "objective", "converge")],
               fpn_res[c("X_est", "objective", "converge")],
               tolerance = 1e-5)
})
