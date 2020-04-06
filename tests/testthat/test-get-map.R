context("gwa_get_map")

test_that("gwa_get_map stops and prints error for invalid pop", {
  expect_error(gwa_get_map("NON"),
               "Invalid population code. See https://www.internationalgenome.org/faq/which-populations-are-part-your-study/")
})
