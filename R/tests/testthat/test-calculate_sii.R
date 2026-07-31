testthat::test_that("calculate_sii returns one result per grouping level", {
  
  input <- tibble::tibble(
    indicator_id = rep(1L, 5),
    start_date = rep(as.Date("2025-01-01"), 5),
    end_date = rep(as.Date("2025-12-31"), 5),
    aggregation_id = rep(1L, 5),
    age_group_code = rep(1L, 5),
    sex_code = rep(1L, 5),
    ethnicity_code = rep(999L, 5),
    source_code = rep(1L, 5),
    imd_code = 1:5,
    numerator = c(40, 45, 50, 55, 60),
    denominator = rep(100, 5)
  )
  
  result <- calculate_sii(
    data = input,
    quintile_col = imd_code,
    numerator_col = numerator,
    denominator_col = denominator
  )
  
  testthat::expect_equal(nrow(result), 1)
  testthat::expect_true(
    is.numeric(result$sii_signed_percentage_points)
  )
})