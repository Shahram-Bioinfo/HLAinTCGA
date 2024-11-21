test_that("allele_freq works correctly with package data", {
  # Load the dataset from the package's data folder
  data("data", package = "HLAinTCGA")  # Ensure the dataset is loaded
  
  # Check if the dataset is properly loaded and is a data frame
  expect_true(exists("data"), "The data object should exist in the package.")
  expect_true(is.data.frame(data), "The data object should be a data frame.")
  
  # Test the allele_freq function with an example allele
  result <- allele_freq("A*01:01", data)
  
  # Verify the output is a list
  expect_type(result, "list")
  
  # Check if the output list contains the expected components
  expect_true("results_table" %in% names(result), "'results_table' should be in the result list.")
  expect_true("plot" %in% names(result), "'plot' should be in the result list.")
  
  # Ensure the results_table is a data frame
  expect_true(is.data.frame(result$results_table), "'results_table' should be a data frame.")
  
  # Ensure the plot is a ggplot object
  expect_s3_class(result$plot, "ggplot")  # Removed the `info` argument
})

