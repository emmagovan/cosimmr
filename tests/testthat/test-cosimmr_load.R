co <- function(expr) capture.output(expr, file = NULL)

test_that("simplest example", {
  data("geese_data_day1")
  
  # Load in with simmr_load
  cosimmr_1 <- with(
    geese_data_day1,
    cosimmr_load(
      formula = mixtures ~ 1,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds,
      correction_means = correction_means,
      correction_sds = correction_sds,
      concentration_means = concentration_means
    )
  )
  expect_s3_class(cosimmr_1, "cosimmr_input")
  expect_true(is.matrix(cosimmr_1$source_means))
  expect_true(is.matrix(cosimmr_1$source_sds))
  expect_true(is.matrix(cosimmr_1$correction_means))
  expect_true(is.matrix(cosimmr_1$correction_sds))
  expect_true(is.matrix(cosimmr_1$concentration_means))
})

test_that("test without tefs and concentration", {
  data("geese_data_day1")
  
  # Load in with simmr_load
  cosimmr_2 <- with(
    geese_data_day1,
    cosimmr_load(
      formula = mixtures ~ 1,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds
    )
  )
  expect_s3_class(cosimmr_2, "cosimmr_input")
  expect_true(is.matrix(cosimmr_2$source_means))
  expect_true(is.matrix(cosimmr_2$source_sds))
  expect_true(is.matrix(cosimmr_2$correction_means))
  expect_true(is.matrix(cosimmr_2$correction_sds))
  expect_true(is.matrix(cosimmr_2$concentration_means))
})

# Test if you're accidentally giving it text data

test_that("test without tefs and concentration", {
  data("geese_data_day1")
  geese_data_day1_tmp <- geese_data_day1
  geese_data_day1_tmp$mixtures <- matrix(as.character(geese_data_day1$mixtures), nrow = 9, ncol = 2)
  
  # Load in with simmr_load
  expect_error(with(
    geese_data_day1_tmp,
    cosimmr_load(
      formula = mixtures ~ 1,
      source_names = source_names,
      source_means = source_means,
      source_sds = source_sds
    )
  ))
})

# test_that("test it works from excel", {
#   library(readxl)
#   path <- system.file("extdata", "geese_data.xls", package = "simmr")
#   geese_data <- lapply(excel_sheets(path), read_excel, path = path)
#   
#   targets <- as.matrix(geese_data[[1]])
#   sources <- as.matrix(geese_data[[2]])
#   TEFs <- as.matrix(geese_data[[3]])
#   concdep <- as.matrix(geese_data[[4]])
#   
#   geese_cosimmr <- cosimmr_load(
#     formula = (targets[, 1:2]) ~as.factor(paste("Day", targets[,3])),
#     source_names = sources[,1],
#     source_means = sources[, 2:3],
#     source_sds = sources[, 4:5],
#     correction_means = TEFs[, 2:3],
#     correction_sds = TEFs[, 4:5],
#     concentration_means = concdep[, 2:3]
#   )
#   co(cosimmr_out <- cosimmr_ffvb(geese_cosimmr))
#   expect_s3_class(cosimmr_out, "cosimmr_output")
# })