rm(list = ls())
cat("\014")

x_train <- readr::read_table2(file = "data-raw/WingWeight/X50",
                              col_names  = FALSE)

y_train <- readr::read_table2(file = "data-raw/WingWeight/Ysim50",
                              col_names  = FALSE)

x_test <- readr::read_table2(file = "data-raw/WingWeight/xpred150",
                             col_names  = FALSE)

y_test <- readr::read_table2(file = "data-raw/WingWeight/ypreds150",
                             col_names  = FALSE)


WingWeight <- list(x = x_train, y = as.vector(y_train),
                   x_test = x_test, y_test = as.vector(y_test))

usethis::use_data(WingWeight, overwrite = TRUE)
