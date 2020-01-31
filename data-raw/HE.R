rm(list = ls())
cat("\014")

x_train <- readr::read_table2(file = "data-raw/HeatExchange/X40_orig_scale",
                              col_names  = FALSE) %>%
  as.matrix()

y_train <- readr::read_table2(file = "data-raw/HeatExchange/Ysim40",
                              col_names  = FALSE) %>%
  as.matrix()

x_test <- readr::read_table2(file = "data-raw/HeatExchange/xpred24",
                             col_names  = FALSE) %>%
  as.matrix()

y_test <- readr::read_table2(file = "data-raw/HeatExchange/ypred24",
                             col_names  = FALSE) %>%
  as.matrix()


HeatExchange <- list(x = x_train, y = as.vector(y_train),
                     x_test = x_test, y_test = as.vector(y_test))

usethis::use_data(HeatExchange, overwrite = TRUE)
