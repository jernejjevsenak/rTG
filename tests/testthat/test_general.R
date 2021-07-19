library(rTG)
library(testthat)

data(parameters)
data(data_trees)

simulation_1 <- XPSgrowth(data_trees = data_trees,
      parameters = parameters,
      ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
      fitting_method = c("gompertz", "gam", "brnn"),
      fitted_save = FALSE,
      search_initial_gom = TRUE,
      add_zeros = TRUE,
      post_process = TRUE)

expect_type(simulation_1, "list")
expect_type(simulation_1$fitted, "list")
expect_type(simulation_1$gompertz_grid_search, "list")
expect_type(simulation_1$gompertz_grid_search_errors, "character")
