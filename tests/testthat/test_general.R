library(rTG)
library(testthat)

data(parameters)
data(data_trees)

simulation_1 <- XPSgrowth(data_trees = data_trees,
      parameters = parameters,
      ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
      fitting_method = c("gompertz", "gam", "brnn"),
      fitted_save = FALSE,
      search_initial_gom = FALSE,
      add_zeros = TRUE,
      post_process = TRUE)


# 1 Check the correct outputs
expect_type(simulation_1, "list")
expect_type(simulation_1$fitted, "list")
expect_type(simulation_1$gompertz_initial_parameters, "list")
expect_type(simulation_1$gompertz_initial_parameters_errors, "character")

# 2  If unified_parameters is used, you must provide initial parameters
expect_error(
XPSgrowth(data_trees = data_trees,
          unified_parameters = TRUE,
          ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
          fitting_method = c("gompertz", "gam", "brnn"),
          fitted_save = FALSE,
          search_initial_gom = FALSE,
          add_zeros = TRUE,
          post_process = TRUE))

