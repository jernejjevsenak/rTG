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

expect_is(simulation_1, "xpsg")
expect_is(simulation_1[[1]], "data.frame")
expect_is(simulation_1[[2]], "logical")
expect_is(simulation_1[[3]], "character")
