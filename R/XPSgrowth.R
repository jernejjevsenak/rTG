#' XPSgrowth
#'
#' XylemPhloemSeasonalGrowth: This Function fits and compares the selected
#' methods for modeling seasonal xylem and phloem data.
#' @param data_trees a data frame with ID variables and wood formation data with
#' columns doy and width
#' @param parameters a data frame with ID variables and initial parameter values
#' for the selected methods
#' @param search_initial_gom logical, should the algorithm to search initial
#' Gompertz parameters be applied? This argument also overwrites manually
#' defined Gompertz parameter values
#' @param search_initial_double_gom logical, should the algorithm to search
#' initial parameters for double Gompertz function be applied? This argument
#' also overwrites manually defined parameter values for double Gompertz
#' @param fitting_method vector of one or more methods to be compared:
#' "gompertz", "double_gompertz", "gam", "brnn"
#' @param ID_vars character vector of variables which indicate column names of
#' ID variables
#' @param fitted_save logical, should the fitted curves be saved in current
#' working directory?
#' @param add_zeros logical, should zero observations at the beginning of
#' growing season be added?
#' @param add_zeros_before if 'min' (character) then zeros will be added prior
#' to the first observation in each year. Alternatively, users can specify
#' absolute doy prior which zeros will be added.
#' @param post_process logical, should the post-process algorithm be applied?
#' @param unified_parameters logical, if TRUE, the algorithm will use only
#' manually selected function parameters. See the arguments 'gom_a', 'gom_b',
#' 'd_gom_k', 'd_gom_a1', 'd_gom_a2', 'd_gom_b1', 'd_gom_b2', 'd_gom_k1',
#' 'd_gom_k2', 'brnn_neurons', 'gam_k' and 'gam_sp'. Default is FALSE
#' @param gom_a numeric, the parameter a for the Gompertz function
#' @param gom_b numeric, the parameter b for the Gompertz function
#' @param gom_k numeric, the parameter k for the Gompertz function
#' @param d_gom_a1 numeric, the parameter a1 for the double Gompertz function
#' @param d_gom_a2 numeric, the parameter a2 for the double Gompertz function
#' @param d_gom_b1 numeric, the parameter b1 for the double Gompertz function
#' @param d_gom_b2 numeric, the parameter b2 for the double Gompertz function
#' @param d_gom_k1 numeric, the parameter k1 for the double Gompertz function
#' @param d_gom_k2 numeric, the parameter k2 for the double Gompertz function
#' @param brnn_neurons positive integer, the number of neurons to be used by
#' the BRNN method
#' @param gam_k numeric, the parameter k for General Additive Model (GAM)
#' @param gam_sp numeric, the parameter sp for General Additive Model (GAM)
#' @param gom_a_range a numerical vector of the possible values of the
#' parameter a, which is considered in the search for the initial Gompertz
#' parameter values
#' @param gom_b_range a numerical vector of the possible values of the
#' parameter b, which is considered in the search for the initial Gompertz
#' parameter values
#' @param gom_k_range a numerical vector of the possible values of the
#' parameter k, which is considered in the search for the initial Gompertz
#' parameter values
#' @param d_gom_a1_range A numerical vector representing the range of potential
#' values for the 'a1' parameter within the double Gompertz function.
#' @param d_gom_a2_range A numerical vector representing the range of potential
#' values for the 'a2' parameter within the double Gompertz function.
#' @param d_gom_b1_range A numerical vector representing the range of potential
#' values for the 'b1' parameter within the double Gompertz function.
#' @param d_gom_b2_range A numerical vector representing the range of potential
#' values for the 'b2' parameter within the double Gompertz function.
#' @param d_gom_k1_range A numerical vector representing the range of potential
#' values for the 'k1' parameter within the double Gompertz function.
#' @param d_gom_k2_range A numerical vector representing the range of potential
#' values for the 'k2' parameter within the double Gompertz function.
#'
#' @return a list with the following elements:
#' \enumerate{
#'  \item $fitted - a data frame with fitted values
#'  \item $gompertz_initial_parameters - a data frame that contains a curated selection of initial parameter values for the Gompertz function.
#'  \item $gompertz_model_parameters - a data frame with final model coefficients for the Gompertz function.
#'  \item $gompertz_initial_parameters_errors - a data frame with unsuccessful cases of Gompertz grid search.
#'  \item $double_gompertz_initial_parameters - a data frame that contains a curated selection of initial parameter values for the double Gompertz function.
#'  \item $double_gompertz_initial_parameters_errors - a data frame with unsuccessful cases of double Gompertz grid search.
#'  \item $doble_gompertz_model_parameters - a data frame with final model coefficients for the double Gompertz function.
#'}
#'
#' @export
#'
#' @examples
#' library(rTG)
#'
#' # 1 Example on xylem and phloem data
#' data(parameters)
#' data(data_trees)
#'
#' # subset data_trees
#' data_trees <- data_trees[c(1:27),]
#'
#' # 1a Example using neural network
#' simulation_1a <- XPSgrowth(data_trees = data_trees,
#'      parameters = parameters,
#'      ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
#'      fitting_method = c("brnn"),
#'      fitted_save = FALSE,
#'      search_initial_gom = FALSE,
#'      add_zeros = TRUE,
#'      add_zeros_before = 'min',
#'      post_process = TRUE)
#'
#' \dontrun{
#' #' # 1b Example on double Gompertz function
#' simulation_1b <- XPSgrowth(data_trees = data_trees,
#'      parameters = parameters,
#'      ID_vars = c("Species", "Tissue", "Site", "Year"),
#'      fitting_method = c("double_gompertz"),
#'      fitted_save = FALSE,
#'      search_initial_double_gom = FALSE,
#'      unified_parameters = TRUE,
#'      add_zeros = TRUE,
#'      add_zeros_before = 'min',
#'      d_gom_a1 = 0.204, d_gom_a2 = 0.240,
#'      d_gom_b1 = 2.433, d_gom_b2 = 2.900,
#'      d_gom_k1 = 0.974, d_gom_k2 = 0.963,
#'      post_process = TRUE)
#'
#' # 1b Example on Double Gompertz function without initial parameters
#' simulation_1c <- XPSgrowth(data_trees = data_trees,
#'      parameters = parameters,
#'      ID_vars = c("Species", "Tissue", "Site", "Year"),
#'      fitting_method = c("double_gompertz"),
#'      fitted_save = FALSE,
#'      search_initial_double_gom = TRUE,
#'      post_process = TRUE)
#'
#' # Obtain model parameters
#' simulation_1c$double_gompertz_model_parameters
#' }
#'
#' # 2 Example on dendrometer data
#' data("data_dendrometers")
#'
#' simulation_2 <- XPSgrowth(data_dendrometers, unified_parameters = TRUE,
#'                   ID_vars = c("site", "species", "year", "tree"),
#'                   fitting_method = c("brnn", "gam"),
#'                   brnn_neurons = 2, gam_k = 9, gam_sp = 0.5,
#'                   search_initial_gom = TRUE, add_zeros = FALSE,
#'                   post_process = TRUE)


XPSgrowth <- function(data_trees, parameters = NULL,
                 search_initial_gom = FALSE,
                 search_initial_double_gom = FALSE,
                 fitting_method = c("gompertz", "GAM", "brnn", "double_gompertz"),
                 ID_vars = NULL,
                 fitted_save = FALSE,
                 add_zeros = TRUE,
                 add_zeros_before = 'min',
                 post_process = TRUE,
                 unified_parameters = FALSE,
                 gom_a = NA, gom_b = NA, gom_k = NA,
                 d_gom_a1 = NA, d_gom_a2 = NA,
                 d_gom_b1 = NA, d_gom_b2 = NA,
                 d_gom_k1 = NA, d_gom_k2 = NA,
                 brnn_neurons = NA,
                 gam_k = NA, gam_sp = NA,
                 gom_a_range = seq(0, 3000, by = 500),
                 gom_b_range = seq(0.01, 1000, by = 50),
                 gom_k_range = seq(0, 500, by = 2),
                 d_gom_a1_range = seq(0, 1, by = 0.001),
                 d_gom_a2_range = seq(0, 5, by = 0.01),
                 d_gom_b1_range = seq(0, 5, by = 0.001),
                 d_gom_b2_range = seq(0, 10, by = 0.1),
                 d_gom_k1_range = seq(0, 1, by = 0.001),
                 d_gom_k2_range = seq(0, 1, by = 0.001)
                 ){

  # Defining global variables
  width <- NULL
  doy <- NULL
  key <- NULL
  first_diff <- NULL
  width_pred <- NULL
  neg_ind <- NULL
  errors_grid <- NA
  txtProgressBar <- NULL
  setTxtProgressBar <- NULL
  lines <- NULL
  note <- NULL
  nls.lm.control <- NULL
  nlsLM <- NULL

  # Progress bar
  pb <- txtProgressBar(min = 0, max = length(fitting_method), style = 3)

  # make sure you have doy and width in lowercase
  if ("DOY" %in% colnames(data_trees)){data_trees <- rename(data_trees, "doy" = "DOY")}
  if ("Doy" %in% colnames(data_trees)){data_trees <- rename(data_trees, "doy" = "Doy")}
  if ("WIDTH" %in% colnames(data_trees)){data_trees <- rename(data_trees, "width" = "WIDTH")}
  if ("Width" %in% colnames(data_trees)){data_trees <- rename(data_trees, "width" = "Width")}

  # If ID_vars is null, we assume that all variables other than doy and width in data_trees are ID_vars
  if (is.null(ID_vars)){
    ID_vars <- colnames(data_trees)
    ID_vars <- ID_vars[!(ID_vars %in% c("doy", "width"))]
  }


  if (sum(is.na(data_trees)) > 0){
    stop("Please remove missing values from data_trees")
  }

  # Columns doy & width must be present in data_trees
  if (!("doy" %in% colnames(data_trees))){
    stop("Column 'doy' is missing in data_trees")
  }

  if (!("width" %in% colnames(data_trees))){
    stop("Column 'width' is missing in data_trees")
  }


  # Just in case, convert fitting methods to lowercase
  fitting_method <- tolower(fitting_method)

  # If parameters table is provided, we make sure that all combinations are present
  if (!is.null(parameters)){
  distinct_combinations_data <- data_trees %>% select(all_of(ID_vars)) %>% distinct()
  distinct_combinations_parameters <- parameters %>% select(all_of(ID_vars)) %>% distinct()

  missing_rows <- anti_join(distinct_combinations_data, distinct_combinations_parameters,
                            by = colnames(distinct_combinations_data))

  if (nrow(missing_rows) > 0){

    stop(
      "data_trees and parameters do not match\n
      The following rows are missing:\n",
      paste(apply(missing_rows, 1, paste, collapse = " | "), collapse = "\n")
    )

  }

  }

  # In case parameters table is missing, we simulate one
  if (is.null(parameters)){

    parameters <- data_trees %>% select(all_of(ID_vars)) %>% distinct() %>%

      mutate(
        gom_a =  3000,
        gom_b = 2000,
        gom_k = 10,

        d_gom_a1 = 3000,
        d_gom_a2 = 3000,
        d_gom_b1 = 2000,
        d_gom_b2 = 2000,
        d_gom_k1 = 10,
        d_gom_k2 = 10,

        brnn_neurons = 3,
        gam_k = 10,
        gam_sp = 0.5)

    if ("gam" %in% fitting_method & unified_parameters == FALSE){

      warning("No parameters for GAM were provided. The following default values will be used: gam_k = 10 and gam_sp = 0.5")

    }

    if ("brnn" %in% fitting_method & unified_parameters == FALSE){

      warning("No parameters for brnn were provided. The following default values will be used: brnn_n = 3")

    }

    if ("gompertz" %in% fitting_method & unified_parameters == FALSE & search_initial_gom == FALSE){

      warning("the parameter search_initial_gompertz was set to TRUE")

      search_initial_gom = TRUE

    }


    if ("double_gompertz" %in% fitting_method & unified_parameters == FALSE & search_initial_double_gom == FALSE){

      warning("the parameter search_initial_double_gom was set to TRUE")

      search_initial_double_gom = TRUE

    }




  }


  # 1 are all ID_vars in both tables?
  if(unified_parameters == FALSE){

    if (sum(ID_vars %in% colnames(parameters) == FALSE)>0){
      stop("check your ID_vars, they don't exist in your parameters")
    }
  }

  if (sum(ID_vars %in% colnames(data_trees) == FALSE)>0){
    stop("check your ID_vars, they don't exist in your data_treses")
  }


  # If you use unified parameters, you must provide them
  if (unified_parameters == TRUE){

    if ("gompertz" %in% fitting_method){

      tm_g_p <- c(gom_a, gom_b, gom_k)

      if (sum(is.na(tm_g_p) > 0)   & search_initial_gom == FALSE){

        stop("If unified_parameters is used, you must provide Gompertz parameters")

      }
    }


    if ("double_gompertz" %in% fitting_method){

      tm_d_g_p <- c(d_gom_a1, d_gom_a2, d_gom_b1, d_gom_b2, d_gom_k1, d_gom_k2)

      if (sum(is.na(tm_d_g_p) > 0) & search_initial_double_gom == FALSE){

        stop("If unified_parameters is used, you must provide double Gompertz parameters")

      }
    }


    if ("brnn" %in% fitting_method){

      if (sum(is.na(brnn_neurons) > 0)){

        stop("If unified_parameters is used, you must provide brnn_neurons parameter")

      }
    }


    if ("gam" %in% fitting_method){

      tm_g_p <- c(gam_k, gam_sp)

      if (sum(is.na(tm_g_p) > 0)){

        stop("If unified_parameters is used, you must provide GAM parameters")


      }
    }

  }

  # Create data key

  if (length(ID_vars) > 1){

    data_trees$key <- apply(data_trees[, ID_vars], 1,
                            paste0, sep = "", collapse = "_")

    data_trees$key <- gsub(" ", "", data_trees$key) # remove empty characters

  } else {

    data_trees$key <- data_trees[, ID_vars]

  }

  if (unified_parameters == TRUE){

    parameters <- data.frame(key = unique(data_trees$key))

    parameters$gom_a <- gom_a
    parameters$gom_b <- gom_b
    parameters$gom_k <- gom_k
    parameters$d_gom_a1 <- d_gom_a1
    parameters$d_gom_a2 <- d_gom_a2
    parameters$d_gom_b1 <- d_gom_b1
    parameters$d_gom_b2 <- d_gom_b2
    parameters$d_gom_k1 <- d_gom_k1
    parameters$d_gom_k2 <- d_gom_k2
    parameters$brnn_neurons <- brnn_neurons
    parameters$gam_k <- gam_k
    parameters$gam_sp <- gam_sp

  } else {

    if (length(ID_vars) > 1){

      parameters$key <- apply(parameters[, ID_vars], 1,
                              paste0, sep = "", collapse = "_")

      parameters$key <- gsub(" ", "", parameters$key) # remove empty characters

    } else {

      parameters$key <- parameters[, ID_vars]

    }

  }

  # Gompertz solution list
  list_temps <- list()
  list_errors <- list()
  list_solutions <- list()
  list_solution_a <- list()
  list_solution_b <- list()
  list_solution_k <- list()

  # Gompertz list for extracted parameters
  list_extracted <- list()
  list_extracted_a <- list()
  list_extracted_b <- list()
  list_extracted_k <- list()

  # double Gompertz
  list_errors_dg <- list()
  list_solutions_dg <- list()
  list_solution_a1 <- list()
  list_solution_a2 <- list()
  list_solution_b1 <- list()
  list_solution_b2 <- list()
  list_solution_k1 <- list()
  list_solution_k2 <- list()
  p3 = 1; p4 = 1; p3a = 1

  # empty list for double Gompertz extracted parameters
  list_extracted_dg <- list()
  list_extracted_a1 <- list()
  list_extracted_a2 <- list()
  list_extracted_b1 <- list()
  list_extracted_b2 <- list()
  list_extracted_k1 <- list()
  list_extracted_k2 <- list()

  # I define those two objects as NA
  # If grid search is used, they are overwritten
  errors_grid <- NA
  final_parameters <- NA
  extracted_parameters <- NA

  errors_grid_dg <- NA
  final_parameters_dg <- NA
  extracted_parameters_dg<- NA

  b_holder = 1
  p = 1
  p2 = 1
  p2a = 1
  pbar_holder = 1

  unique_keys <- unique(data_trees$key)

for (ut in fitting_method){

  setTxtProgressBar(pb, pbar_holder)

     current_fitting_method = ut

if (current_fitting_method == "gompertz"){

  for (i in unique_keys){

    par_grid <- expand.grid(gom_a = gom_a_range,
                            gom_b = gom_b_range,
                            gom_k = gom_k_range)

    temp_parameters <- parameters[parameters$key == i,]

    gom_a <- temp_parameters$gom_a
    gom_b <- temp_parameters$gom_b
    gom_k <- temp_parameters$gom_k

    temp_data <- data_trees[data_trees$key == i,]

    # specify the measurement type - so we can distinguish from added zeros
    temp_data$note <- "raw measurement"

    # add zeros at the beginning
    if(add_zeros == TRUE){

      if (add_zeros_before == 'min'){

        min_doy <- min(temp_data$doy) - 1

      } else {

        if (!is.numeric(add_zeros_before) & add_zeros_before < 0){

          stop("The argument 'add_zeros_before' should be numeric and greater than 0")

        }

        min_doy <- add_zeros_before

      }

      row_list <- list()
      for (J in 1:min_doy){
        temp_row <- temp_data[1,]
        row_list[[J]] <- temp_row
      }

      new_rows <- do.call(rbind, row_list)
      new_rows$doy <- c(1:min_doy)
      new_rows$width <- 0
      new_rows$note <- "added zero"
      temp_data <- rbind(new_rows, temp_data)
    }

    capture.output(output <- try(nls(width ~ a*exp(-exp(b - k*doy)), data = temp_data,
                      start=list(a = gom_a, b = gom_b, k = gom_k),
                      nls.control(maxiter = 1000, tol = 1e-05,
                                  minFactor = 1/1024, printEval = FALSE,
                                  warnOnly = FALSE)), silent = TRUE))

    # When all parameters are null, the class can still be nls. Additional check is needed
    parm_test <- c(gom_a, gom_b, gom_k)

    if (is(output, "nls") & !is.null(parm_test)){

      # In some cases, when data is very noise, straight line fits
      # Here I test for 0 prediction
      test_zero <- round(mean(predict(output)),2)

      if ((test_zero) < 0.01){
        next()
      }

      temp_data$width_pred <- predict(output)

      extracted_a <- as.numeric(coef(output)[1])
      extracted_b <- as.numeric(coef(output)[2])
      extracted_k <- as.numeric(coef(output)[3])

      # Here we save the extracted parameters from Gompertz model
      list_extracted[[p2a]] <- i
      list_extracted_a[[p2a]] <- extracted_a
      list_extracted_b[[p2a]] <- extracted_b
      list_extracted_k[[p2a]] <- extracted_k
      p2a <- p2a + 1

      temp_data <- dplyr::arrange(temp_data, doy)

      # all what is below 0.1 goes to 0
      if (post_process == TRUE){
      temp_data$width_pred <- ifelse(temp_data$width_pred  < 0.01, 0,
                                   temp_data$width_pred)
      }

      temp_data$method <- "gompertz"

      list_temps[[b_holder]] <- temp_data
      b_holder = b_holder + 1

      if (fitted_save == TRUE){

        ggplot(temp_data, aes(x = doy, y = width_pred)) + geom_line() +
          geom_point(temp_data, mapping = aes(x = doy, y = width, alpha = note)) +
          ylab("width predicted") + theme_light() + guides(alpha = "none")

        ggsave(paste0("gom_", i, ".png"), width = 7, height = 6)
      }

    } else {

      list_errors[[p]] <- i
      p = p + 1

      if (search_initial_gom == TRUE){

        # print("Searching for initial Gompertz parameters... This might take some time")

      for (ii in 1:nrow(par_grid)){

        gom_a <- par_grid$gom_a[ii]
        gom_b <- par_grid$gom_b[ii]
        gom_k <- par_grid$gom_k[ii]

        capture.output(output <- try(nls(width ~ a*exp(-exp(b - k*doy)), data = temp_data,
                          start=list(a = gom_a, b = gom_b, k = gom_k),
                          nls.control(maxiter = 1000, tol = 1e-05,
                                      minFactor = 1/1024, printEval = FALSE,
                                      warnOnly = FALSE)), silent=TRUE))



        if (is(output, "nls")){

          # In some cases, when data is very noise, straight line fits
          # Here I test for 0 prediction
          test_zero <- round(mean(predict(output)),2)

          if ((test_zero) < 0.01){
            next()
          }

          extracted_a <- as.numeric(coef(output)[1])
          extracted_b <- as.numeric(coef(output)[2])
          extracted_k <- as.numeric(coef(output)[3])

          temp_data$width_pred <- predict(output)

          temp_data <- dplyr::arrange(temp_data, doy)

          temp_data$method <- "gompertz"
          temp_data$first_diff <- NULL

          list_temps[[b_holder]] <- temp_data
          b_holder = b_holder + 1

          if (fitted_save == TRUE){

            ggplot(temp_data, aes(x = doy, y = width_pred)) + geom_line() +
              geom_point(temp_data, mapping = aes(x = doy, y = width, alpha = note)) +
              ylab("width predicted") + theme_light() + guides(alpha = "none")

            ggsave(paste0("gom_", i, ".png"), width = 7, height = 6)
          }

          # Here we save the extracted parameters from Gompertz model
          list_extracted[[p2]] <- i
          list_extracted_a[[p2]] <- extracted_a
          list_extracted_b[[p2]] <- extracted_b
          list_extracted_k[[p2]] <- extracted_k

          # These are the initial input parameters
          list_solution_a[[p2]] <- gom_a
          list_solution_b[[p2]] <- gom_b
          list_solution_k[[p2]] <- gom_k
          list_solutions[[p2]] <- i
          p2 <- p2 + 1

          break()
        }
      }
    }
   }
  }

  errors_grid <- c(do.call(rbind, list_errors))
  solutions <- c(do.call(rbind, list_solutions))

  final_parameters <- data.frame(
    ID = c(do.call(rbind, list_solutions)),
    solution_a = c(do.call(rbind, list_solution_a)),
    solution_b = c(do.call(rbind, list_solution_b)),
    solution_k = c(do.call(rbind, list_solution_k))
  )

  extracted_parameters <- data.frame(
    ID = c(do.call(rbind, list_extracted)),
    a = c(do.call(rbind, list_extracted_a)),
    b = c(do.call(rbind, list_extracted_b)),
    k = c(do.call(rbind, list_extracted_k))
  )

  # remove solutions from errors grid
  errors_grid <- errors_grid[!(errors_grid %in% solutions)]

} else if (current_fitting_method == "double_gompertz"){

  # here we reduce the number of parameter combinations
  d_gom_a1_range <- d_gom_a1_range[sample(length(d_gom_a1_range), 25)]
  d_gom_a2_range <- d_gom_a2_range[sample(length(d_gom_a2_range), 25)]
  d_gom_b1_range <- d_gom_b1_range[sample(length(d_gom_b1_range), 25)]
  d_gom_b2_range <- d_gom_b2_range[sample(length(d_gom_b2_range), 25)]
  d_gom_k1_range <- d_gom_k1_range[sample(length(d_gom_k1_range), 25)]
  d_gom_k2_range <- d_gom_k2_range[sample(length(d_gom_k2_range), 25)]

  par_grid <- expand.grid(d_gom_a1 = d_gom_a1_range,
                          d_gom_a2 = d_gom_a2_range,
                          d_gom_b1 = d_gom_b1_range,
                          d_gom_b2 = d_gom_b2_range,
                          d_gom_k1 = d_gom_k1_range,
                          d_gom_k2 = d_gom_k2_range)

  # I select randomly 10,000 rows to reduce the computation time
  random_rows <- sample(nrow(par_grid), min(10000, nrow(par_grid)))
  par_grid <- par_grid[random_rows, ]


  for (i in unique_keys){

  temp_parameters <- parameters[parameters$key == i,]

  d_gom_a1 <- temp_parameters$d_gom_a1
  d_gom_a2 <- temp_parameters$d_gom_a2
  d_gom_b1 <- temp_parameters$d_gom_b1
  d_gom_b2 <- temp_parameters$d_gom_b2
  d_gom_k1 <- temp_parameters$d_gom_k1
  d_gom_k2 <- temp_parameters$d_gom_k2

  temp_data <- data_trees[data_trees$key == i,]

  # specify the measurement type - so we can distinguish from added zeros
  temp_data$note <- "raw measurement"

  # add zeros at the beginning
  if(add_zeros == TRUE){

    if (add_zeros_before == 'min'){

      min_doy <- min(temp_data$doy) - 1

    } else {

      if (!is.numeric(add_zeros_before) & add_zeros_before < 0){

        stop("The argument 'add_zeros_before' should be numeric and greater than 0")

      }

      min_doy <- add_zeros_before

    }

    row_list <- list()
    for (J in 1:min_doy){
      temp_row <- temp_data[1,]
      row_list[[J]] <- temp_row
    }

    new_rows <- do.call(rbind, row_list)
    new_rows$doy <- c(1:min_doy)
    new_rows$width <- 0
    new_rows$note <- "added zero"
    temp_data <- rbind(new_rows, temp_data)
  }

  capture.output(output <- try(nlsLM(y ~ double_gompertz(x, a1, b1, k1, a2, b2, k2),
                   start = c(a1 = d_gom_a1, a2 = d_gom_a2,
                             b1 = d_gom_b1, b2 = d_gom_b2,
                             k1 = d_gom_k1, k2 = d_gom_k2),
                   data = temp_data,
                   control = nls.lm.control(maxiter = 1000)), silent = TRUE) )

  # When all parameters are null, the class can still be nls. Additional check is needed
  parm_test <- c(d_gom_a1, d_gom_a2,
                 d_gom_b1, d_gom_b2,
                 d_gom_k1, d_gom_k2)

  if (is(output, "nls") & !is.null(parm_test)){

    # In some cases, when data is very noise, straight line fits
    # Here I test for 0 prediction
    test_zero <- round(mean(predict(output)),2)

    if ((test_zero) < 0.01){
      next()
    }

    temp_data$width_pred <- predict(output)

    extracted_a1 <- as.numeric(coef(output)[1])
    extracted_a2 <- as.numeric(coef(output)[2])
    extracted_b1 <- as.numeric(coef(output)[3])
    extracted_b2 <- as.numeric(coef(output)[4])
    extracted_k1 <- as.numeric(coef(output)[5])
    extracted_k2 <- as.numeric(coef(output)[6])

    # Here we save the extracted parameters from Gompertz model
    list_extracted_dg[[p3a]] <- i
    list_extracted_a1[[p3a]] <- extracted_a1
    list_extracted_b1[[p3a]] <- extracted_b1
    list_extracted_k1[[p3a]] <- extracted_k1
    list_extracted_a2[[p3a]] <- extracted_a2
    list_extracted_b2[[p3a]] <- extracted_b2
    list_extracted_k2[[p3a]] <- extracted_k2
    p3a <- p3a + 1

    temp_data <- dplyr::arrange(temp_data, doy)

    # all what is below 0.1 goes to 0
    if (post_process == TRUE){
      temp_data$width_pred <- ifelse(temp_data$width_pred  < 0.01, 0,
                                     temp_data$width_pred)
    }

    temp_data$method <- "double_gompertz"

    list_temps[[b_holder]] <- temp_data
    b_holder = b_holder + 1

    if (fitted_save == TRUE){

      ggplot(temp_data, aes(x = doy, y = width_pred)) + geom_line() +
        geom_point(temp_data, mapping = aes(x = doy, y = width, alpha = note)) +
        ylab("width predicted") + theme_light() + guides(alpha = "none")

      ggsave(paste0("d_gom_", i, ".png"), width = 7, height = 6)
    }

  } else {

    list_errors_dg[[p4]] <- i
    p4 = p4 + 1

    if (search_initial_double_gom == TRUE){

      # print("Searching for initial Gompertz parameters... This might take some time")

      for (ii in 1:nrow(par_grid)){

        d_gom_a1 <- par_grid$d_gom_a1[ii]
        d_gom_a2 <- par_grid$d_gom_a2[ii]
        d_gom_b1 <- par_grid$d_gom_b1[ii]
        d_gom_b2 <- par_grid$d_gom_b2[ii]
        d_gom_k1 <- par_grid$d_gom_k1[ii]
        d_gom_k2 <- par_grid$d_gom_k2[ii]

        capture.output(output <- try(minpack.lm::nlsLM(width ~ double_gompertz(doy, a1, b1, k1, a2, b2, k2),
                                        start = c(a1 = d_gom_a1, a2 = d_gom_a2,
                                                  b1 = d_gom_b1, b2 = d_gom_b2,
                                                  k1 = d_gom_k1, k2 = d_gom_k2),
                                        data = temp_data,
                                        control = minpack.lm::nls.lm.control(maxiter = 1000)), silent = TRUE) )

        if (is(output, "nls")){

          # In some cases, when data is very noise, straight line fits
          # Here I test for 0 prediction
          test_zero <- round(mean(predict(output)),2)

          if ((test_zero) < 0.01){
            next()
          }

          temp_data$width_pred <- predict(output)

          extracted_a1 <- as.numeric(coef(output)[1])
          extracted_a2 <- as.numeric(coef(output)[2])
          extracted_b1 <- as.numeric(coef(output)[3])
          extracted_b2 <- as.numeric(coef(output)[4])
          extracted_k1 <- as.numeric(coef(output)[5])
          extracted_k2 <- as.numeric(coef(output)[6])

          temp_data <- dplyr::arrange(temp_data, doy)

          if (post_process == TRUE){

            avg_fit <- mean(temp_data$width_pred)
            max_fit <- max(temp_data$width_pred)

            lagged_width_pred <- temp_data$width_pred[-length(temp_data$width_pred)]
            lagged_width_pred <- c(NA, lagged_width_pred)

            temp_data$first_diff <- temp_data$width_pred - lagged_width_pred
            temp_data[1, "first_diff"]  <- temp_data[2, "first_diff"]

            # In case of negative predictions
            if (any(temp_data$width_pred < 0)){

              shortcut1 <- temp_data
              shortcut1$neg_ind <- ifelse(shortcut1$width_pred < 0, TRUE, FALSE)

              th_doy <- as.numeric(max(shortcut1[shortcut1$neg_ind == TRUE, ][, "doy"]))

              temp_data$width_pred <- ifelse(temp_data$doy <= th_doy, 0,
                                             temp_data$width_pred)
            }

            temp_data$width_pred <- ifelse(temp_data$first_diff < 0 &
                                             (temp_data$width_pred < avg_fit), 0,
                                           ifelse(temp_data$width_pred < 0, 0,
                                                  temp_data$width_pred))

            test <- temp_data[temp_data$width_pred > avg_fit, ]
            test <- test[test$first_diff < 0, ]

            if (sum(test$first_diff < 0, na.rm = TRUE) > 0){

              for (J in 1:nrow(temp_data)){

                if (temp_data[J, "width_pred"] < 0.0001){

                  next()

                } else if (temp_data[J, "first_diff"] < 0){

                  temp_data[J, "width_pred"] <- temp_data[J - 1, "width_pred"]

                } else if (J != 1 && (temp_data[J, "width_pred"] - temp_data[J - 1, "width_pred"]) < 0){

                  temp_data[J, "width_pred"] <- temp_data[J - 1, "width_pred"]

                } else {

                  temp_data[J, "width_pred"] <- temp_data[J, "width_pred"]

                }

              }

            }
          }

          temp_data$method <- "double_gompertz"
          temp_data$first_diff <- NULL

          list_temps[[b_holder]] <- temp_data
          b_holder = b_holder + 1

          if (fitted_save == TRUE){

            ggplot(temp_data, aes(x = doy, y = width_pred)) + geom_line() +
              geom_point(temp_data, mapping = aes(x = doy, y = width, alpha = note)) +
              ylab("width predicted") + theme_light() + guides(alpha = "none")

            ggsave(paste0("d_gom_", i, ".png"), width = 7, height = 6)
          }

          list_solution_a1[[p3]] <- d_gom_a1
          list_solution_a2[[p3]] <- d_gom_a2
          list_solution_b1[[p3]] <- d_gom_b1
          list_solution_b2[[p3]] <- d_gom_b2
          list_solution_k1[[p3]] <- d_gom_k1
          list_solution_k2[[p3]] <- d_gom_k2
          list_solutions_dg[[p3]] <- i

          # Here we save the extracted parameters from Gompertz model
          list_extracted_dg[[p3]] <- i
          list_extracted_a1[[p3]] <- extracted_a1
          list_extracted_b1[[p3]] <- extracted_b1
          list_extracted_k1[[p3]] <- extracted_k1
          list_extracted_a2[[p3]] <- extracted_a2
          list_extracted_b2[[p3]] <- extracted_b2
          list_extracted_k2[[p3]] <- extracted_k2

          p3 <- p3 + 1

          break()
        }
      }
    }
  }
}

     errors_grid_dg <- c(do.call(rbind, list_errors_dg))
     solutions_dg <- c(do.call(rbind, list_solutions_dg))

     final_parameters_dg <- data.frame(
       solutions_double_gompertz = c(do.call(rbind, list_solutions_dg)),
       solution_a1 = c(do.call(rbind, list_solution_a1)),
       solution_a2 = c(do.call(rbind, list_solution_a2)),
       solution_b1 = c(do.call(rbind, list_solution_b1)),
       solution_b2 = c(do.call(rbind, list_solution_b2)),
       solution_k1 = c(do.call(rbind, list_solution_k1)),
       solution_k2 = c(do.call(rbind, list_solution_k2))
     )

     extracted_parameters_dg <- data.frame(
       ID = c(do.call(rbind, list_extracted_dg)),
       a1 = c(do.call(rbind, list_extracted_a1)),
       b1 = c(do.call(rbind, list_extracted_b1)),
       k1 = c(do.call(rbind, list_extracted_k1)),

       a2 = c(do.call(rbind, list_extracted_a2)),
       b2 = c(do.call(rbind, list_extracted_b2)),
       k2 = c(do.call(rbind, list_extracted_k2))
     )


     # remove solutions from errors grid
     errors_grid_dg <- errors_grid_dg[!(errors_grid_dg %in% solutions_dg)]

  } else if (current_fitting_method == "brnn"){

    data_neurons <- parameters[,c("key", "brnn_neurons")][!is.na(parameters$key),]

    for (i in unique_keys){

    temp_data <- data_trees[data_trees$key == i,]
    temp_neurons <- data_neurons[data_neurons$key == i, ]

    # specify the measurement type - so we can distinguish from added zeros
    temp_data$note <- "raw measurement"

      # add zeros at the beginning
    if (add_zeros == TRUE){

      if (add_zeros_before == 'min'){

        min_doy <- min(temp_data$doy) - 1

      } else {

        if (!is.numeric(add_zeros_before) & add_zeros_before < 0){

          stop("The argument 'add_zeros_before' should be numeric and greater than 0")

        }

        min_doy <- add_zeros_before

      }

    row_list <- list()

      for (J in 1:min_doy){
        temp_row <- temp_data[1,]
        row_list[[J]] <- temp_row
      }
      new_rows <- do.call(rbind, row_list)
      new_rows$doy <- c(1:min_doy)
      new_rows$width <- 0
      new_rows$note <- "added zero"
      temp_data <- rbind(new_rows, temp_data)

    }

      capture.output(output <- try(brnn(width ~ doy, data = temp_data,
                                        neurons = temp_neurons$brnn_neurons[1]), silent=TRUE))

      temp_data$width_pred <- predict(output)

      temp_data <- dplyr::arrange(temp_data, doy)

      if (post_process == TRUE){

      avg_fit <- mean(temp_data$width_pred)
      max_fit <- max(temp_data$width_pred)

      lagged_width_pred <- temp_data$width_pred[-length(temp_data$width_pred)]
      lagged_width_pred <- c(NA, lagged_width_pred)

      temp_data$first_diff <- temp_data$width_pred - lagged_width_pred
      temp_data[1, "first_diff"]  <- temp_data[2, "first_diff"]

      # In case of negative predictions
      if (any(temp_data$width_pred < 0)){

        shortcut1 <- temp_data
        shortcut1$neg_ind <- ifelse(shortcut1$width_pred < 0, TRUE, FALSE)

        th_doy <- as.numeric(max(shortcut1[shortcut1$neg_ind == TRUE, ][, "doy"]))

        temp_data$width_pred <- ifelse(temp_data$doy <= th_doy, 0,
                                     temp_data$width_pred)
        }

        temp_data$width_pred <- ifelse(temp_data$first_diff < 0 &
                                       (temp_data$width_pred < avg_fit), 0,
                              ifelse(temp_data$width_pred < 0, 0,
                                     temp_data$width_pred))

        test <- temp_data[temp_data$width_pred > avg_fit, ]
        test <- test[test$first_diff < 0, ]

      if (sum(test$first_diff < 0, na.rm = TRUE) > 0){

        for (J in 1:nrow(temp_data)){

          if (temp_data[J, "width_pred"] < 0.0001){

            next()

          } else if (temp_data[J, "first_diff"] < 0){

            temp_data[J, "width_pred"] <- temp_data[J - 1, "width_pred"]

          } else if (J != 1 && (temp_data[J, "width_pred"] - temp_data[J - 1, "width_pred"]) < 0){

            temp_data[J, "width_pred"] <- temp_data[J - 1, "width_pred"]

          } else {

            temp_data[J, "width_pred"] <- temp_data[J, "width_pred"]

          }

        }

      }
    }

      temp_data$first_diff <- NULL
      temp_data$method <- "brnn"

      list_temps[[b_holder]] <- temp_data
      b_holder = b_holder + 1

      if (fitted_save == TRUE){

        ggplot(temp_data, aes(x = doy, y = width_pred)) + geom_line() +
          geom_point(temp_data, mapping = aes(x = doy, y = width, alpha = note)) +
          ylab("width predicted") + theme_light() + guides(alpha = "none")

        ggsave(paste0("brnn_", i, ".png"), width = 7, height = 6)
      }

    }

  } else if (current_fitting_method == "gam"){

    for (i in unique_keys){

      temp_data <- data_trees[data_trees$key == i, ]

      temp_parameters <- parameters[parameters$key == i, ]

      gam_k <- temp_parameters$gam_k
      gam_sp <- temp_parameters$gam_sp

      # specify the measurement type - so we can distinguish from added zeros
      temp_data$note <- "raw measurement"

      # add zeros at the beggining
    if(add_zeros == TRUE){

      if (add_zeros_before == 'min'){

        min_doy <- min(temp_data$doy) - 1

      } else {

        if (!is.numeric(add_zeros_before) & add_zeros_before < 0){

          stop("The argument 'add_zeros_before' should be numeric and greater than 0")

        }

        min_doy <- add_zeros_before

      }

      row_list <- list()
      for (J in 1:min_doy){
        temp_row <- temp_data[1,]
        row_list[[J]] <- temp_row
      }
      new_rows <- do.call(rbind, row_list)
      new_rows$doy <- c(1:min_doy)
      new_rows$width <- 0
      new_rows$note <- "added zero"

      temp_data <- rbind(new_rows, temp_data)
      temp_data$width <- as.numeric(temp_data$width)
    }

      output <- gam(width ~ s(doy, k = gam_k[1], bs ="ps", sp = gam_sp[1]),
                    data = temp_data, method = "REML")

      temp_data$width_pred <- predict(output)

      temp_data <- dplyr::arrange(temp_data, doy)

      if (post_process == TRUE){

      avg_fit <- ifelse(
        add_zeros == TRUE, mean(temp_data$width_pred),
        mean(temp_data$width_pred)/2 ) # That is a rule of thumb, but works nice for all the data


      max_fit <- max(temp_data$width_pred)

      lagged_width_pred <- temp_data$width_pred[-length(temp_data$width_pred)]
      lagged_width_pred <- c(NA, lagged_width_pred)

      temp_data$first_diff <- temp_data$width_pred - lagged_width_pred
      temp_data[1, "first_diff"]  <- temp_data[2, "first_diff"]

      # In case of negative predictions
      if (sum(temp_data$width_pred < 0, na.rm = TRUE) > 0){

        shortcut1 <- temp_data
        shortcut1$neg_ind <- ifelse(shortcut1$width_pred < 0, TRUE, FALSE)

        th_doy <- as.numeric(max(shortcut1[shortcut1$neg_ind == TRUE, ][, "doy"]))
        temp_data$width_pred <- ifelse(temp_data$doy <= th_doy, 0,
                                     temp_data$width_pred)

        }

      temp_data$width_pred <- ifelse(temp_data$first_diff < 0 &
                                     (temp_data$width_pred < avg_fit), 0,
                                   ifelse(temp_data$width_pred < 0, 0,
                                          temp_data$width_pred))

      test <- temp_data[temp_data$width_pred > avg_fit, ]
      test <- test[test$first_diff < 0, ]

      if (sum(test$first_diff < 0, na.rm = TRUE) > 0){}

      for (J in 1:nrow(temp_data)){

        J_shift <- ifelse(J == 1, 0, 1)

        if (temp_data[J, "width_pred"] < 0.0001){

          next()

        } else if (temp_data[J, "first_diff"] < 0){

          temp_data[J, "width_pred"] <- temp_data[J - J_shift, "width_pred"]

        } else if ((temp_data[J, "width_pred"] - temp_data[J - J_shift, "width_pred"]) < 0){

          temp_data[J, "width_pred"] <- temp_data[J - J_shift, "width_pred"]

        } else {

          temp_data[J, "width_pred"] <- temp_data[J, "width_pred"]

          }

        }

       # plot(y = temp_data$width, x = temp_data$doy, main = i)
       # lines(y = temp_data$width_pred, x = temp_data$doy, type = "l")

      }

      temp_data$method <- "gam"
      temp_data$first_diff <- NULL
      list_temps[[b_holder]] <- temp_data
      b_holder = b_holder + 1

      if (fitted_save == TRUE){

        ggplot(temp_data, aes(x = doy, y = width_pred)) + geom_line() +
          geom_point(temp_data, mapping = aes(x = doy, y = width, alpha = note)) +
          ylab("width predicted") + theme_light() + guides(alpha = "none")

        ggsave(paste0("gam_", i, ".png"), width = 7, height = 6)
      }
    }

  } else {

      stop("Select proper fitting_method!")
  }

     pbar_holder = pbar_holder + 1

} # end of ut


  if (is.null(do.call(rbind, list_temps))){



    if ("double_gompertz" %in% fitting_method & search_initial_double_gom == FALSE){

      stop("all trials to fit the data failed. \ntry to set the parameter search_initial_double_gom = TRUE ")

    }

    if ("gompertz" %in% fitting_method & search_initial_gom == FALSE){

      stop("all trials to fit the data failed. \ntry to set the parameter search_initial_gom = TRUE ")

    }


  }


  output_list <- list(fitted = do.call(rbind, list_temps),
                      gompertz_initial_parameters = final_parameters,
                      gompertz_initial_parameters_errors = errors_grid,
                      gompertz_model_parameters = extracted_parameters,
                      double_gompertz_initial_parameters = final_parameters_dg,
                      double_gompertz_initial_parameters_errors = errors_grid_dg,
                      double_gompertz_model_parameters = extracted_parameters_dg
                      )

  class(output_list) <- "xpsg"

  return(output_list)

}






