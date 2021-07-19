#' XPSgrowth
#'
#' XylemPhloemSeasonalGrowth: This Function fits and compares the selected
#' methods for modeling seasonal xylem and phloem data.
#' @param data_trees a data frame with ID variables and wood formation data with
#' columns DOY and Width
#' @param parameters a data frame with ID variables and initial parameter values
#' for the selected methods
#' @param search_initial_gom logical, should the algorithm to search initial
#' Gompertz parameters be applied?
#' @param fitting_method vector of one or more methods to be compared:
#' "gompertz", "gam", "brnn"
#' @param ID_vars character vector of variables which indicate column names of
#' ID variables
#' @param fitted_save logical, should the fitted curves be saved in current
#' working directory?
#' @param add_zeros logical, should zero observations at the beginning of
#' growing season be added?
#' @param add_zeros_before if 'min' (character) then zeros will be added prior
#' to the first observation in each year. Alternatively, users can specify
#' absolute DOY prior which zeros will be added.
#' @param post_process logical, should the post-process algorithm be applied?
#' @param unified_parameters logical, if FALSE, the algorithm will use only
#' manually selected function parameters. See the arguments 'gom_a', 'gom_b',
#' 'gom_k', 'brnn_neurons', 'gam_k' and 'gam_sp'. Default is FALSE
#' @param gom_a numeric, the parameter a for the Gompertz function
#' @param gom_b numeric, the parameter b for the Gompertz function
#' @param gom_k numeric, the parameter k for the Gompertz function
#' @param brnn_neurons positive integer, the number of neurons to be used by
#' the BRNN method
#' @param gam_k numeric, the parameter k for General Additive Model (GAM)
#' @param gam_sp numeric, the parameter sp for General Additive Model (GAM)
#'
#' @return a list with the following elements:
#' \enumerate{
#'  \item $fitted - a data frame with fitted wood formation data
#'  \item $gompertz_grid_search - a data frame with selected initial parameter values
#'  \item $gompertz_grid_search_errors - a data frame with unsuccessful cases of gompertz grid search
#'}
#'
#' @export
#'
#' @examples
#' library(rTG)
#'
#' # Load data
#' data(parameters)
#' data(data_trees)
#' simulation_1 <- XPSgrowth(data_trees = data_trees,
#'      parameters = parameters,
#'      ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
#'      fitting_method = c("gompertz", "gam", "brnn"),
#'      fitted_save = FALSE,
#'      search_initial_gom = TRUE,
#'      add_zeros = TRUE,
#'      post_process = TRUE)

XPSgrowth <- function(data_trees, parameters = NULL,
                 search_initial_gom = FALSE,
                 fitting_method = c("gompertz", "gam", "brnn"),
                 ID_vars = NULL,
                 fitted_save = FALSE,
                 add_zeros = TRUE,
                 add_zeros_before = 'min',
                 post_process = TRUE,
                 unified_parameters = FALSE,
                 gom_a = NA, gom_b = NA, gom_k = NA,
                 brnn_neurons = 3,
                 gam_k = 8, gam_sp = 0.5){

  # Defining global variables
  Width <- NULL
  DOY <- NULL
  key <- NULL
  first_diff <- NULL
  Width_pred <- NULL
  neg_ind <- NULL
  errors_grid <- NA
  txtProgressBar <- NULL
  setTxtProgressBar <- NULL
  # add_zeros_before <- NULL

  # Progress bar
  pb <- txtProgressBar(min = 0, max = length(fitting_method), style = 3)

  # 1 are all ID_vars in both tables?

  if(unified_parameters == FALSE){

    if (sum(ID_vars %in% colnames(parameters) == FALSE)>0){
      stop("check your ID_vars, they don't exist in your parameters")
    }
  }

  if (sum(ID_vars %in% colnames(data_trees) == FALSE)>0){
    stop("check your ID_vars, they don't exist in your data_treses")
  }

  if (sum(is.na(data_trees)) > 0){
    stop("Please remove missing values from data_trees")
  }

  # Columns DOY & Width must be present in data_trees
  if (!("DOY" %in% colnames(data_trees))){
    stop("Column 'DOY' is missing in data_trees")
  }

  if (!("Width" %in% colnames(data_trees))){
    stop("Column 'Width' is missing in data_trees")
  }

  # Create data key
  data_trees$key <- apply(data_trees[, ID_vars], 1,
                          paste0, sep = "", collapse = "_")

  data_trees$key <- gsub(" ", "", data_trees$key) # remove empty characters

  if (unified_parameters == TRUE){

    parameters <- data.frame(key = unique(data_trees$key))

    parameters$gom_a <- gom_a
    parameters$gom_b <- gom_b
    parameters$gom_k <- gom_k
    parameters$brnn_neurons <- brnn_neurons
    parameters$gam_k <- gam_k
    parameters$gam_sp <- gam_sp

  } else {

    parameters$key <- apply(parameters[, ID_vars], 1,
                            paste0, sep = "", collapse = "_")

    parameters$key <- gsub(" ", "", parameters$key) # remove empty characters

  }

  list_temps <- list()
  list_errors <- list()
  list_solutions <- list()
  list_solution_a <- list()
  list_solution_b <- list()
  list_solution_k <- list()

  b_holder = 1
  p = 1
  p2 = 1
  pbar_holder = 1

  unique_keys <- unique(data_trees$key)


for (ut in fitting_method){

  setTxtProgressBar(pb, pbar_holder)

     current_fitting_method = ut

if (current_fitting_method == "gompertz"){

  for (i in unique_keys){

    gom_a_range <- 3000
    gom_b_range <- seq(50, 1000, 50)
    gom_k_range <- seq(1, 500, 2)

    par_grid <- expand.grid(gom_a = gom_a_range,
                            gom_b = gom_b_range,
                            gom_k = gom_k_range)

    temp_parameters <- parameters[parameters$key == i,]

    gom_a <- temp_parameters$gom_a
    gom_b <- temp_parameters$gom_b
    gom_k <- temp_parameters$gom_k

    temp_data <- data_trees[data_trees$key == i,]

    # add zeros at the beginning
    if(add_zeros == TRUE){

      if (add_zeros_before == 'min'){

        min_DOY <- min(temp_data$DOY)

      } else {

        if (!is.numeric(add_zeros_before) & add_zeros_before < 0){

          stop("The argument 'add_zeros_before' should be numeric and greater than 0")

        }

        min_DOY <- add_zeros_before

      }


      row_list <- list()
      for (J in 1:min_DOY){
        temp_row <- temp_data[1,]
        row_list[[J]] <- temp_row
      }

      new_rows <- do.call(rbind, row_list)
      new_rows$DOY <- c(1:min_DOY)
      new_rows$Width <- 0
      temp_data <- rbind(new_rows, temp_data)
    }

    output <- try(nls(Width ~ a*exp(-exp(b - k*DOY)), data = temp_data,
                      start=list(a = gom_a, b = gom_b, k = gom_k),
                      nls.control(maxiter = 1000, tol = 1e-05,
                                  minFactor = 1/1024, printEval = FALSE,
                                  warnOnly = FALSE)), silent = TRUE)

    if (class(output) == "nls"){

      temp_data$Width_pred <- predict(output)

      # all what is below 0.1 goes to 0
      if (post_process == TRUE){
      temp_data$Width_pred <- ifelse(temp_data$Width_pred  < 0.01, 0,
                                   temp_data$Width_pred)
      }

      temp_data$method <- "gompertz"

      list_temps[[b_holder]] <- temp_data
      b_holder = b_holder + 1

      if (fitted_save == TRUE){
        ggplot(temp_data, aes(x = DOY, y = Width_pred)) + geom_line() +
          geom_point(temp_data, mapping = aes(x = DOY, y = Width)) +
          ylab("Width predicted") + theme_light()
        ggsave(paste0("gom_", i, ".png"))
      }

    } else {

      list_errors[[p]] <- i
      p = p + 1

      if (search_initial_gom == TRUE){

      for (ii in 1:nrow(par_grid)){

        gom_a <- par_grid$gom_a[ii]
        gom_b <- par_grid$gom_b[ii]
        gom_k <- par_grid$gom_k[ii]

        output <- try(nls(Width ~ a*exp(-exp(b - k*DOY)), data = temp_data,
                          start=list(a = gom_a, b = gom_b, k = gom_k),
                          nls.control(maxiter = 1000, tol = 1e-05,
                                      minFactor = 1/1024, printEval = FALSE,
                                      warnOnly = FALSE)))


        if (class(output) == "nls"){

          temp_data$Width_pred <- predict(output)

          temp_data$method <- "gompertz"
          temp_data$first_diff <- NULL

          list_temps[[b_holder]] <- temp_data
          b_holder = b_holder + 1

          if (fitted_save == TRUE){
            ggplot(temp_data, aes(x = DOY, y = Width_pred)) + geom_line() +
              geom_point(temp_data, mapping = aes(x = DOY, y = Width)) +
              ylab("Width predicted") + theme_light()
            ggsave(paste0("gom_", i, ".png"))
          }

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
    solutions = c(do.call(rbind, list_solutions)),
    solution_a = c(do.call(rbind, list_solution_a)),
    solution_b = c(do.call(rbind, list_solution_b)),
    solution_k = c(do.call(rbind, list_solution_k))
  )

  errors <- errors[!(errors %in% solutions)]

  } else if (current_fitting_method == "brnn"){

    data_neurons <- parameters[,c("key", "brnn_neurons")][!is.na(parameters$key),]

    for (i in unique_keys){

    temp_data <- data_trees[data_trees$key == i,]
    temp_neurons <- data_neurons[data_neurons$key == i, ]

      # add zeros at the beginning
    if (add_zeros == TRUE){

      if (add_zeros_before == 'min'){

        min_DOY <- min(temp_data$DOY)

      } else {

        if (!is.numeric(add_zeros_before) & add_zeros_before < 0){

          stop("The argument 'add_zeros_before' should be numeric and greater than 0")

        }

        min_DOY <- add_zeros_before

      }

    row_list <- list()

      for (J in 1:min_DOY){
        temp_row <- temp_data[1,]
        row_list[[J]] <- temp_row
      }
      new_rows <- do.call(rbind, row_list)
      new_rows$DOY <- c(1:min_DOY)
      new_rows$Width <- 0
      temp_data <- rbind(new_rows, temp_data)

    }

      capture.output(output <- brnn(Width ~ DOY, data = temp_data,
                     neurons = temp_neurons$brnn_neurons))

      temp_data$Width_pred <- predict(output)

      # plot(y = temp_data$Width, x = temp_data$DOY)
      # lines(y = temp_data$Width_pred, x = temp_data$DOY, type = "l")

      if (post_process == TRUE){

      avg_fit <- mean(temp_data$Width_pred)
      max_fit <- max(temp_data$Width_pred)

      lagged_Width_pred <- temp_data$Width_pred[-length(temp_data$Width_pred)]
      lagged_Width_pred <- c(NA, lagged_Width_pred)

      temp_data$first_diff <- temp_data$Width_pred - lagged_Width_pred
      temp_data[1, "first_diff"]  <- temp_data[2, "first_diff"]

      # In case of negative predictions
      if (any(temp_data$Width_pred < 0)){

        shortcut1 <- temp_data
        shortcut1$neg_ind <- ifelse(shortcut1$Width_pred < 0, TRUE, FALSE)

        th_DOY <- as.numeric(max(shortcut1[shortcut1$neg_ind == TRUE, ][, "DOY"]))

        temp_data$Width_pred <- ifelse(temp_data$DOY <= th_DOY, 0,
                                     temp_data$Width_pred)
        }

        temp_data$Width_pred <- ifelse(temp_data$first_diff < 0 &
                                       (temp_data$Width_pred < avg_fit), 0,
                              ifelse(temp_data$Width_pred < 0, 0,
                                     temp_data$Width_pred))

        test <- temp_data[temp_data$Width_pred > avg_fit, ]
        test <- test[test$first_diff < 0, ]

      if (sum(test$first_diff < 0, na.rm = TRUE) > 0){

        for (J in 1:nrow(temp_data)){


          if (temp_data[J, "Width_pred"] < 0.0001){

            next()

          } else if (temp_data[J, "first_diff"] < 0){

            temp_data[J, "Width_pred"] <- temp_data[J - 1, "Width_pred"]

          } else if ((temp_data[J, "Width_pred"] - temp_data[J - 1, "Width_pred"]) < 0){

            temp_data[J, "Width_pred"] <- temp_data[J - 1, "Width_pred"]

          } else {

            temp_data[J, "Width_pred"] <- temp_data[J, "Width_pred"]

          }

        }

      }
    }
      #plot(y = temp_data$Width, x = temp_data$DOY, type = "p")
      #lines(y = temp_data$Width_pred, x = temp_data$DOY)

      errors <- NA
      final_parameters <- NA

      temp_data$first_diff <- NULL
      temp_data$method <- "brnn"

      list_temps[[b_holder]] <- temp_data
      b_holder = b_holder + 1

      if (fitted_save == TRUE){
        ggplot(temp_data, aes(x = DOY, y = Width_pred)) + geom_line() +
          geom_point(temp_data, mapping = aes(x = DOY, y = Width)) +
          ylab("Width predicted") + theme_light()
        ggsave(paste0("brnn_", i, ".png"))
      }

    }
  } else if (current_fitting_method == "gam"){

    for (i in unique_keys){

      temp_data <- data_trees[data_trees$key == i, ]

      temp_parameters <- parameters[parameters$key == i, ]

      gam_k <- temp_parameters$gam_k
      gam_sp <- temp_parameters$gam_sp

      # add zeros at the beggining
    if(add_zeros == TRUE){

      if (add_zeros_before == 'min'){

        min_DOY <- min(temp_data$DOY)

      } else {

        if (!is.numeric(add_zeros_before) & add_zeros_before < 0){

          stop("The argument 'add_zeros_before' should be numeric and greater than 0")

        }

        min_DOY <- add_zeros_before

      }

      row_list <- list()
      for (J in 1:min_DOY){
        temp_row <- temp_data[1,]
        row_list[[J]] <- temp_row
      }
      new_rows <- do.call(rbind, row_list)
      new_rows$DOY <- c(1:min_DOY)
      new_rows$Width <- 0

      temp_data <- rbind(new_rows, temp_data)
      temp_data$Width <- as.numeric(temp_data$Width)
    }

      output <- gam(Width ~ s(DOY, k = gam_k, bs ="ps", sp = gam_sp),
                    data = temp_data, method = "REML")

      temp_data$Width_pred <- predict(output)

      # plot(y = temp_data$Width, x = temp_data$DOY, main = i)
      # lines(y = temp_data$Width_pred, x = temp_data$DOY, type = "l")

      if (post_process == TRUE){

      avg_fit <- mean(temp_data$Width_pred)
      max_fit <- max(temp_data$Width_pred)

      lagged_Width_pred <- temp_data$Width_pred[-length(temp_data$Width_pred)]
      lagged_Width_pred <- c(NA, lagged_Width_pred)

      temp_data$first_diff <- temp_data$Width_pred - lagged_Width_pred
      temp_data[1, "first_diff"]  <- temp_data[2, "first_diff"]

      # In case of negative predictions
      if (sum(temp_data$Width_pred < 0, na.rm = TRUE) > 0){

        shortcut1 <- temp_data
        shortcut1$neg_ind <- ifelse(shortcut1$Width_pred < 0, TRUE, FALSE)

        th_DOY <- as.numeric(max(shortcut1[shortcut1$neg_ind == TRUE, ][, "DOY"]))
        temp_data$Width_pred <- ifelse(temp_data$DOY <= th_DOY, 0,
                                     temp_data$Width_pred)

        }

      temp_data$Width_pred <- ifelse(temp_data$first_diff < 0 &
                                     (temp_data$Width_pred < avg_fit), 0,
                                   ifelse(temp_data$Width_pred < 0, 0,
                                          temp_data$Width_pred))

      test <- temp_data[temp_data$Width_pred > avg_fit, ]
      test <- test[test$first_diff < 0, ]

      if (sum(test$first_diff < 0, na.rm = TRUE) > 0){}

      for (J in 1:nrow(temp_data)){


        if (temp_data[J, "Width_pred"] < 0.0001){

          next()

        } else if (temp_data[J, "first_diff"] < 0){

          temp_data[J, "Width_pred"] <- temp_data[J - 1, "Width_pred"]

        } else if ((temp_data[J, "Width_pred"] - temp_data[J - 1, "Width_pred"]) < 0){

          temp_data[J, "Width_pred"] <- temp_data[J - 1, "Width_pred"]

        } else {

          temp_data[J, "Width_pred"] <- temp_data[J, "Width_pred"]

          }

        }

     # plot(y = temp_data$Width, x = temp_data$DOY, main = i)
     # lines(y = temp_data$Width_pred, x = temp_data$DOY, type = "l")

      }

      errors <- NA
      final_parameters <- NA

      temp_data$method <- "gam"
      temp_data$first_diff <- NULL
      list_temps[[b_holder]] <- temp_data
      b_holder = b_holder + 1

      if (fitted_save == TRUE){
        ggplot(temp_data, aes(x = DOY, y = Width_pred)) +
          geom_line() + geom_point(temp_data, mapping = aes(x = DOY, y = Width)) +
          ylab("Width predicted") + theme_light()
        ggsave(paste0("gam_", i, ".png"))
      }
    }

  } else {

      stop("Select proper fitting_method!")
  }

     pbar_holder = pbar_holder + 1

} # end of ut

  output_list <- list(fitted = do.call(rbind, list_temps),
                      gompertz_grid_search = final_parameters,
                      gompertz_grid_search_errors = errors_grid

                      )

  class(output_list) <- "xpsg"

  return(output_list)

}







