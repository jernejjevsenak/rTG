#' XPSgrowth
#'
#' XylemPhloemSeasonalGrowth: This Function fits and compares the selected methods for modeling seasonal xylem and phloem data.
#' @param data_trees a data frame with ID variables and wood formation data with columns DY and EWM
#' @param data_site a data frame with ID variables and initial parameter values for the selected methods
#' @param search_initial_gom logical, should the algorithm to search initial Gompertz parameters be applied?
#' @param fitting_method vector of one or more methods to be compared: "gompertz", "gam", "brnn"
#' @param ID_vars character vector of variables which indicate column names of ID variables
#' @param fitted_save logical, should the fitted curves be saved in current working directory?
#' @param add_zeros logical, should zero observations at the beginning of growing season be added?
#' @param post_process logical, should the post-process algorithm be applied?
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
#' \dontrun{
#' library(rTG)
#'
#' # Load data
#' data(data_site)
#' data(data_trees)
#' simulation_1 <- XPSgrowth(data_trees = data_trees,
#'      data_site = data_site,
#'      ID_vars = c("Species", "Tissue", "Site", "Year", "Tree"),
#'      fitting_method = c("gompertz", "gam", "brnn"),
#'      fitted_save = TRUE,
#'      search_initial_gom = TRUE,
#'      add_zeros = TRUE,
#'      post_process = TRUE)
#'
#' }

XPSgrowth <- function(data_trees, data_site,
                 search_initial_gom = FALSE,
                 fitting_method = c("gompertz", "gam", "brnn"),
                 ID_vars = NULL,
                 fitted_save = FALSE,
                 add_zeros = TRUE,
                 post_process = TRUE){

  # Defining global variables
  EWM <- NULL
  DY <- NULL
  key <- NULL
  first_diff <- NULL
  EWM_pred <- NULL
  brnn_neurons <- NULL
  neg_ind <- NULL
  errors_grid <- NA

  # 1 are all ID_vars in both tables?
  if (sum(ID_vars %in% colnames(data_site) == FALSE)>0){

    stop("check you ID_vars, they don't exist in your data_site")
  }

  if (sum(ID_vars %in% colnames(data_trees) == FALSE)>0){

    stop("check you ID_vars, they don't exist in your data_treses")
  }

  # Create data key
  data_site$key <- apply(data_site[, ID_vars], 1,
                         paste0, sep = "", collapse = "_")

  data_trees$key <- apply(data_trees[, ID_vars], 1,
                          paste0, sep = "", collapse = "_")

  list_temps <- list()
  list_errors <- list()
  list_solutions <- list()
  list_solution_a <- list()
  list_solution_b <- list()
  list_solution_k <- list()

  b_holder = 1
  p = 1
  p2 = 1

  unique_keys <- unique(data_trees$key)


for (ut in fitting_method){

     current_fitting_method = ut

if (current_fitting_method == "gompertz"){

  for (i in unique_keys){

    gom_a_range <- 3000
    gom_b_range <- seq(50, 1000, 50)
    gom_k_range <- seq(1, 500, 2)

    par_grid <- expand.grid(gom_a = gom_a_range,
                            gom_b = gom_b_range,
                            gom_k = gom_k_range)

    parameters <- data_site[data_site$key == i,]

    gom_a <- parameters$gom_a
    gom_b <- parameters$gom_b
    gom_k <- parameters$gom_k

    temp_data <- data_trees[data_trees$key == i,]

    # add zeros at the beginning
    if(add_zeros == TRUE){

      min_DY <- min(temp_data$DY)
      row_list <- list()
      for (J in 1:min_DY){
        temp_row <- temp_data[1,]
        row_list[[J]] <- temp_row
      }

      new_rows <- do.call(rbind, row_list)
      new_rows$DY <- c(1:min_DY)
      new_rows$EWM <- 0
      temp_data <- rbind(new_rows, temp_data)
    }

    output <- try(nls(EWM ~ a*exp(-exp(b - k*DY)), data = temp_data,
                      start=list(a = gom_a, b = gom_b, k = gom_k),
                      nls.control(maxiter = 1000, tol = 1e-05,
                                  minFactor = 1/1024, printEval = FALSE,
                                  warnOnly = FALSE)), silent = TRUE)

    if (class(output) == "nls"){

      temp_data$EWM_pred <- predict(output)

      # all what is below 0.1 goes to 0
      if (post_process == TRUE){
      temp_data$EWM_pred <- ifelse(temp_data$EWM_pred  < 0.01, 0,
                                   temp_data$EWM_pred)
      }

      temp_data$method <- "gompertz"

      list_temps[[b_holder]] <- temp_data
      b_holder = b_holder + 1

      if (fitted_save == TRUE){
        ggplot(temp_data, aes(x = DY, y = EWM_pred)) + geom_line() +
          geom_point(temp_data, mapping = aes(x = DY, y = EWM)) +
          ylab("EWM predicted") + theme_light()
        ggsave(paste0("gom_", i, ".png"))
      }

    } else {

      list_errors[[p]] <- i
      p = p + 1

      if (search_initial_gom == TRUE){

      for (ii in 1:nrow(par_grid)){

        print(ii)
        gom_a <- par_grid$gom_a[ii]
        gom_b <- par_grid$gom_b[ii]
        gom_k <- par_grid$gom_k[ii]

        output <- try(nls(EWM ~ a*exp(-exp(b - k*DY)), data = temp_data,
                          start=list(a = gom_a, b = gom_b, k = gom_k),
                          nls.control(maxiter = 1000, tol = 1e-05,
                                      minFactor = 1/1024, printEval = FALSE,
                                      warnOnly = FALSE)))


        if (class(output) == "nls"){

          temp_data$EWM_pred <- predict(output)

          temp_data$method <- "gompertz"
          temp_data$first_diff <- NULL

          list_temps[[b_holder]] <- temp_data
          b_holder = b_holder + 1

          if (fitted_save == TRUE){
            ggplot(temp_data, aes(x = DY, y = EWM_pred)) + geom_line() +
              geom_point(temp_data, mapping = aes(x = DY, y = EWM)) +
              ylab("EWM predicted") + theme_light()
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

    data_neurons <- data_site[,c("key", "brnn_neurons")][!is.na(data_site$key),]

    for (i in unique_keys){

    temp_data <- data_trees[data_trees$key == i,]
    temp_neurons <- data_neurons[data_neurons$key == i, ]

      # add zeros at the beginning
    if (add_zeros == TRUE){
    min_DY <- min(temp_data$DY)
      row_list <- list()

      for (J in 1:min_DY){
        temp_row <- temp_data[1,]
        row_list[[J]] <- temp_row
      }
      new_rows <- do.call(rbind, row_list)
      new_rows$DY <- c(1:min_DY)
      new_rows$EWM <- 0
      temp_data <- rbind(new_rows, temp_data)

    }

      output <- brnn(EWM ~ DY, data = temp_data,
                     neurons = temp_neurons$brnn_neurons)

      temp_data$EWM_pred <- predict(output)

      # plot(y = temp_data$EWM, x = temp_data$DY)
      # lines(y = temp_data$EWM_pred, x = temp_data$DY, type = "l")

      if (post_process == TRUE){

      avg_fit <- mean(temp_data$EWM_pred)
      max_fit <- max(temp_data$EWM_pred)

      lagged_EWM_pred <- temp_data$EWM_pred[-length(temp_data$EWM_pred)]
      lagged_EWM_pred <- c(NA, lagged_EWM_pred)

      temp_data$first_diff <- temp_data$EWM_pred - lagged_EWM_pred
      temp_data[1, "first_diff"]  <- temp_data[2, "first_diff"]

      # In case of negative predictions
      if (any(temp_data$EWM_pred < 0)){

        shortcut1 <- temp_data
        shortcut1$neg_ind <- ifelse(shortcut1$EWM_pred < 0, TRUE, FALSE)

        th_DY <- as.numeric(max(shortcut1[shortcut1$neg_ind == TRUE, ][, "DY"]))

        temp_data$EWM_pred <- ifelse(temp_data$DY <= th_DY, 0,
                                     temp_data$EWM_pred)
        }

        temp_data$EWM_pred <- ifelse(temp_data$first_diff < 0 &
                                       (temp_data$EWM_pred < avg_fit), 0,
                              ifelse(temp_data$EWM_pred < 0, 0,
                                     temp_data$EWM_pred))

        test <- temp_data[temp_data$EWM_pred > avg_fit, ]
        test <- test[test$first_diff < 0, ]

      if (sum(test$first_diff < 0, na.rm = TRUE) > 0){

        for (J in 1:nrow(temp_data)){


          if (temp_data[J, "EWM_pred"] < 0.0001){

            next()

          } else if (temp_data[J, "first_diff"] < 0){

            temp_data[J, "EWM_pred"] <- temp_data[J - 1, "EWM_pred"]

          } else if ((temp_data[J, "EWM_pred"] - temp_data[J - 1, "EWM_pred"]) < 0){

            temp_data[J, "EWM_pred"] <- temp_data[J - 1, "EWM_pred"]

          } else {

            temp_data[J, "EWM_pred"] <- temp_data[J, "EWM_pred"]

          }


        }

      }
    }
      #plot(y = temp_data$EWM, x = temp_data$DY, type = "p")
      #lines(y = temp_data$EWM_pred, x = temp_data$DY)

      errors <- NA
      final_parameters <- NA

      temp_data$first_diff <- NULL
      temp_data$method <- "brnn"

      list_temps[[b_holder]] <- temp_data
      b_holder = b_holder + 1

      if (fitted_save == TRUE){
        ggplot(temp_data, aes(x = DY, y = EWM_pred)) + geom_line() +
          geom_point(temp_data, mapping = aes(x = DY, y = EWM)) +
          ylab("EWM predicted") + theme_light()
        ggsave(paste0("brnn_", i, ".png"))
      }

    }
  } else if (current_fitting_method == "gam"){

    i = "QUPE_XYLEM_KRAS_2017_5"

    for (i in unique_keys){

      temp_data <- data_trees[data_trees$key == i, ]

      parameters <- data_site[data_site$key == i, ]

      gam_k <- parameters$gam_k
      gam_sp <- parameters$gam_sp

      # add zeros at the beggining
    if(add_zeros == TRUE){

      min_DY <- min(temp_data$DY) - 1
      row_list <- list()
      for (J in 1:min_DY){
        temp_row <- temp_data[1,]
        row_list[[J]] <- temp_row
      }
      new_rows <- do.call(rbind, row_list)
      new_rows$DY <- c(1:min_DY)
      new_rows$EWM <- 0

      temp_data <- rbind(new_rows, temp_data)
      temp_data$EWM <- as.numeric(temp_data$EWM)
    }

      output <- gam(EWM ~ s(DY, k = gam_k, bs ="ps", sp = gam_sp),
                    data = temp_data, method = "REML")

      temp_data$EWM_pred <- predict(output)

      # plot(y = temp_data$EWM, x = temp_data$DY, main = i)
      # lines(y = temp_data$EWM_pred, x = temp_data$DY, type = "l")

      if (post_process == TRUE){

      avg_fit <- mean(temp_data$EWM_pred)
      max_fit <- max(temp_data$EWM_pred)

      lagged_EWM_pred <- temp_data$EWM_pred[-length(temp_data$EWM_pred)]
      lagged_EWM_pred <- c(NA, lagged_EWM_pred)

      temp_data$first_diff <- temp_data$EWM_pred - lagged_EWM_pred
      temp_data[1, "first_diff"]  <- temp_data[2, "first_diff"]

      # In case of negative predictions
      if (sum(temp_data$EWM_pred < 0, na.rm = TRUE) > 0){

        shortcut1 <- temp_data
        shortcut1$neg_ind <- ifelse(shortcut1$EWM_pred < 0, TRUE, FALSE)

        th_DY <- as.numeric(max(shortcut1[shortcut1$neg_ind == TRUE, ][, "DY"]))
        temp_data$EWM_pred <- ifelse(temp_data$DY <= th_DY, 0,
                                     temp_data$EWM_pred)

        }


      temp_data$EWM_pred <- ifelse(temp_data$first_diff < 0 &
                                     (temp_data$EWM_pred < avg_fit), 0,
                                   ifelse(temp_data$EWM_pred < 0, 0,
                                          temp_data$EWM_pred))

      test <- temp_data[temp_data$EWM_pred > avg_fit, ]
      test <- test[test$first_diff < 0, ]

      if (sum(test$first_diff < 0, na.rm = TRUE) > 0){}

      for (J in 1:nrow(temp_data)){


        if (temp_data[J, "EWM_pred"] < 0.0001){

          next()

        } else if (temp_data[J, "first_diff"] < 0){

          temp_data[J, "EWM_pred"] <- temp_data[J - 1, "EWM_pred"]

        } else if ((temp_data[J, "EWM_pred"] - temp_data[J - 1, "EWM_pred"]) < 0){

          temp_data[J, "EWM_pred"] <- temp_data[J - 1, "EWM_pred"]

        } else {

          temp_data[J, "EWM_pred"] <- temp_data[J, "EWM_pred"]

          }


        }

     # plot(y = temp_data$EWM, x = temp_data$DY, main = i)
     # lines(y = temp_data$EWM_pred, x = temp_data$DY, type = "l")

      }

      errors <- NA
      final_parameters <- NA

      temp_data$method <- "gam"
      temp_data$first_diff <- NULL
      list_temps[[b_holder]] <- temp_data
      b_holder = b_holder + 1

      if (fitted_save == TRUE){
        ggplot(temp_data, aes(x = DY, y = EWM_pred)) +
          geom_line() + geom_point(temp_data, mapping = aes(x = DY, y = EWM)) +
          ylab("EWM predicted") + theme_light()
        ggsave(paste0("gam_", i, ".png"))
      }
    }

  } else {

      stop("Select proper fitting_method!")
    }

} # end of ut

  output_list <- list(fitted = do.call(rbind, list_temps),
                      gompertz_grid_search = final_parameters,
                      gompertz_grid_search_errors = errors_grid

                      )

  class(output_list) <- "wf"

  return(output_list)

}







