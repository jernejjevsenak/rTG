#' @method summary xpsg
#' @export

summary.xpsg <- function(object, ...){

  key <- NULL
  method <- NULL
  width_pred <- NULL
  width <- NULL

  s1 <-  group_by(object[[1]], key, method) %>%
    summarise(RMSE =sqrt(sum((width_pred - width)^2)/n()))

  return(s1)

}
