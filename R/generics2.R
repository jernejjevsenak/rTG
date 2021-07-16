#' @method summary xpsg
#' @export

summary.xpsg <- function(object, ...){

  key <- NULL
  method <- NULL
  Width_pred <- NULL
  Width <- NULL

s1 <-  group_by(object[[1]], key, method) %>%
    summarise(RMSE =sqrt(sum((Width_pred - Width)^2)/n()))

return(s1)

}
