#' @method plot xpsg
#' @export

plot.xpsg <- function(x, ...){

  DOY <- NULL
  Width_pred <- NULL
  key <- NULL

p1 <- ggplot(x[[1]], aes(x = DOY, y = Width_pred, col = key)) +
        geom_line() + facet_grid(.~method)

return(p1)

}
