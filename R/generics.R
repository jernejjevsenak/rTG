#' @method plot xpsg
#' @export

plot.xpsg <- function(x, ...){

  doy <- NULL
  width_pred <- NULL
  key <- NULL

  p1 <- ggplot(x[[1]], aes(x = doy, y = width_pred, col = key)) +
        geom_line() + facet_grid(.~method)

  return(p1)

}
