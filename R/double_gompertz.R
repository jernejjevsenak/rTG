#' double_gompertz
#'
#' A method to fit the double Gompertz function
#'
#' @param x independent variable
#' @param a1 parameter a1
#' @param a2 parameter a2
#' @param b1 parameter b1
#' @param b2 parameter b2
#' @param k1 parameter k1
#' @param k2 parameter k2
#'
#' @return a nls object with the function fit
#'
#' @keywords internal

# define the double_gompertz function
double_gompertz <- function(x, a1, b1, k1, a2, b2, k2) {
  y1 <- a1 * exp(-b1 * exp(-k1 * x))
  y2 <- a2 * exp(-b2 * exp(-k2 * x))
  return(y1 + y2)
}

