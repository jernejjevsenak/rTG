#' data_site
#'
#' A dataset with model fitting parameters for different regression methods.
#'
#' @format A data frame with 79 rows and 2 variables:
#' \describe{
#'   \item{Tissue}{XYLEM or PHLOEM}
#'   \item{Species}{Fagus sylvatica (FASY), Picea abies (PIAB), Quercus pubescens (QUPE)}
#'   \item{Site}{Panska reka (PAN), Karst (KRAS)}
#'   \item{Year}{2011, 2017}
#'   \item{gom_a}{The initial value for the Gompertz parameter a}
#'   \item{gom_b}{The initial value for the Gompertz parameter b}
#'   \item{gom_c}{The initial value for the Gompertz parameter c}
#'   \item{brnn_neurons}{The number of neurons for BRNN method}
#'   \item{gam_k}{The k parameter value for GAM method}
#'   \item{gam_sp}{The sp parameter value for GAM method}
#' }
#'
#' @export
"data_site"
