% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameters.R
\docType{data}
\name{parameters}
\alias{parameters}
\title{parameters}
\format{
A data frame with 79 rows and 2 variables:
\describe{
  \item{Tissue}{XYLEM or PHLOEM}
  \item{Species}{Fagus sylvatica (FASY), Picea abies (PIAB), Quercus pubescens (QUPE)}
  \item{Site}{Panska reka (PAN), Karst (KRAS)}
  \item{Year}{2011, 2017}
  \item{Tree}{Tree ID indicators ranging from 1 to 6}
  \item{gom_a}{The initial value for the Gompertz parameter a}
  \item{gom_b}{The initial value for the Gompertz parameter b}
  \item{gom_k}{The initial value for the Gompertz parameter c}
  \item{brnn_neurons}{The number of neurons for BRNN method}
  \item{gam_k}{The k parameter value for GAM method}
  \item{gam_sp}{The sp parameter value for GAM method}
}
}
\usage{
parameters
}
\description{
data frame with model fitting parameters for different regression methods.
}
\keyword{datasets}
