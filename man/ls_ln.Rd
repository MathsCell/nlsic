% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlsic.R
\name{ls_ln}
\alias{ls_ln}
\title{Linear Least Squares, least norm solution}
\usage{
ls_ln(a, b, rcond = 1e+10, mnorm = NULL, x0 = NULL)
}
\arguments{
\item{a}{matrix or its QR decomposition}

\item{b}{vector or matrix of right hand sides}

\item{rcond}{maximal condition number for rank definition}

\item{mnorm}{norm matrix (can be dense or sparse) for which \code{\%*\%} operation with a dense vector is defined}

\item{x0}{optional vector from which a least norm distance is searched for}
}
\value{
solution vector
}
\description{
Linear Least Squares, least norm solution
}
