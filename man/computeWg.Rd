% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Infer.R
\name{computeWg}
\alias{computeWg}
\title{Updates the variational edges weights inside the VEM.}
\usage{
computeWg(Rho, Omega, W, r, n, alpha, hist = FALSE)
}
\arguments{
\item{Rho}{Correlation matrix.}

\item{Omega}{Matrix filled with precision values common to all spanning trees.}

\item{W}{Edge weight matrix.}

\item{r}{Number of missing actors.}

\item{n}{Number of samples.}

\item{alpha}{Tempering parameter.}

\item{hist}{Prints the histogram of the log-values of OO and OH blocs if TRUE.}
}
\value{
The variational edge weights matrix.
}
\description{
Updates the variational edges weights inside the VEM.
}
