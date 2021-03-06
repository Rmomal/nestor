% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Initialize.R
\name{boot_FitSparsePCA}
\alias{boot_FitSparsePCA}
\title{Finds initial cliques using a sparse PCA on bootstraps sub-samples}
\usage{
boot_FitSparsePCA(M, B, r, min.size = 1, cores = 1, unique = TRUE)
}
\arguments{
\item{M}{Gaussian proxy for the original dataset. Possibly obtained using PLNmodels.}

\item{B}{Number of bootstrap samples.}

\item{r}{Number of missing actors.}

\item{min.size}{Minimum number of neighbors of missing actors.}

\item{cores}{Number of cores for parallel computation (uses mclapply, not available for Windows).}

\item{unique}{Boolean for keeping only unique results.}
}
\value{
\itemize{
\item{cliqueList:}{ list of all cliques found. Each element is itself a list of size r.}
\item{nb_occ:}{ vector of the number of times each cliques has been found by sPCA.}}
}
\description{
Finds initial cliques using a sparse PCA on bootstraps sub-samples
}
\examples{
data=generate_missing_data(n=100,p=10,r=1,type="scale-free", plot=TRUE)
PLNfit=norm_PLN(data$Y)
boot_FitSparsePCA(PLNfit$MO, B=100, r=1)
}
