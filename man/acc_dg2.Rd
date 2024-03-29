% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PVBcorrect_functions.R
\name{acc_dg2}
\alias{acc_dg2}
\title{PVB correction by Begg and Greenes' method 2 (deGroot et al, one covariate)}
\usage{
acc_dg2(data, test, disease, covariate, description = TRUE)
}
\arguments{
\item{data}{A data frame, with at least "Test" and "Disease" variables.}

\item{test}{The "Test" variable name, i.e. the test result. The variable must be in binary; positive = 1, negative = 0 format.}

\item{disease}{The "Disease" variable name, i.e. the true disease status. The variable must be in binary; positive = 1, negative = 0 format.}

\item{covariate}{The name(s) of covariate(s), i.e. other variables associated with either test or disease status.
Specify as name vector, e.g. c("X1", "X2") for two or more variables. The variables must be in formats acceptable to GLM.}

\item{description}{Print the name of this analysis. The default is \code{TRUE}. This can be turned off for repeated analysis, for example in bootstrapped results.}
}
\value{
A data frame object containing the accuracy results.
}
\description{
Perform PVB correction by Begg and Greenes' method 2 as described in deGroot et al (2011),
  in which it also includes PPV and NPV calculation. This is limited to only one covariate.
}
\examples{
acc_dg2(data = cad_pvb, test = "T", disease = "D", covariate = "X3")
  # equivalent to acc_ebg(), saturated_model
}
\references{
\enumerate{
  \item{de Groot, J. A. H., Janssen, K. J. M., Zwinderman, A. H., Bossuyt, P. M. M., Reitsma, J. B., & Moons, K. G. M. (2011). Correcting for partial verification bias: a comparison of methods. Annals of Epidemiology, 21(2), 139–148.}
}
}
