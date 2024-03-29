% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PVBcorrect_functions.R
\name{acc_em}
\alias{acc_em}
\title{PVB correction by EM-based logistic regression method}
\usage{
acc_em(
  data,
  test,
  disease,
  covariate = NULL,
  mnar = TRUE,
  description = TRUE,
  show_t = TRUE,
  t_max = 500,
  cutoff = 1e-04,
  ci = FALSE,
  ci_level = 0.95,
  ci_type = "basic",
  R = 999,
  seednum = NULL,
  t_print_freq = 100,
  return_t = FALSE,
  r_print_freq = 100
)
}
\arguments{
\item{data}{A data frame, with at least "Test" and "Disease" variables.}

\item{test}{The "Test" variable name, i.e. the test result. The variable must be in binary; positive = 1, negative = 0 format.}

\item{disease}{The "Disease" variable name, i.e. the true disease status. The variable must be in binary; positive = 1, negative = 0 format.}

\item{covariate}{The name(s) of covariate(s), i.e. other variables associated with either test or disease status.
Specify as name vector, e.g. c("X1", "X2") for two or more variables. The variables must be in formats acceptable to GLM.}

\item{mnar}{The default is assuming missing not at random (MNAR) missing data mechanism, \code{MNAR = TRUE}.
Set this to \code{FALSE} to obtain results assuming missing at random (MAR) missing data mechanism.
This will be equivalent to using \code{\link{acc_ebg}}.}

\item{description}{Print the name of this analysis. The default is \code{TRUE}. This can be turned off for repeated analysis, for example in bootstrapped results.}

\item{show_t}{Print the current EM iteration number t. The default is \code{TRUE}.}

\item{t_max}{The maximum iteration number for EM. Default \code{t_max = 500}. It is recommended to increase the
number when covariates are included.}

\item{cutoff}{The cutoff value for the minimum change between iteration.
This defines the convergence of the EM procedure. Default \code{cutoff = 0.0001}. This can be set to a larger value
to test the procedure.}

\item{ci}{View confidence interval (CI). The default is \code{FALSE}.}

\item{ci_level}{Set the CI width. The default is 0.95 i.e. 95\% CI.}

\item{ci_type}{Set confidence interval (CI) type. Acceptable types are "norm", "basic", "perc", and "bca",
for bootstrapped CI. See \code{\link[boot]{boot.ci}} for details.}

\item{R}{The number of bootstrap samples. Default \code{R = 999}.}

\item{seednum}{Set the seed number for the bootstrapped CI. The default is not set, so it depends on the user
to set it outside or inside the function.}

\item{t_print_freq}{Print the current EM iteration number t at each specified interval.
Default \code{t_print_freq = 100}.}

\item{return_t}{Return the final EM iteration number t.
This can be used for the purpose of checking the EM convergence. The default is \code{FALSE}, but is set to
TRUE when \code{ci = TRUE}.}

\item{r_print_freq}{Print the current bootstrap sample number at each specified interval.
Default \code{r_print_freq = 100}.}
}
\value{
A list object containing:
\describe{
  \item{boot_data}{An object of class "boot" from \code{\link[boot]{boot}}.
  Contains Sensitivity, Specificity, PPV, NPV and t (i.e. EM iteration taken for convergence).
  Use acc_em_object$boot_data$t[,5] to check the t.}
  \item{boot_ci_data}{A list of objects of type "bootci" from from \code{\link[boot]{boot.ci}}.
  Contains Sensitivity, Specificity, PPV, and NPV.}
  \item{acc_results}{The accuracy results.}
}
}
\description{
Perform PVB correction by EM-based logistic regression method.
}
\examples{
# For sample run, test with low R boot number, will take very long time to finish for large R

# Here we save to R objects for detailed analysis later
em_out = acc_em(data = cad_pvb, test = "T", disease = "D", ci = TRUE, seednum = 12345,
                R = 2, t_max = 1000, cutoff = 0.0005)
em_out$acc_results
em_out$boot_data$t  # bootstrapped data, 1:5 columns are Sn, Sp, PPV, NPV,
                    # t (i.e. EM iteration taken for convergence)
em_out$boot_ci_data

# With covariate, will take some time
# Also check the time taken to finish
start_time = proc.time()
em_outx = acc_em(data = cad_pvb, test = "T", disease = "D", covariate = "X3", ci = TRUE,
                 seednum = 12345, R = 2, t_max = 10000, cutoff = 0.0005)
                 # with covariate, better set larger t_max
elapsed_time = proc.time() - start_time; elapsed_time  # view elapsed time in seconds
em_outx$acc_results
em_outx$boot_data$t
}
\references{
\enumerate{
  \item{Kosinski, A. S., & Barnhart, H. X. (2003). Accounting for nonignorable verification bias in assessment of diagnostic tests. Biometrics, 59(1), 163–171.}
}
}
