% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova.R
\name{summarize_coef}
\alias{summarize_coef}
\title{Summarize coefficients from simulation results}
\usage{
summarize_coef(
  dat,
  svars = c("estimate", "bias"),
  funs = list(mean = mean, sd = sd, n = function(.x, ...) length(na.omit(.x)))
)
}
\arguments{
\item{dat}{data.frame; simulation results}

\item{svars}{character; which variables to summarize}

\item{funs}{list; list of functions to apply to \code{svars}}
}
\description{
Summarize coefficients from simulation results
}
