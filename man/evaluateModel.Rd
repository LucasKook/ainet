% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-utils.R
\name{evaluateModel}
\alias{evaluateModel}
\title{Evaluate model at new observations with \code{loss}}
\usage{
evaluateModel(m, newx, y_true, loss, ...)
}
\arguments{
\item{m}{model object}

\item{newx}{newdata}

\item{y_true}{observed outcome in newdata}

\item{loss}{function; loss to evaluate model with}

\item{...}{see methods}
}
\value{
evaluated loss at predictions and observed outcomes
}
\description{
Evaluate model at new observations with \code{loss}
}
