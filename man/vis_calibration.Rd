% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova.R
\name{vis_calibration}
\alias{vis_calibration}
\title{Viz calibration from simulation results}
\usage{
vis_calibration(
  pdat,
  metric = c("cslope", "clarge"),
  save = TRUE,
  lim = c(-100, 100),
  only_one = FALSE,
  theme_fun = theme(text = element_text(size = 13.5)),
  outdir = "."
)
}
\arguments{
\item{pdat}{data.frame; data for plotting}

\item{metric}{character; calibration slope or calibration in the large}

\item{save}{whether to save plots}

\item{lim}{numeric; y-axis limits}

\item{only_one}{logical; return only one plot}

\item{theme_fun}{custom ggplot2 theme to apply to each plot}

\item{outdir}{character; directory in which to save results and plots}
}
\value{
ggplot2 objects
}
\description{
Viz calibration from simulation results
}
