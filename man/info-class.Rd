% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TdAvgClass.R
\docType{class}
\name{info-class}
\alias{info-class}
\alias{info}
\title{An S4 class to represent information of peaks to keep, to order, and window, step size and coordinates of coverage value sampled.}
\description{
Store
order for plotting and filtering (integer ord.idx),
to use for return index, filter order of peaks (ord.idx) by whether to keep (kpt.idx)
width of up/down-stream window (integer w),
step size within the window (integer s),
coordinates of data collected with reference to window center (numeric points).
}
\section{Slots}{

\describe{
\item{\code{ord.idx}}{integer object represeting order of (subset of)ChIP-seq peaks to plot in the Tornado plot.}

\item{\code{w}}{integer object representing number of basepairs flanking the middle of intervals for plotting.}

\item{\code{s}}{integer object representing every n'th basepairs to sample data values.}

\item{\code{points}}{numeric object representing coordinates in basepairs with reference to the mid-point of interval for data collected.}
}}

