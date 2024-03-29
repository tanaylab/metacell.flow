% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mctnetflow.r
\name{mctnetflow_get_type_flows}
\alias{mctnetflow_get_type_flows}
\title{Compute matrix of flows over cell types}
\usage{
mctnetflow_get_type_flows(mcf, time, max_time)
}
\arguments{
\item{mcf}{mcf network object}

\item{max_time}{maximum time point}

\item{min_time}{minimum time point}
}
\value{
a list of matrices show flows from type t to t'
}
\description{
Compute matrix of flows over cell types
}
