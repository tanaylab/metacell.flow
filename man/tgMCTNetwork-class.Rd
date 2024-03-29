% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mctnetwork.r
\docType{class}
\name{tgMCTNetwork-class}
\alias{tgMCTNetwork-class}
\alias{tgMCTNetwork}
\title{temporal netwrok over metacells}
\description{
Splitting metacells over a discrete time axis, defining manifold connections and estimated flows over them
}
\section{Slots}{

\describe{
\item{\code{mc_id}}{id of the metacell object we represent as a network}

\item{\code{mgraph_id}}{id of the mgraph object defining manifold structure}

\item{\code{times_nms}}{names of the time points (Default 1:T)}

\item{\code{mc_t}}{distribution of metacells (rows) over time points (cols)}

\item{\code{mc_mgraph}}{a data frame defining triplets mc1, mc2, distance.}

\item{\code{mc_manif_p}}{probability of moving between mcs - matrix}

\item{\code{mc_manif_cost}}{probability of moving between mcs}

\item{\code{network}}{- a data frame defining the network structure}
}}

