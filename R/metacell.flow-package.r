#' metacell.flow
#'
#' @import dplyr
#' @import ggplot2
#' @import tgconfig
#' @import tgstat
#' @import tgutil
#' @import Matrix
#' @import lpsymphony
#' @importClassesFrom Matrix Matrix dgCMatrix dgeMatrix
#' @importFrom tidyr pivot_longer pivot_wider spread gather
#' @importFrom dbscan dbscan
#' @importFrom parallel mclapply
#' @importFrom zoo rollmean rollmedian
#' @importFrom pdist pdist
#' @importFrom cluster silhouette
#' @importFrom igraph cluster_louvain graph_from_data_frame
#' @importClassesFrom graph graphNEL
#' @importFrom graph plot addEdge addNode nodeRenderInfo
#' @name metacell.flow
#' @docType package
NULL


# this is important for running this package from scripts
#' @import methods
NULL
