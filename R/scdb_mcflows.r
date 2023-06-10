#' Add mc2tnetwork to single-cell database
#'
#' @param net_id - id of mc2tnetwork
#' @param mct - mc2tnetwork object
#' @param scdb_dir - directory of single-cell database
#' 
#' @export 
#'  
scdb_add_mc2tnetwork = function(mct,net_id,scdb_dir) {

    fn = sprintf("%s/mc2tnetwork.%s.rds",scdb_dir,net_id)
    saveRDS(object = mct,file = fn)
}



#' get mc2tnetwork object
#'
#' @param net_id - id of mc2tnetwork
#' @param scdb_dir - directory of single-cell database
#' 
#' @export 
#'  
scdb_mc2tnetwork = function(net_id,scdb_dir) {

    fn = sprintf("%s/mc2tnetwork.%s.rds",scdb_dir,net_id)
    mct = readRDS(file = fn)
	return(mct)
}


#' Add mc2tnetflow to single-cell database
#'
#' @param flow_id - id of mc2tnetflow
#' @param mcf - mc2tnetflow object
#' @param scdb_dir - directory of single-cell database
#' 
#' @export 
#'  
scdb_add_mc2tnetflow = function(mcf,flow_id,scdb_dir) {

    fn = sprintf("%s/mc2tnetflow.%s.rds",scdb_dir,flow_id)
    saveRDS(object = mcf,file = fn)
}



#' get mc2tnetflow object
#'
#' @param flow_id - id of mc2tnetflow
#' @param scdb_dir - directory of single-cell database
#' 
#' @export 
#'  
scdb_mc2tnetflow = function(flow_id,scdb_dir) {

    fn = sprintf("%s/mc2tnetflow.%s.rds",scdb_dir,flow_id)
    mcf = readRDS(file = fn)
	return(mcf)
}



#' Add mgraph to single-cell database
#'
#' @param mgraph_id - id of mgraph
#' @param mgraph - mgraph object
#' @param scdb_dir - directory of single-cell database
#' 
#' @export 
#'  
scdb_add_mgraph2 = function(mgraph,mgraph_id,scdb_dir) {

    fn = sprintf("%s/mgraph2.%s.tsv",scdb_dir,mgraph_id)
    write.table(x = mgraph,file = fn,sep ='\t')

}


#' Get mgraph from single-cell database
#'
#' @param mgraph_id - id of mgraph
#' @param scdb_dir - directory of single-cell database
#' 
#' @export 
#'  
scdb_mgraph2 = function(mgraph_id,scdb_dir) {

    fn = sprintf("%s/mgraph2.%s.tsv",scdb_dir,mgraph_id)
    mgraph = read.table(file = fn,sep ='\t',h = T)
	return(mgraph)
}
