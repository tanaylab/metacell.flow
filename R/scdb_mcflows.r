#' Initializing scdb_flow
#'
#' This initialize scdb_flow to an already initialed scdb. This will extend the db with mctnetwork objects
#' 
#'
#' @examples
#' \dontrun{
#' # we first initialize a db object
#' scdb_flow_init()
#  # we can see which objects are available
#' net = scdb_mctnetwork("mynet")
#' }
#'
#' @export
scdb_flow_init = function()
{
	if(!exists(".scdb")) {
		stop("first initialized scdb using scdb_init")
	}
	base_dir = .scdb_base
	.scdb[["mctnetwork"]] <<- list()
	.scdb[["mctnetflow"]] <<- list()
}

#' scdb_mctnetwork - get a mctnetwork object
#'
#' @param id - id of mctnetwork
#'
#' @export
#'
scdb_mctnetwork = function(id) 
{
	return(.scdb_get_obj(id, "mctnetwork"))
}

#' scdb_add_mctnetwork - add mctnetwork to the DB and cahce
#'
#' @param id - id of mctnetwork
#' @param mctnetwork - mctnetwork data frame
#'
#' @export
#'
scdb_add_mctnetwork = function(id, mctnetwork) 
{
	if(class(mctnetwork)[1] != "tgMCTNetwork") {
		stop("Cannot add non tgMCTNetwork object as a mctnetwork in scdb")
	}
	.scdb_add_obj(id, "mctnetwork", mctnetwork);
}

#' scdb_del_mctnetwork - del mctnetwork from the DB and cahce
#'
#' @param id - id of mctnetwork
#'
#' @export
scdb_del_mctnetwork = function(id)
{
	.scdb_del_obj(id, "mctnetwork");
}

#' scdb_mctnetflow - get a mctnetflow object
#'
#' @param id - id of mctnetflow
#'
#' @export
#'
scdb_mctnetflow = function(id) 
{
	return(.scdb_get_obj(id, "mctnetflow"))
}

#' scdb_add_mctnetflow - add mctnetflow to the DB and cahce
#'
#' @param id - id of mctnetflow
#' @param mctnetflow - mctnetflow data frame
#'
#' @export
#'
scdb_add_mctnetflow = function(id, mctnetflow) 
{
	if(class(mctnetflow)[1] != "tgMCTNetFlow") {
		stop("Cannot add non tgMCTNetwork object as a mctnetflow in scdb")
	}
	.scdb_add_obj(id, "mctnetflow", mctnetflow);
}

#' scdb_del_mctnetflow - del mctnetflow from the DB and cahce
#'
#' @param id - id of mctnetflow
#'
#' @export
scdb_del_mctnetflow = function(id)
{
	.scdb_del_obj(id, "mctnetflow");
}
