#' Flow solution over a metacell temporal network
#'
#' Relating to a mctnetwork object, this store on particular flow solution
#'
#' @slot net_id - network model ID
#' @slot edge_flows - flow of mincost solution
#' @slot mc_t_infer - inferred fraction ofcells per mc and time
#' @slot mc_forward - forward flow distribution (mc_forward[[t]][mc1,mc2] = how much of the flow from mc1 is going to mc2 at time [t,t+1]
#' @slot mc_backward - backward flow distribution (mc_backward[[t]][mc1,mc2] = how much of the flow into mc1 came from  mc2 at time [t,t+1]
#'
#' @export tgMCTNetFlow
#' @exportClass tgMCTNetFlow
tgMCTNetFlow <- setClass(
   "tgMCTNetFlow",
	slots = c(
	  net_id = "character",
	  edge_flows = "vector",
	  mc_t_infer = "matrix",
	  mc_forward = "list",
	  mc_backward = "list")
)

#' Construct a flow object given an mctnetwork
#'
#'
#'
#' @param net_id metacell object id 
#' @param flow named vector flows on each net edge (name is edge ID)
#' @export

setMethod(
  "initialize",
  signature = "tgMCTNetFlow",
  definition =
    function(.Object, net_id, flows=NULL) {
		.Object@net_id = net_id
		if(!is.null(flows)) {
			.Object@edge_flows= flows
		}
#		.Object@mc_t_infer = NULL
#		.Object@mc_forward = NULL
#		.Object@mc_backward= NULL
		mcnet = scdb_mctnetwork(net_id)
		if(is.null(mcnet)) {
			stop("MC-ERR unkown mctnetwork id ", net_id, " when building a flow object")
		}
		if(!is.null(flows)) {
			if(mcnet@net$ID != names(flows)) {
				stop("initializing flows must be done with a flow vector over the associated network object edge IDs, current flows IDs are not consistent")
			}
		}
      return(.Object)
    }
)

#' Generate a new flow object in scdb
#'
#' This constructs a flow object and put in the flow scdb. Flow themselves are NULL by default, use gen_mincost to compute them.
#'
#' @param flow_id id of scdb network object ot be added
#' @param net_id id of scdb mc time network object 
#' @param init_mincost TRUE if you want to compute flows using the mincost algorith
#' @param flows flows per network edge (NULL by default, if not NULL it will supress the init_mincost flag)
#' @export
mcell_new_mctnetflow= function(flow_id, net_id, 
							init_mincost = F, flow_tolerance=0.01,
							max_flow_tolerance = 0.05,
							flows=NULL)
{
	mcf = tgMCTNetFlow(net_id, flows = flows)
	if(is.null(flows) & init_mincost) {
		mcf = mctnetflow_gen_mincost(mcf, flow_tolerance,max_flow_tolerance)
	}
	scdb_add_mctnetflow(flow_id, mcf)
}

#' Compute mincost flow and update the mct objects with flows and 
#' inferred mc capacity per time 
#'
#' @param mcf flow object
#' @param flow_tolerance how much flow we should miss per time
#'
#' @return an updated flow object
#' @export
mctnetflow_gen_mincost = function(mcf, flow_tolerance=0.01,max_flow_tolerance = 0.05)
{
	if (max_flow_tolerance < flow_tolerance | max_flow_tolerance >= 1) {
	  stop("max_flow_tolerance must be larger than flow_tolerance (and smaller than 1)")
	}
  mctnet = scdb_mctnetwork(mcf@net_id)
	net = mctnet@network
	ncnstr = mctnetwork_lp_constraints(mctnet, 1-flow_tolerance, 100,1 - max_flow_tolerance)

	edge_costs = net$cost[order(net$ID)]
	sol = lpsymphony::lpsymphony_solve_LP(obj = edge_costs,
										mat = ncnstr$lhs, 
										rhs = ncnstr$rhs, 
										dir = ncnstr$dir)
  
	flows = sol$solution
	mcf@edge_flows = flows[net$ID]
	names(mcf@edge_flows) = net$ID	

	net$flow = flows[net$ID]
	f_cap = net$mc1 == net$mc2 & net$time1 == net$time2
	net_cap = net[f_cap,]

	mc_time_weight = summarize(group_by(net_cap,mc1,time1),tot_weight = sum(flow))
	mc_t_post = sparseMatrix(i = as.numeric(mc_time_weight$mc1),
									j = as.numeric(mc_time_weight$time1), 
									x = mc_time_weight$tot_weight)
	mc_t_post = as.matrix(mc_t_post)
	rownames(mc_t_post) = c(1:nrow(mc_t_post))
	colnames(mc_t_post) = c(1:ncol(mc_t_post))

	mcf@mc_t_infer = mc_t_post
	return(mcf)
}


#' Compute matrices for forward and backwrd propagation using the flows
#'
#' @param mcf mcf flow object
#'
#' @return an updated mcf object
#' @export

mctnetflow_comp_propagation= function(mcf)
{
	mctnet = scdb_mctnetwork(mcf@net_id)
#flows mc,mc,t
	max_t = ncol(mctnet@mc_t)
	for(t in 1:(max_t-1)) {
		mc_flows = mctnetflow_get_flow_mat(mcf, mctnet, t)

		mc_forward = mc_flows/rowSums(mc_flows)
		mc_backward= t(t(mc_flows)/colSums(mc_flows))
		mc_backward[is.nan(mc_backward)] = 0
		mc_forward[is.nan(mc_forward)] = 0
		mcf@mc_forward[[t]] = mc_forward
		mcf@mc_backward[[t]] = mc_backward
	}
	return(mcf)
}

   
#' Compute matrix of flows between MCs in a given time point
#'
#' @param mcf flow object
#' @param mctnet mctnet network object
#' @param time flows will be computed for the (time,time+1) interval. Time=-1 will generate total flow over all times
#'
#' @return a matrix on metacells
#' @export

mctnetflow_get_flow_mat = function(mcf, mctnet, time)
{
	net = mctnet@network
	if(time == -1) {
		f_t = net$type1 == "norm_f" & net$type2 == "norm_b"
	} else {
		f_t = net$time1 == time & net$time2==time+1

	}

	net$flow  = mcf@edge_flows
	net_t = net[f_t,] 
   flow = as.data.frame(summarize(group_by(net_t, mc1, mc2),
													tot_flow = sum(flow)))
    
   mc_mat = pivot_wider(data = flow, 
				names_from = mc2, 
				values_from = tot_flow,
				values_fill = list(tot_flow = 0))

   mc_mat = as.data.frame(mc_mat)
   rownames(mc_mat) = mc_mat$mc1
   mc_mat = mc_mat[,-1]
	max_mc = ncol(mc_mat)
	if(time == -2) {
		mc_mat = mc_mat[as.character(c(-1, 1:max_mc)), as.character(1:max_mc)]
		mc_mat = cbind(rep(0,nrow(mc_mat)), mc_mat)
	} else {
		mc_mat = mc_mat[as.character(1:max_mc), as.character(1:max_mc)]
	}
   mc_mat = as.matrix(mc_mat)

	return(mc_mat)
}

#' Compute backward and forward flow propgation of metacell probability from time t
#'
#' @param mct mct network object
#' @param t flows will be computed for the (time,time+1) interval
#' @param mc_p probabilities at time t
#'
#' @return a list with two elements: probs is a matrix of probabilities over metacells (rows)  and time (columns). step_m is a list of sparse matrices inferred flows between metacells per time.
#' @export

mctnetflow_propogate_from_t = function(mcf, t, mc_p)
{
	mctnet = scdb_mctnetwork(mcf@net_id)
	max_t = ncol(mctnet@mc_t)
	step_m = list()
	probs = matrix(0, nrow = nrow(mctnet@mc_t), ncol=ncol(mctnet@mc_t))
	probs[,t] = mc_p
	if(t > 1) {
		for(i in (t-1):1) {
			step_m[[i]] = Matrix(t(t(as.matrix(mcf@mc_backward[[i]])) * probs[,i+1]), sparse=T)
			
			probs[,i] = as.matrix(mcf@mc_backward[[i]]) %*% probs[,i+1]
		}
	}
	if(t < max_t) {
		for(i in (t+1):max_t) {
			step_m[[i-1]] = Matrix(as.matrix(mcf@mc_forward[[i-1]]) * probs[,i-1], sparse=T)
			probs[,i] = t(probs[,i-1]) %*% as.matrix(mcf@mc_forward[[i-1]])
		}
	}
	return(list(probs=probs, step_m=step_m))
}

#' Given a list of genes and a matrix, this compute mean umi per metacell (rows) and time (column).
#'
#' @param mcf_id  flow object id
#' @param mat_id  umi matrix object
#' @param genes list of gene names 
#' @param min_percentile percentile value to use as minimum threshold (this is used to avoid very small values and regularize color scale, default 0.05)
#'
#' @return A matrix on metacells and time with mean umi fraction over the gene module
#' @export
#'
mctnetflow_gen_gmod_mc_t = function(mcf_id, mat_id,  genes, min_percentile=0.05)
{
	mcf = scdb_mctnetflow(mcf_id)
	if(is.null(mcf)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	mctnet = scdb_mctnetwork(mcf@mct_id)
	mc = scdb_mc(mctnet@mc_id)
	mat = scdb_mat(mat_id)
	genes = intersect(genes, rownames(mat@mat))
	gmod_tot = colSums(mat@mat[genes,names(mc@mc)])/colSums(mat@mat[,names(mc@mc)])
	gmod_df = data.frame(mc = mc@mc, 
								t = mctnet@cell_time[names(mc@mc)], 
								gmod = gmod_tot)
	mc_time_gmod = summarise(group_by(gmod_df, mc, t), gmod=mean(gmod))
	mc_t_gmod_m = sparseMatrix(i=mc_time_gmod$mc, 
									j=mc_time_gmod$t, 
									x=mc_time_gmod$gmod)
	mc_t_gmod_m = pmax(mc_t_gmod_m, quantile(mc_time_gmod$gmod, min_percentile))
	return(mc_t_gmod_m)
}

#' This generate two large matrices showing expression of genes per metacell ad time point, as well as the expression of the gene in inferred ancestral states given the flow model
#'
#' @param mcf_id  mcf network object
#' @param mat_id  umi matrix object
#' @param genes list of gene names 
#'
#' @return A matrix on metacells and time with mean umi fraction over the gene module
#' @export
#'
mctnetflow_gen_e_gmt_p = function(mcf_id, mat_id,  genes)
{
	mcf = scdb_mctnetflow(mcf_id)
	if(is.null(mcf)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	mctnet = scdb_mctnetwork(mcf@mct_id)
	if(is.null(mctnet)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	mc = scdb_mc(mctnet@mc_id)
	mat = scdb_mat(mat_id)
	genes = intersect(genes, rownames(mat@mat))

	max_t = length(mcf@mc_backward) + 1

	csize = colSums(mat@mat)
	e_gm_t = list()
	tot_m_t = list()
	for(t in 1:max_t) {
		cell_t = intersect(names(mc@mc), names(which(mctnet@cell_time==t)))
		message("at ", t, " with ", length(cell_t), " cells")
		mgt = as.matrix(mat@mat[genes, cell_t])
#		tot_gm = t(tgs_matrix_tapply(mgt, mc@mc[cell_t], sum))
		tot_gm = t(apply(mgt,1, function(x) tapply(x,mc@mc[cell_t], sum)))
		tot_c = tapply(csize[cell_t], mc@mc[cell_t], sum)
		tot_m_t[[t]] = tot_c
		e_gm_t[[t]] = t(t(tot_gm)/as.vector(tot_c))
		rownames(e_gm_t[[t]]) = genes
	}
	e_gm_t_prev = list()
	for(t in 2:max_t) {
		t_mcs = colnames(e_gm_t[[t-1]])
		t_mcs2= colnames(e_gm_t[[t]])
		e_gm_t_prev[[t]] = Matrix(e_gm_t[[t-1]] %*% as.matrix(mcf@mc_backward[[t-1]][t_mcs, t_mcs2]), sparse=T)
	}
	
	e_gm_t_prev2 = list()
	for(t in 3:max_t) {
	  t_mcs0 = colnames(e_gm_t[[t-2]])
	  t_mcs = colnames(e_gm_t[[t-1]])
	  t_mcs2= colnames(e_gm_t[[t]])
	  e_gm_t_prev2[[t]] = Matrix(e_gm_t[[t-2]] %*% as.matrix(mcf@mc_backward[[t-2]][t_mcs0, t_mcs]) %*% as.matrix(mcf@mc_backward[[t-1]][t_mcs, t_mcs2]), sparse=T)
	}
	
	return(list(e_gm_t = e_gm_t, e_gm_t_prev = e_gm_t_prev, e_gm_t_prev2 = e_gm_t_prev2, tot_m_t = tot_m_t))
}

#' Compute matrix of flows over cell types
#'
#' @param mcf mcf network object
#' @param min_time minimum time point
#' @param max_time maximum time point
#'
#' @return a list of matrices show flows from type t to t'
#' @export

mctnetflow_get_type_flows = function(mcf, time, max_time)
{
	mctnet = scdb_mctnetwork(mcf@mct_id)
	net = mctnet@network
	if(is.null(mcf@edge_flows)) {
		stop("trying to query uninitialized flows")
	}
	net$flow = mcf@edge_flows
	mc = scdb_mc(mctnet@mc_id)

	all_types = unique(mc@colors)
	mct_mats = list()
	for(t in time:(max_time-1)) {
		f_t = net$time1 == t & net$time2 == t+1 &
					net$type1 != "growth" & net$type2!="growth"

		net_t = net[f_t,] 
		net_t$mc_t1 = mc@colors[as.numeric(net_t$mc1)]
		net_t$mc_t2 = mc@colors[as.numeric(net_t$mc2)]
		flow = as.data.frame(summarize(group_by(net_t, mc_t1, mc_t2),
													tot_flow = sum(flow)))
    
	   mct_mat = pivot_wider(data = flow, 
				names_from = mc_t2, 
				values_from = tot_flow,
				values_fill = list(tot_flow = 0))

	   mct_mat = as.data.frame(mct_mat)
		rownames(mct_mat) = mct_mat$mc_t1
	   mct_mat = mct_mat[,-1]
		max_mc = ncol(mct_mat)
		mct_mat = mct_mat[all_types, all_types]
		mct_mat = as.matrix(mct_mat)
		mct_mats[[t]] = mct_mat
	}

	return(mct_mats)
}

mctnetflow_get_egc_on_cluster_transition = function(mcf, min_time, max_time, type1, type2, mc_type=NULL)
{
	mctnet = scdb_mctnetwork(mcf@mct_id)
	mc = scdb_mc(mctnet@mc_id)
	e_gc = mc@e_gc
	net = mctnet@network
	net$flow = mcf@edge_flows

	if(is.null(mc_type)) {
		mc_type = mc@colors
		names(mc_type) = as.character(1:length(mc_type))
	}

#	flow_mm = mctnetwork_get_flow_mat(mct, time, max_time=time)

	f_t = net$time1 >= min_time & net$time2 <= max_time &
					net$time1 == net$time2-1 &
					net$type1 != "growth" & net$type2!="growth"

	net = net[f_t,]
	f_types = mc_type[as.numeric(net$mc1)]==type1 & mc_type[as.numeric(net$mc2)]==type2

	net = net[f_types,]

	src_mc_wgt = tapply(net$flow, net$mc1, sum)	
	targ_mc_wgt = tapply(net$flow, net$mc2, sum)
	src_mc_wgt_n = as.vector(src_mc_wgt/sum(src_mc_wgt))
	names(src_mc_wgt_n) = names(src_mc_wgt)
	targ_mc_wgt_n = as.vector(targ_mc_wgt/sum(targ_mc_wgt))
	names(targ_mc_wgt_n) = names(targ_mc_wgt)

	src_e_gc = colSums(t(e_gc[,names(src_mc_wgt_n)]) * src_mc_wgt_n)
	targ_e_gc = colSums(t(e_gc[,names(targ_mc_wgt_n)]) * targ_mc_wgt_n)

	return(data.frame(src = src_e_gc, targ = targ_e_gc, lf = log2(1e-5+targ_e_gc)-log2(1e-5+src_e_gc)))
}
