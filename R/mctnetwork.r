#' temporal netwrok over metacells
#'
#' Splitting metacells over a discrete time axis, defining manifold connections and estimated flows over them
#'
#' @slot mc_id id of the metacell object we represent as a network
#' @slot mgraph_id id of the mgraph object defining manifold structure
#' @slot times_nms names of the time points (Default 1:T)
#' @slot mc_t distribution of metacells (rows) over time points (cols)
#' @slot mc_mgraph a data frame defining triplets mc1, mc2, distance.
#' @slot mc_manif_p probability of moving between mcs - matrix
#' @slot mc_manif_cost  probability of moving between mcs
#' @slot network - a data frame defining the network structure
#'
#' @export tgMCTNetwork
#' @exportClass tgMCTNetwork
tgMCTNetwork <- setClass(
   "tgMCTNetwork",
	slots = c(
	  mc_id = "character",
	  mgraph_id = "character",
	  cell_time = "vector",
	  mc_t = "matrix",
	  mc_mgraph = "data.frame",
	  mc_cost_mat = "matrix",
	  mc_manif_cost = "data.frame",
	  network = "data.frame")
)

#' Construct a meta cell time network
#'
#'
#'
#' @param mc_id metacell object id 
#' @param mgraph_id metacell graph object id 
#' @param cell_time vector assigning a time (number) to each cell
#' @param mc_prolif factor per mc
#' @export

setMethod(
  "initialize",
  signature = "tgMCTNetwork",
  definition =
    function(.Object, mc_id, mgraph_id, cell_time) {
		.Object@mc_id = mc_id
		.Object@mgraph_id = mgraph_id
		.Object@cell_time = cell_time
		mc = scdb_mc(mc_id)
		if(is.null(mc)) {
			stop("MC-ERR unkown mc_id ", mc_id, " when building network")
		}
		mc_t = table(mc@mc, cell_time[names(mc@mc)])
		mc_t = matrix(mc_t, ncol=ncol(mc_t))
		mc_t = t(t(mc_t)/colSums(mc_t))

		.Object@mc_t = matrix(mc_t, ncol=ncol(mc_t))

		mgraph = scdb_mgraph(mgraph_id)
		if(is.null(mgraph)) {
			stop("MC-ERR unkown manifold id ", mgraph_id, " when building network")
		} 
		if(mgraph@mc_id != mc_id) {
			stop("MC-ERR mismatch of mc_id ", mc_id, " and mgraph ", mgraph_id, " when building network")
		}
		.Object@mc_mgraph = mgraph@mgraph

      return(.Object)
    }
)

#' Generate a new network in scdb
#'
#' This constructs a meta cell time network object without actually initializing it. 
#'
#' @param net_id id of scdb network object ot be added
#' @param mc_id metacell object
#' @param mgraph_id metacell manfold graph object
#' @param cell_time assigning of time (discrete) to cells
#' @export
mcell_new_mctnetwork = function(net_id, mc_id, mgraph_id, cell_time)
{
	scdb_add_mctnetwork(net_id, tgMCTNetwork(mc_id, mgraph_id, cell_time))
}

#' Construct and initialize a mct time network from a graph and leak table
#'
#' This constructs a meta cell time network object from an mc and mgraph objects and assignments of cells to graphs
#'
#' @param net_id id of scdb network object ot be added
#' @param mgraph_id metacell manfold graph object
#' @param cell_time assigning of time (discrete) to cells
#' @param mc_proliferation_rate_fn a table holding the proliferation rate parameter per metacell. All zero if you don't know what to expect. Proliferation rates should be given in real time units e.g. 4 corresponds to approx 4 divisions per time unit. Age group times are specified in the time_age_groups parameter. 
#' @param time_age_groups a vector containing the real times of each age group bin - used to estimate growth rate per metacell and time point. If not specified, pseudotimes c(1:n_age_groups) are used
#' @param t_exp the time period for the matrix expoenential used when computing manifold costs. This is default to 1, and should usually stay this way
#' @param T_cost the threshold of manifold edges costs that are removed from the network to make the resulting optimization problem smaller.
#' @param capacity_var_factor optionally a vector specifying how tight should the capacity constraint on each metacell is modeled. 
#' @param off_capacity_cost1 the slope of the piecewise linear capacity cost around 1. Typically a smaller number at an order of 1
#' @param off_capacity_cost2 the slope of the piecewise linear cpaacity cost away from 1 - tpyically a large number representing strong constraint for not meeting the capacity constraints at all (either too little or too much flow)
#' @export
mcell_mctnet_from_mgraph = function(net_id,
								mgraph_id,
								cell_time,
								mc_proliferation_rate_fn = NULL,
								time_age_groups = NULL,
								t_exp = 1,
								T_cost = 1e+5,
								capacity_var_factor = NULL,
								off_capacity_cost1 = 1,
								off_capacity_cost2 = 1000)
{
	
  mgraph = scdb_mgraph(mgraph_id)
	if(is.null(mgraph)) {
		stop("missing mgraph id ", mgraph_id, " when building mct net")
	}
	mc_id = mgraph@mc_id
	mc = scdb_mc(mc_id)
  
	
	if(!is.null(mc_proliferation_rate_fn)) {
	  if(!file.exists(mc_proliferation_rate_fn)) {
	    stop("cannot open proliferation rate table", mc_proliferation_rate_fn, " when constructing network")
	  }
	  df_prolif = read.table(mc_proliferation_rate_fn, h=T)
	  df_prolif = df_prolif[order(df_prolif$mc),]
	  proliferation_rate = df_prolif$proliferation_rate
	  names(proliferation_rate) = 1:length(proliferation_rate)
	} else {
	  proliferation_rate = rep(0,ncol(mc@e_gc))
	  names(proliferation_rate) = 1:length(proliferation_rate)
	}

#build an empty mct object
	mcell_new_mctnetwork(net_id = net_id,
                       mc_id = mc_id,
                       mgraph_id = mgraph_id,
                       cell_time = cell_time)

	mct = scdb_mctnetwork(net_id)

#initialize manifold costs
	mct = mctnetwork_comp_manifold_costs(mct, 
										t_exp=t_exp, 
										T_cost=T_cost)

	mct = mctnetwork_gen_network(mct = mct, 
	                             mc_proliferation_rate = proliferation_rate,time_age_groups = time_age_groups,
	                             capacity_var_factor = capacity_var_factor,
	                             off_capacity_cost1 = off_capacity_cost1,
	                             off_capacity_cost2 = off_capacity_cost2)

	scdb_add_mctnetwork(net_id, mct)
}

#' Compute and update manifold costs for a mc time network
#'
#'
#' @param mct mct network object
#' @param t_exp time parameter for the markov exponential
#' @param T_cost threshold on the cost for retan
#'
#' @return an updated mct object
#' @export

mctnetwork_comp_manifold_costs = function(mct, t_exp = 1, T_cost = 1000)
{
  mgraph = mct@mc_mgraph	
  df = mgraph[,c(2,1,3)]
  colnames(df) = colnames(mgraph)
  
  mgraph2 = rbind(mgraph, df)
  mgraph2 = unique(mgraph2[,c(1,2,3)])
  mgraph2$dist_inv = 1/mgraph2$dist
  
  adj_mat = sparseMatrix(i = mgraph2$mc1,j = mgraph2$mc2,x = mgraph2$dist_inv)
  adj_mat = as.matrix(adj_mat)
  
  diag(adj_mat) = 0
  
  row_max = apply(adj_mat,1,max)
  median_rate = median(row_max)
  
  adj_mat = adj_mat/median_rate
  
  row_sums = rowSums(adj_mat)
  diag(adj_mat) = -row_sums
  
  trans_mat = expm(t_exp*adj_mat)
  
  trans_mat = as.matrix(trans_mat)
  trans_mat = trans_mat/apply(trans_mat,1,max)

#  diag(trans_mat) = rowMaxs(trans_mat)

  cost_mat = round(10/trans_mat)
  cost_mat = as.matrix(cost_mat)
  rownames(cost_mat) = 1:nrow(cost_mat)
  colnames(cost_mat) = 1:ncol(cost_mat)
  
  manifold = as.data.frame.table(cost_mat)
  colnames(manifold) = c("mc1","mc2","cost")
  manifold = manifold[is.finite(manifold$cost),]
  manifold = manifold[manifold$cost < T_cost,]
  
  colnames(manifold) = c("mc1","mc2","cost")

	mct@mc_cost_mat = cost_mat
	mct@mc_manif_cost = manifold
	return(mct)
}

#' Compute and update network strcuture from mc, manifodld costs
#'
#' @param mct mct network object
#' @param mc_proliferation_rate proliferation rate per metacell per time unit. Positive rates correspond to expansion, negative rates to shrinkage/apoptosis. All zero if you don't know what to expect.
#' @param time_age_groups vector of real time for each age group. Real values are used for translating proliferation rates into expansion/shrinkage rates between subsequent age groups
#'
#' @return an updated mct object
#' @export

mctnetwork_gen_network = function(mct, mc_proliferation_rate,
                                  time_age_groups = NULL,
                                  capacity_var_factor = NULL,
                                  off_capacity_cost1 = 1,
                                  off_capacity_cost2 = 1000)
{
	k_inf = 100000
  

	
	manifold = mct@mc_manif_cost
	mc_t_freq = mct@mc_t
	max_t = ncol(mc_t_freq)
	if(is.null(time_age_groups)) {
	  time_age_groups = c(1:ncol(mct@mc_t))
	}

	# calculate normalized growth rates per metacell and time point
	growth_rates = mctnetwork_calculate_growth_rates(mct = mct,mc_proliferation_rate = mc_proliferation_rate,time_age_groups = time_age_groups)
	# add a column of 1's for the last time point
	growth_rates = cbind(growth_rates,rep(1,nrow(growth_rates)))
	

	if(is.null(capacity_var_factor)) {
		capacity_var_factor = rep(0.25, nrow(mc_t_freq))
	}

	edges = data.frame(from = c(), to = c(), 
						ID = c(), capacity = c(), cost = c(), 
						time1 = c(), time2 = c(), min_capacity = c(), 
						type1 = c(), type2 = c(),
						growth_rate = c())
  #types are: norm, src ,sink

	n_mc = nrow(mc_t_freq)

  #generate lists of nodes per network layer

	source_id = 1

	normal_nodes_t_back = list()
	normal_nodes_t_front = list()	

	next_n_id = 2
	for(t in 1:max_t) {
		normal_nodes_t_back[[t]] = next_n_id:(next_n_id + n_mc - 1)
		next_n_id = next_n_id + n_mc
		normal_nodes_t_front[[t]] = next_n_id:(next_n_id + n_mc - 1)
		next_n_id = next_n_id + n_mc
	}
	
  #last node id is the sink
	sink_id = next_n_id
	
  #add edges from source to back time 1
	next_e_id = 1
	edges_src = data.frame(
						from = rep(source_id, n_mc),
						to = normal_nodes_t_back[[1]],
						ID = next_e_id:(next_e_id + n_mc - 1),
						capacity = rep(k_inf, n_mc),
						min_capacity = rep(0,n_mc),
						cost = rep(0, n_mc),
						mc1 = rep(-1, n_mc),
						mc2 = 1:n_mc,
						time1 = rep(0, n_mc), time2 = rep(1, n_mc),
						type1 = rep("src", n_mc), type2 = rep("norm_b", n_mc),
						growth_rate = rep(1,n_mc))
	next_e_id = next_e_id + n_mc

	edges_capacity = list()
	edges_manifold = list()
	alpha = capacity_var_factor
	for(t in 1:max_t) {
		edges_capacity[[t]] = data.frame(
						from = rep(normal_nodes_t_back[[t]],4),
						to = rep(normal_nodes_t_front[[t]],4),
						ID = next_e_id:(next_e_id + 4*n_mc - 1),
						capacity = c((1/(1+alpha))*mc_t_freq[,t],
						             (1-1/(1+alpha))*mc_t_freq[,t],
										     (alpha)*mc_t_freq[,t],
										     ifelse(mc_t_freq[,t]>0, k_inf,0)),
						min_capacity = rep(0, 4*n_mc),
						cost = rep(c(-off_capacity_cost2, -off_capacity_cost1, 
										off_capacity_cost1, off_capacity_cost2), each=n_mc),
						mc1 = rep(1:n_mc, 4),
						mc2 = rep(1:n_mc, 4),
						time1 = rep(t, 4*n_mc), time2 = rep(t, 4*n_mc),
						type1 = rep("norm_b", 4*n_mc),
						type2 = rep("norm_f", 4*n_mc),
						growth_rate = rep(growth_rates[,t],4))

		next_e_id = next_e_id + 4*n_mc

		if(t != max_t) {
			n_e = nrow(manifold)
			norm_norm_from = normal_nodes_t_front[[t]][manifold$mc1]
			norm_norm_to = normal_nodes_t_back[[t+1]][manifold$mc2]

			edges_manifold[[t]] = data.frame(
						from = norm_norm_from,
						to = norm_norm_to,
						ID = next_e_id:(next_e_id+ n_e-1),
						capacity = rep(k_inf, n_e),
						min_capacity = rep(0,n_e),
						cost = manifold$cost,
						mc1 = manifold$mc1,
						mc2 = manifold$mc2,
						time1 = rep(t, n_e),
						time2 = rep(t+1, n_e),
						type1 = rep("norm_f", n_e),
						type2 = rep("norm_b", n_e),
						growth_rate = rep(1,n_e))

			next_e_id = next_e_id + n_e
		}	
	}
	
	edges_sink = data.frame(
						from = normal_nodes_t_front[[max_t]],
						to = rep(sink_id, n_mc),
						ID = next_e_id:(next_e_id + n_mc - 1),
						capacity = rep(k_inf, n_mc),
						min_capacity = rep(0, n_mc),
						cost = rep(0, n_mc),
						mc1 = 1:n_mc,
						mc2 = -1,
						time1 = rep(max_t, n_mc), time2 = rep(max_t+1, n_mc),
						type1 = rep("norm_f", n_mc), type2 = rep("sink", n_mc),
						growth_rate = rep(1,n_mc))
	next_e_id = next_e_id + n_mc

	ec = do.call("rbind", edges_capacity)
	em = do.call("rbind", edges_manifold)

	network = rbind(edges_src, ec)
	network = rbind(network, em)
	network = rbind(network, edges_sink)
	
	# construct dataframe with all node IDs
	# AT THE MOMENT WE DO NOT EXPORT THIS DATA.FRAME
	t_node = c()
	for (i in 1:max_t) {
	  t_node = c(t_node,rep(i,n_mc))
	}
	
	df_norm_nodes = data.frame(node = c(unlist(normal_nodes_t_back),unlist(normal_nodes_t_front)),
	                               t = rep(t_node,2),
	                               mc = rep(c(1:n_mc),2*max_t),
	                               type = c(rep("norm_b",n_mc*max_t),rep("norm_f",n_mc*max_t)))
	
	
	df_src_sink_nodes = data.frame(node = c(source_id,sink_id),t = c(0,max_t +1),mc = c(-1,-1),type = c("src","sink"))
	
	nodes = rbind(df_src_sink_nodes,df_norm_nodes)
	nodes = nodes[order(nodes$node),]
  rownames(nodes) = c(1:nrow(nodes))
  
  # ----------------------------------------------
  
	mct@network = network

	return(mct)
}

mctnetwork_build_growth_mats = function(mct, mc_proliferation_rate)
{
	max_mc = nrow(mct@mc_t)
	min_t = 1
	max_t = ncol(mct@mc_t)
	mc_t_freq = mct@mc_t

	growth_leak = matrix(0, nrow=max_mc, ncol=(max_t-min_t+1))
	for(t in min_t:max_t) {
		growth_leak[,t-(min_t-1)] = mc_t_freq[,t]*mc_proliferation_rate
	}

	tot_mass = colSums(mc_t_freq)

	lost_mass = colSums(growth_leak)/tot_mass

	growth_compensate = t(t(mc_t_freq)*c(0,lost_mass[-max_t]))
	return(list(leak=growth_leak, compensate=growth_compensate))
}

mctnetwork_calculate_growth_rates = function(mct,mc_proliferation_rate,time_age_groups) {
  
  mean_rate = mean(mc_proliferation_rate)
  
  t_diff = diff(time_age_groups)
  prolif_rate = exp(outer(mc_proliferation_rate - mean_rate,t_diff)*log(2))
  
  mean_expansion_rate = colSums(mct@mc_t[,-ncol(mct@mc_t),drop = F] * prolif_rate) 
  
  prolif_rate = t(t(prolif_rate)/mean_expansion_rate)
  
  return(prolif_rate)
}


#return lp contstrains in the std (e.g. symphony) style
mctnetwork_lp_constraints <- function(mct, total_flow, k_inf_cap,min_total_flow) 
{
  net = mct@network

  f_capped_edges = net$capacity < k_inf_cap
  n_cap = sum(f_capped_edges)
  lhs_cap_constr = data.frame(id=1:n_cap,
									edge= net$ID[f_capped_edges],
									coef=rep(1, n_cap))
  
  rhs_cap_constr = net$capacity[f_capped_edges]
  dir_cap_constr = rep('<=', times = n_cap)

  next_cnstr_id = nrow(lhs_cap_constr) + 1
  
  f_min_capped_edges = net$min_capacity > 0
  n_min_cap = sum(f_min_capped_edges)
  if(n_min_cap > 0) {
		lhs_min_cap_constr = data.frame(id = next_cnstr_id:(next_cnstr_id + n_min_cap - 1),
                                  edge = net$ID[f_min_capped_edges],
                                  coef = rep(1,n_min_cap))
  
		rhs_min_cap_constr = net$min_capacity[f_min_capped_edges]
		dir_min_cap_constr = rep('>=', times = n_min_cap)
  
		next_cnstr_id = next_cnstr_id + nrow(lhs_min_cap_constr)
		lhs_cap_constr = rbind(lhs_cap_constr, lhs_min_cap_constr)
		rhs_cap_constr = c(rhs_cap_constr, rhs_min_cap_constr)
		dir_cap_constr = c(dir_cap_constr, dir_min_cap_constr)
  }
  
  #no source constraint
  f = net$from != 1
  lhs_flow_constr_o = data.frame(id = (next_cnstr_id - 2 + net$from[f]),
									edge = net$ID[f],
									coef = rep(-1,sum(f))) 

  #no sink consraint
  f = net$to != max(net$to)
  lhs_flow_constr_i = data.frame(id =next_cnstr_id -2 + net$to[f],
									edge = net$ID[f],
									coef=net$growth_rate[f])
  rhs_flow_constr = rep(0, max(net$to)- 2)
  dir_flow_constr = rep('==', max(net$to) - 2)

  next_cnstr_id = next_cnstr_id + max(net$to) - 2
  
  # constraint on the total flow starting at the source
  src_adj = net[net$from == 1,"ID"]
  
  lhs_src_flow_constr = data.frame(id = rep(next_cnstr_id, length(src_adj)),
                                   edge = src_adj,
                                   coef=rep(1, length(src_adj)))
  rhs_src_flow_constr = total_flow
  dir_src_flow_constr = '=='
  next_cnstr_id = next_cnstr_id + 1

  # no constraint on total flow entering the sink
  if(0) {
    next_cnstr_id = next_cnstr_id + 1
    sink_adj = net[net$to == max(net$to),"ID"]
    lhs_sink_constr = data.frame(id = rep(next_cnstr_id, length(sink_adj)),
                                 edge = sink_adj,
                                 coef=rep(1, length(sink_adj)))
    rhs_sink_constr = total_flow
    dir_sink_constr = '=='
  }
  
  # inequality constraints on the total flow per capacity layer per time point
  # should be >= min_total_flow
  
  lhs_tot_flow_constr = data.frame()
  for (t in 1:max(net$time1)) {
    
    cap_adj = net[( net$time1 == t ) & (net$type1 == "norm_b"),"ID"]
    lhs_tot_flow_constr_tmp = data.frame(id = rep(next_cnstr_id,length(cap_adj)),
                                         edge = cap_adj,
                                         coef = rep(1,length(cap_adj)))
    
    lhs_tot_flow_constr = rbind(lhs_tot_flow_constr,lhs_tot_flow_constr_tmp)
    next_cnstr_id = next_cnstr_id + 1
  }
  
  rhs_tot_flow_constr = rep(min_total_flow,max(net$time1))
  dir_tot_flow_constr = rep('>=',max(net$time1))
 

  # Build constraints matrix
  lhs_df = rbind(rbind(lhs_cap_constr, lhs_flow_constr_o), lhs_flow_constr_i)
  lhs_df = rbind(lhs_df, lhs_src_flow_constr, lhs_tot_flow_constr)

#  lhs = sparseMatrix(i = lhs_df$id, lhs_df$edge, x= lhs_df$coef)
  lhs = slam::simple_triplet_matrix(i = lhs_df$id, lhs_df$edge, v= lhs_df$coef)
  constraints <- list(
    lhs = lhs,
    dir = c(dir_cap_constr, dir_flow_constr, dir_src_flow_constr, dir_tot_flow_constr),
    rhs = c(rhs_cap_constr, rhs_flow_constr, rhs_src_flow_constr, rhs_tot_flow_constr))

	return(constraints)
}

