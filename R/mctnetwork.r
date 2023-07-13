#' temporal netwrok over metacells
#'
#' Splitting metacells over a discrete time axis, defining manifold connections and estimated flows over them
#'
#' @slot mc_t distribution of metacells (rows) over time points (cols)
#' @slot mc_manif_cost  probability of moving between mcs
#' @slot mc_t_expansion_rate matrix metacells over time group. Expansion/reduction rate per metacell and time bin. Can included effects of proliferation, ingress/egress & cell death
#' @slot temporal_bin_time data.frame with columns bin and t giving for each temporal bin its actual time
#' @slot network - a data frame defining the network structure
#' @slot metacell_names - vector with metacell names
#' 
#'
#' @export tgMCTNetwork
#' @exportClass tgMCTNetwork
tgMCTNetwork <- setClass(
   "tgMCTNetwork",
	slots = c(
	  mc_t = "matrix",
	  mc_t_raw = "matrix",
	  mc_manif_cost = "data.frame",
	  mc_t_expansion_rate = 'matrix',
	  temporal_bin_time = 'data.frame',
	  metacell_names = 'vector',
	  network = "data.frame")
)

#' Construct a metacell time network
#'
#'
#' @param mc_t distribution of metacells (rows) over time points (cols)
#' @export
setMethod(
  "initialize",
  signature = "tgMCTNetwork",
  definition =
    function(.Object, mc_t,mc_t_raw,temporal_bin_time = NULL,metacell_names = NULL) {

		if(is.null(metacell_names)) {
			if(is.null(rownames(mc_t))) {
				warning("please provide metacell names")
			}
			metacell_names = rownames(mc_t)
		}
		.Object@metacell_names = metacell_names

		if(is.table(mc_t)) {
			mc_t = matrix(mc_t, ncol=ncol(mc_t),nrow =nrow(mc_t),dimnames = list(rownames(mc_t),colnames(mc_t)))
		}
		if(is.table(mc_t_raw)) {
			mc_t_raw = matrix(mc_t_raw, ncol=ncol(mc_t_raw),nrow =nrow(mc_t_raw),dimnames = list(rownames(mc_t_raw),colnames(mc_t_raw)))
		}


		mc_t = mc_t[metacell_names,]
		mc_t = t(t(mc_t)/colSums(mc_t))

		.Object@mc_t = as.matrix(mc_t)
		.Object@mc_t_raw = as.matrix(mc_t_raw[metacell_names,])

		if(is.null(temporal_bin_time)) {
			temporal_bin_time = data.frame(bin = c(1:ncol(mc_t)),t = as.double(c(1:ncol(mc_t))))
		}
		.Object@temporal_bin_time = temporal_bin_time

      return(.Object)
    }
)


#' Generate a new network in scdb
#'
#' This constructs a metacell time network object with mc_t entry and metacell_names and temporal_bin_time_entry 
#'
#' @param mc_t matrix distribution of metacells (rows) over time points (cols)
#' @param mc_t_raw matrix containing number of cells per metacell (rows) over time points (cols), not normalized
#' @param temporal_bin_time data.frame with columns bin and t giving for each temporal bin its actual time
#' @param metacell_names vector with metacell names corresponding to rownames of mc_t
#' 
#' @export
mc2_new_mctnetwork = function(mc_t,mc_t_raw,temporal_bin_time = NULL,metacell_names = NULL)
{
	mct = tgMCTNetwork(mc_t,mc_t_raw,temporal_bin_time = temporal_bin_time,metacell_names = metacell_names)
	return(mct)
}


#' Construct and initialize a mct time network from a graph and leak table
#'
#' This constructs a meta cell time network object from an mc and mgraph objects and assignments of cells to graphs
#'
#' @param net_id id of scdb network object ot be added
#' @param mgraph_id metacell manfold graph object
#' @param cell_time assigning of time (discrete) to cells
#' @param mc_proliferation_rate a named vector holding the proliferation rate parameter per metacell, with names being metacell names. All zero if you don't know what to expect. Proliferation rates should be given in real time units e.g. 4 corresponds to approx 4 divisions per time unit. Timbe bin times are specified in the temporal_bin_time parameter. 
#' @param temporal_bin_time a data.frame with columns bin and t containing the real times of each time bin - used to estimate growth rate per metacell and time point. If not specified, pseudotimes c(1:ncol(mc_t)) are used
#' @param t_exp the time period for the matrix expoenential used when computing manifold costs. This is default to 1, and should usually stay this way
#' @param T_cost the threshold of manifold edges costs that are removed from the network to make the resulting optimization problem smaller.
#' 
#' @return mct object with computed manifold costs and mctnetwork
#' @export

mctnetwork_from_mgraph_and_mc_proliferation_rate = function(mgraph,
								    					  mc_t,
														  mc_t_raw,
														  mc_proliferation_rate = NULL,
														  temporal_bin_time = NULL,
														  t_exp = 1,
														  T_cost = 1e+5,
														  metacell_names = NULL) {	
	

	#build an empty mct object
	mct = mc2_new_mctnetwork(mc_t,mc_t_raw,temporal_bin_time = temporal_bin_time,metacell_names = metacell_names)

	#initialize manifold costs
	mct = mctnetwork_comp_manifold_costs_from_mgraph(mct,
													 mgraph = mgraph,
													 t_exp=t_exp, 
													 T_cost=T_cost)

	
	if(is.null(mc_proliferation_rate)) {
		mc_proliferation_rate = rep(0,length(mct@metacell_names))
		names(mc_proliferation_rate) = metacell_names
	}
	# compute mc_t_expansion_rate from mc_proliferation_rate
	mc_t_expansion_rate = mctnetwork_calculate_relative_growth_rates_from_mc_expansion_rate(mct@mc_t,mc_expansion_rate = mc_proliferation_rate,temporal_bin_time = mct@temporal_bin_time)
	mct@mc_t_expansion_rate = mc_t_expansion_rate

	mct = mctnetwork_generate_network(mct = mct)

	return(mct)
}

#' Compute and update manifold costs for a mc time network
#'
#'
#' @param mct mct network object
#' @param mgraph data.frame with columns mc1, mc2 and dist
#' @param t_exp time parameter for the markov exponential
#' @param T_cost threshold on the mc-mc cost for being included
#'
#' @return an updated mct object with the entry mc_manif_cost updatetd
#' @export

mctnetwork_comp_manifold_costs_from_mgraph = function(mct,mgraph, t_exp = 1, T_cost = 1e+5)
{

  metacell_names = mct@metacell_names

  # mgraph might contain more metacells than the flow object (because flow is only computed on a restricted set of metacells)
  # remove all mgraph edges that are leaving the set of included metacell_names
  f = ( as.character(mgraph$mc1) %in% metacell_names ) & ( as.character(mgraph$mc2) %in% metacell_names )
  mgraph = mgraph[f,]
  mgraph$mc1 = as.character(mgraph$mc1)
  mgraph$mc2 = as.character(mgraph$mc2)

  df = mgraph[,c("mc2","mc1","dist")]
  colnames(df) = c("mc1","mc2","dist") 
  
  df_mc_index = data.frame(mc = metacell_names,mc_index = seq_along(metacell_names))

  mgraph2 = rbind(mgraph[,c("mc1","mc2","dist")], df)
  mgraph2 = unique(mgraph2)
  mgraph2$dist_inv = 1/mgraph2$dist

  mgraph2 = dplyr::left_join(mgraph2,dplyr::rename(df_mc_index,mc1 = mc,mc_ind1 = mc_index))
  mgraph2 = dplyr::left_join(mgraph2,dplyr::rename(df_mc_index,mc2 = mc,mc_ind2 = mc_index))
  
  adj_mat = sparseMatrix(i = mgraph2$mc_ind1,
  						 j = mgraph2$mc_ind2,
						 x = mgraph2$dist_inv,
						 dims = c(length(metacell_names),length(metacell_names)),
						 dimnames = list(metacell_names,metacell_names))

  #adj_mat = as.matrix(adj_mat)
  
  diag(adj_mat) = 0
  
  row_max = apply(adj_mat,1,max)
  median_rate = median(row_max)
  
  adj_mat = adj_mat/median_rate
  
  row_sums = rowSums(adj_mat)
  diag(adj_mat) = -row_sums
  
  trans_mat = expm(t_exp*adj_mat)
  
  #trans_mat = as.matrix(trans_mat)
  trans_mat = trans_mat/apply(trans_mat,1,max)

  #diag(trans_mat) = rowMaxs(trans_mat)

  cost_mat = round(10/trans_mat)
  cost_mat = as.matrix(cost_mat)
  
  manifold = as.data.frame.table(cost_mat)
  colnames(manifold) = c("mc1","mc2","cost")
  manifold = manifold[is.finite(manifold$cost),]
  manifold = manifold[manifold$cost < T_cost,]
  
  #colnames(manifold) = c("mc1","mc2","cost")

  #	mct@mc_cost_mat = cost_mat
	mct@mc_manif_cost = manifold
	return(mct)
}




mctnetwork_generate_network = function(mct)
{
	k_inf = 100000
  
	manifold = mct@mc_manif_cost
	mc_t_freq = mct@mc_t
	max_t = ncol(mc_t_freq)

	metacell_names = mct@metacell_names


	# calculate normalized growth rates per metacell and time point
	browser()
	growth_rates = mctnetwork_renormalize_mc_expansion_rates(mct@mc_t_expansion_rate,mc_t_freq = mct@mc_t)
	# add a column of 1's for the last time point
	growth_rates = cbind(growth_rates,rep(1,nrow(growth_rates)))
	
    
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
						mc1 = rep("src", n_mc),
						mc2 = metacell_names,
						time1 = rep(0, n_mc), time2 = rep(1, n_mc),
						type1 = rep("src", n_mc), type2 = rep("norm_b", n_mc),
						growth_rate = rep(1,n_mc))
	next_e_id = next_e_id + n_mc

	edges_capacity = list()
	edges_manifold = list()

	for(t in 1:max_t) {
		edges_capacity[[t]] = data.frame(
						from = rep(normal_nodes_t_back[[t]],1),
						to = rep(normal_nodes_t_front[[t]],1),
						ID = next_e_id:(next_e_id + n_mc - 1),
						capacity = ifelse(mc_t_freq[,t]>0, k_inf,0),
						min_capacity = rep(0, n_mc),
						cost = rep(0, n_mc),
						mc1 = metacell_names,
						mc2 = metacell_names,
						time1 = rep(t, n_mc), time2 = rep(t, n_mc),
						type1 = rep("norm_b", n_mc),
						type2 = rep("norm_f", n_mc),
						growth_rate = growth_rates[,t])

		next_e_id = next_e_id + n_mc

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
						mc1 = metacell_names,
						mc2 = rep("sink",n_mc),
						time1 = rep(max_t, n_mc), time2 = rep(max_t+1, n_mc),
						type1 = rep("norm_f", n_mc), type2 = rep("sink", n_mc),
						growth_rate = rep(1,n_mc))
	next_e_id = next_e_id + n_mc

	ec = do.call("rbind", edges_capacity)
	em = do.call("rbind", edges_manifold)

	network = rbind(edges_src, ec)
	network = rbind(network, em)
	network = rbind(network, edges_sink)

	# df_metacell_index_to_name_dictionary = data.frame(mc_index = as.character(c(-2,-1,c(1:nrow(mc_t_freq)))),
	# 									mc_name = c(c("sink","src"),mct@metacell_names))

	# network = dplyr::left_join(network,dplyr::rename(df_metacell_index_to_name_dictionary,mc1 = mc_index,mc1_name = mc_name))
	# network = dplyr::left_join(network,dplyr::rename(df_metacell_index_to_name_dictionary,mc2 = mc_index,mc2_name = mc_name))

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
	
	
	df_src_sink_nodes = data.frame(node = c(source_id,sink_id),t = c(0,max_t +1),mc = c(-1,-2),type = c("src","sink"))
	
	nodes = rbind(df_src_sink_nodes,df_norm_nodes)
	nodes = nodes[order(nodes$node),]
  rownames(nodes) = c(1:nrow(nodes))
  
  # ----------------------------------------------
  
	mct@network = network

	return(mct)
}




mctnetwork_calculate_relative_growth_rates_from_mc_expansion_rate = function(mc_t_freq,mc_expansion_rate,temporal_bin_time) {
  
  metacell_names = rownames(mc_t_freq)

  mc_expansion_rate = mc_expansion_rate[metacell_names]
  mean_rate = mean(mc_expansion_rate)
  
  t_diff = diff(temporal_bin_time$t)
  prolif_rate = exp(outer(mc_expansion_rate - mean_rate,t_diff)*log(2))
  
  mean_expansion_rate = colSums(mc_t_freq[,-ncol(mc_t_freq),drop = F] * prolif_rate) 
  
  prolif_rate = t(t(prolif_rate)/mean_expansion_rate)
  
  return(prolif_rate)
}

#' normalizes mc_t_expansion_rate given an mc_t distribution such that the total frequency after proliferation per time point still equals one given the frequencies mc_t
#' 
#' @param mc_t_expansion_rate matrix of metacells over time points minus the last time point (i.e. ncol(mc_t_expansion_rate) should equal ncol(mc_t) - 1)
#' @param mc_t matrix of metacell frequency over time points
#' 
#' @return renormalized mc_t_expansion_rate
#' @export
mctnetwork_renormalize_mc_expansion_rates = function(mc_t_expansion_rate,mc_t_freq) {

	common_rownames = intersect(rownames(mc_t_expansion_rate),rownames(mc_t_freq))
	if ( (length(common_rownames) < nrow(mc_t_expansion_rate) ) |  (length(common_rownames) < nrow(mc_t_freq) ) ) {
		stop("rownames of mc_t_expansion_rate and mc_t_freq are not the same")
	}
	
	mean_expansion_rate = colSums(mc_t_freq[common_rownames,-ncol(mc_t_freq),drop = F] * mc_t_expansion_rate[common_rownames,]) 

	mc_t_expansion_rate = t(t(mc_t_expansion_rate)/mean_expansion_rate)
	return(mc_t_expansion_rate)
}


#return lp contstrains in the std (e.g. symphony) style
mctnetwork_linear_constraints <- function(mct, total_flow, min_total_flow) 
{

  net = mct@network

  # add flow_scaling column. This will transform the capacity edge variables into relative quantities giving the amount of flow relativ to the expected mean flow.
  net$flow_scaling = rep(1,nrow(net))

  max_t = max(mct@network$time1)

  for (t in 1:max_t) {

    f_edges = (net$type1 == "norm_b") & (net$type2 == "norm_f") & (net$time1 == t)

    net$flow_scaling[f_edges] = mct@mc_t[as.character(net$mc1[f_edges]),t] * total_flow

  }

  # first construct the equality constraints at every node
  # influx*growth_rate  = outflux
  # constraints are indexed by the node IDs, leaving out the source and the sink for which no constraints are provided

  #no source constraint
  f = net$from != 1
  lhs_flow_constr_o = data.frame(id = (net$from[f] - 1),
									edge = net$ID[f],
									coef = -net$flow_scaling[f]) 

  #no sink consraint
  f = net$to != max(net$to)
  lhs_flow_constr_i = data.frame(id = net$to[f] - 1,
									edge = net$ID[f],
									coef=net$growth_rate[f]*net$flow_scaling[f])
  rhs_flow_constr = rep(0, max(net$to)- 2)
  dir_flow_constr = rep('==', max(net$to) - 2)

  next_cnstr_id = max(net$to) - 1
  
  # constraint on the total flow starting at the source
  src_adj = net[net$from == 1,"ID"]
  
  lhs_src_flow_constr = data.frame(id = rep(next_cnstr_id, length(src_adj)),
                                   edge = src_adj,
                                   coef=rep(1, length(src_adj)))
  rhs_src_flow_constr = total_flow
  dir_src_flow_constr = '=='
  next_cnstr_id = next_cnstr_id + 1
  
  # inequality constraints on the total flow per capacity layer per time point
  # should be >= min_total_flow
  
  lhs_tot_flow_constr = data.frame()
  for (t in 1:max(net$time1)) {
    
    f = ( net$time1 == t ) & (net$type1 == "norm_b")
    cap_adj = net$ID[f]

    lhs_tot_flow_constr_tmp = data.frame(id = rep(next_cnstr_id,length(cap_adj)),
                                         edge = cap_adj,
                                         coef = net$flow_scaling[f])
    
    lhs_tot_flow_constr = rbind(lhs_tot_flow_constr,lhs_tot_flow_constr_tmp)
    next_cnstr_id = next_cnstr_id + 1
  }
  
  rhs_tot_flow_constr = rep(min_total_flow,max(net$time1))
  dir_tot_flow_constr = rep('>=',max(net$time1))
 

  # Build constraints matrix
  lhs_df = rbind(lhs_flow_constr_o, lhs_flow_constr_i,lhs_src_flow_constr, lhs_tot_flow_constr)

  lhs = sparseMatrix(i = lhs_df$id,j =  lhs_df$edge, x= lhs_df$coef)
  #lhs = slam::simple_triplet_matrix(i = lhs_df$id, lhs_df$edge, v= lhs_df$coef)
  constraints <- list(
    lhs = lhs,
    dir = c(dir_flow_constr, dir_src_flow_constr, dir_tot_flow_constr),
    rhs = c(rhs_flow_constr, rhs_src_flow_constr, rhs_tot_flow_constr))

    mct@network = net

 return(list(mct = mct,constraints = constraints))
}

