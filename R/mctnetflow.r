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
	  mct = "tgMCTNetwork",
	  edge_flows = "vector",
	  mc_t_infer = "matrix",
	  mc_forward = "list",
	  mc_backward = "list")
)

#' Construct a flow object given an mctnetwork
#'
#'
#'
#' @param mct mc2tnetwork object 
#' @param flow named vector flows on each net edge (name is edge ID)
#' @export

setMethod(
  "initialize",
  signature = "tgMCTNetFlow",
  definition =
    function(.Object, mct, flows=NULL) {
		.Object@mct = mct
		if(!is.null(flows)) {
			.Object@edge_flows= flows
		}
#		.Object@mc_t_infer = NULL
#		.Object@mc_forward = NULL
#		.Object@mc_backward= NULL
		if(!is.null(flows)) {
			if(mct@network$ID != names(flows)) {
				stop("initializing flows must be done with a flow vector over the associated network object edge IDs, current flows IDs are not consistent")
			}
		}
      return(.Object)
    }
)


generate_piecewise_linear_mosek_problem = function(mct,
                                                   mc_t_sd = NULL,
                                                   total_flow = 0.99,
                                                   min_total_flow = 0.95,
                                                   lambda = 1,
                                                   mu = 1,
                                                   sd_over_mean_upper_threshold = 4,
                                                   sd_over_mean_lower_threshold = 0.0125,
                                                   time_bin_capacity_cost_to_zero = NULL) {

    k_inf = 1e+5
    max_cost = 1e+5
    n_edges = nrow(mct@network)

	mc_t_mean = mct@mc_t



    # linear constraints of the network

    cmp_constraints = mctnetwork_linear_constraints(mct = mct,
                                                       total_flow = total_flow,
                                                       min_total_flow = min_total_flow)

    linear_constraints = cmp_constraints$constraints
    mct  = cmp_constraints$mct
    # additionally we add below positivity constraints on all the variables and equal-zero-constraints on all capacity edges for which mc_t is zero

    max_t = max(mct@network$time1)
	
    net = mct@network

	df_metacell_index_to_name_dictionary = data.frame(mc_index = c(-2,-1,c(1:length(mct@metacell_names))),
										mc_name = c(c("sink","src"),mct@metacell_names))

	net = dplyr::left_join(net,dplyr::rename(df_metacell_index_to_name_dictionary,mc1_ind = mc_index,mc1 = mc_name))
	net = dplyr::left_join(net,dplyr::rename(df_metacell_index_to_name_dictionary,mc2_ind = mc_index,mc2 = mc_name))


    net_f = net[net$type1 == "norm_b" & net$type2 == "norm_f",]

    capacity_edge_dictionary = as.matrix(sparseMatrix(i = as.integer(net_f$mc1_ind),j = as.integer(net_f$time1),x = net_f$ID),dims = c(length(mct@metacell_names),max_t))
	rownames(capacity_edge_dictionary) = mct@metacell_names
    
    # add additional constraint

    if(!is.null(time_bin_capacity_cost_to_zero)) {
        f_pos = mct@mc_t > 0 | col(mct@mc_t) %in% time_bin_capacity_cost_to_zero
        capacity_edges_time_bins_zero = as.vector(capacity_edge_dictionary[,time_bin_capacity_cost_to_zero])
        
    } else {
        f_pos = mct@mc_t > 0 
    }    
    f_zero = !f_pos


    capacity_edges_with_positive_capacity =  capacity_edge_dictionary[f_pos]
    capacity_edges_with_zero_capacity = capacity_edge_dictionary[f_zero]
    n_capacity_edges = length(capacity_edges_with_positive_capacity)
	capacity_edges_mean = mc_t_mean[f_pos]
	names(capacity_edges_mean) = capacity_edges_with_positive_capacity

	if(is.null(mc_t_sd)) {
		capacity_edges_mean_over_sd = rep(1,length(capacity_edges_with_positive_capacity))
		names(capacity_edges_mean_over_sd) = capacity_edges_with_positive_capacity
	} else {
    	capacity_edges_sd_over_mean = mc_t_sd[f_pos]/mc_t_mean[f_pos]
    	names(capacity_edges_sd_over_mean) = capacity_edges_with_positive_capacity
    	capacity_edges_sd_over_mean = pmin(pmax(capacity_edges_sd_over_mean,sd_over_mean_lower_threshold),sd_over_mean_upper_threshold)
    	capacity_edges_mean_over_sd = 1/capacity_edges_sd_over_mean
	}



    n_edges_discrete = 8    
    
    upper_fold_max = 8

    bux_add = c(rep(1/n_edges_discrete,n_edges_discrete),rep((upper_fold_max - 1)/n_edges_discrete,n_edges_discrete),k_inf)
    bux_add = rep(bux_add,n_capacity_edges)

    convex_cost_first_derivative = function(x,lambda = 1) {
        

        # the two convex functions are f_1(x) = exp(lambda*(x -1)) - lambda*x + lambda - 1
        #                              f_2(x) = 0.5*((lambda^2)/x + (lambda^2)*x - 2*(lambda^2))

        # they have the first derivatives f_1'(x) = lambda*exp(lambda*(x-1)) - lambda
        #                                  f_2'(x) = -0.5*lambda^2/x^2

        # they have the second derivatives f_1''(x) = lambda^2 * exp(lambda(x-1))
        #                                  f_2''(x) = lambda^2/x^3

        # because we approximate these functions by piecewise linear functions, we need to use the first derivative here


        #y1 = exp(lambda*(x -1)) - lambda*x + lambda - 1
        
        y1 = lambda*exp(lambda*(x-1)) - lambda
		y1 = pmin(y1,max_cost)
        y2 = 0.5*(lambda^2)*(-1/(x^2) + 1 )
        y2[!is.finite(y2)] = -max_cost

        
        return( 0.5*y1 + 0.5*y2)
    }

    n_additional_edges = n_capacity_edges*(2*n_edges_discrete + 1)
    A_net_pw_linear = sparseMatrix(i = c(1:n_capacity_edges),
                                   j = capacity_edges_with_positive_capacity,
                                   x = rep(-1,n_capacity_edges),
                                   dims = c(n_capacity_edges,n_edges))
    
    A_add_pw_linear = sparseMatrix(i = rep(c(1:n_capacity_edges),each = 2*n_edges_discrete + 1),
                                   j = c(1:(n_capacity_edges*(2*n_edges_discrete + 1))),
                                   x = rep(1,n_capacity_edges*(2*n_edges_discrete + 1)),
                                   dims = c(n_capacity_edges,n_capacity_edges*(2*n_edges_discrete + 1)))

    mean_over_sd = capacity_edges_mean_over_sd[as.character(capacity_edges_with_positive_capacity)]
    names(mean_over_sd) = capacity_edges_with_positive_capacity

    if(is.null(time_bin_capacity_cost_to_zero)) {
        c_add_list = lapply(capacity_edges_with_positive_capacity,function(capacity_edge) {
            
            lambda = mean_over_sd[as.character(capacity_edge)]
            
            mc_mean = capacity_edges_mean[as.character(capacity_edge)]

            c_add_pw_linear = mc_mean*c(-max_cost,
                                        convex_cost_first_derivative(seq(from = 1/n_edges_discrete,to = 1-1/n_edges_discrete,length.out = n_edges_discrete - 1),lambda), 
                                        convex_cost_first_derivative(seq(from = 1 + (upper_fold_max - 1)/n_edges_discrete,to = upper_fold_max,length.out = n_edges_discrete),lambda),
                                        max_cost)
            return(c_add_pw_linear)
        })
        names(c_add_list) = capacity_edges_with_positive_capacity
    } else {

        c_add_list = lapply(capacity_edges_with_positive_capacity,function(cap_edge) {
            lambda = mean_over_sd[as.character(cap_edge)]
            mc_mean = capacity_edges_mean[as.character(cap_edge)]

            if(cap_edge %in% capacity_edges_time_bins_zero) {
                c_add_pw_linear = rep(0,2*n_edges_discrete + 1)
            } else {
                c_add_pw_linear = mc_mean*c(-max_cost,
                                            convex_cost_first_derivative( seq(from = 1/n_edges_discrete,to = 1-1/n_edges_discrete,length.out = n_edges_discrete - 1),lambda), 
                                            convex_cost_first_derivative(seq(from = 1 + (upper_fold_max - 1)/n_edges_discrete,to = upper_fold_max,length.out = n_edges_discrete),lambda),
                                            max_cost)
            }

            return(c_add_pw_linear)
        })
        names(c_add_list) = capacity_edges_with_positive_capacity
    }


    c_add = do.call(what = c,args = c_add_list)

    # linear constraint matrix
    
    A_net_linear = linear_constraints$lhs
    A_add_linear = sparseMatrix( i = integer(0), j = integer(0),x = numeric(0),dims = c(nrow(A_net_linear),ncol(A_add_pw_linear)))

    mosek_A = rbind(cbind(A_net_linear,A_add_linear),cbind(A_net_pw_linear,A_add_pw_linear))

    blc_mosek = ifelse(linear_constraints$dir %in% c("==",">="),linear_constraints$rhs,-Inf)
    buc_mosek = ifelse(linear_constraints$dir %in% c("==","<="),linear_constraints$rhs, Inf)
    
    blc_mosek = c(blc_mosek,rep(0,n_capacity_edges))
    buc_mosek = c(buc_mosek,rep(0,n_capacity_edges))
    
    bc_mosek = rbind(blc = blc_mosek,buc = buc_mosek)

    blx_mosek = rep(0,n_edges + n_additional_edges)
    bux_mosek = rep(Inf,n_edges)
    # set upper bound to zero for capacity edges where mct@mc_t is zero
    bux_mosek[capacity_edges_with_zero_capacity] = 0
    bux_mosek = c(bux_mosek,bux_add)

    bx_mosek = rbind(blx=blx_mosek, 
                     bux=bux_mosek)

    linear_edge_costs = mct@network$cost[order(mct@network$ID)]

    mosek_c = c(lambda*linear_edge_costs,mu*c_add)

    


    prob <- list(sense="min")
    prob$c  <- mosek_c
    prob$A  <- mosek_A
    prob$bc <- bc_mosek
    prob$bx <- bx_mosek    

    return(list(mosek_problem = prob,mct = mct))
}






solve_piecewise_linear_mosek_network = function(mct,
                                                mc_t_sd = NULL,
                                                mu = 1000,
                                                lambda = 1,
                                                total_flow = 0.99,
                                                min_total_flow = 0.95,
                                                sd_over_mean_lower_threshold = 0.0125,
                                                sd_over_mean_upper_threshold = 4,
                                                time_bin_capacity_cost_to_zero = NULL) {
    
    cmp_piecewise_linear_mosek = generate_piecewise_linear_mosek_problem(mct = mct,
                                                                         mc_t_sd = mc_t_sd,
                                                                         mu = mu,
                                                                         lambda = lambda,
                                                                         total_flow = total_flow,
                                                                         min_total_flow = min_total_flow,
                                                                         time_bin_capacity_cost_to_zero = time_bin_capacity_cost_to_zero,
                                                                         sd_over_mean_lower_threshold = sd_over_mean_lower_threshold,
                                                                         sd_over_mean_upper_threshold = sd_over_mean_upper_threshold)
    
    mct = cmp_piecewise_linear_mosek$mct
    mosek_sol <- Rmosek::mosek(cmp_piecewise_linear_mosek$mosek_problem)
    
    if( mosek_sol$sol$itr$solsta != 'OPTIMAL') {
      warning('Solution status not optimal: ',mosek_sol$sol$itr$solsta)
      solved_edge_flows = rep(0,nrow(mct@network))
    } else {
      solved_edge_flows = mosek_sol$sol$itr$xx[1:nrow(mct@network)]
    }

    
    if(0) {
        net = mct@network

        # add flow_scaling column. This will transform the capacity edge variables into relative quantities giving the amount of flow relativ to the expected mean flow.
        net$flow_scaling = rep(1,nrow(net))

        max_t = max(mct@network$time1)

        for (t in 1:max_t) {

          f_edges = (net$type1 == "norm_b") & (net$type2 == "norm_f") & (net$time1 == t)

          net$mc1[f_edges]
          net$flow_scaling[f_edges] = mct@mc_t[as.integer(net$mc1[f_edges]),t] * total_flow

        }
        mct@network = net
    }

    net = mct@network
    mcf = tgMCTNetFlow(mct, flows = NULL)
    mcf@edge_flows = solved_edge_flows*net$flow_scaling[order(net$ID)]

    net$flow = solved_edge_flows[net$ID]*net$flow_scaling
    f_cap = net$mc1 == net$mc2 & net$time1 == net$time2
    net_cap = net[f_cap,]

    mc_time_weight = dplyr::summarize(dplyr::group_by(net_cap,mc1,time1),tot_weight = sum(flow))

	mc_index = seq_along(mct@metacell_names)
	names(mc_index) = mct@metacell_names
    mc_t_post = sparseMatrix(i = mc_index[mc_time_weight$mc1],
                         j = as.numeric(mc_time_weight$time1), 
                         x = mc_time_weight$tot_weight)

    mc_t_post = as.matrix(mc_t_post)
    rownames(mc_t_post) = mct@metacell_names
    colnames(mc_t_post) = c(1:ncol(mc_t_post))
    
    mcf@mc_t_infer = mc_t_post

    return(list(mcf = mcf,mct = mct,mosek_solution_status = mosek_sol$sol$itr$solsta,mosek_output = mosek_sol))
}






#' Compute matrices for forward and backwrd propagation using the flows
#'
#' @param mcf mcf flow object
#'
#' @return an updated mcf object
#' @export

mctnetflow_comp_propagation= function(mcf) {
	mctnet = mcf@mct
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

mctnetflow_get_flow_mat = function(mcf, mctnet, time) {
	net = mctnet@network
	if(time == -1) {
		f_t = net$type1 == "norm_f" & net$type2 == "norm_b"
	} else {
		f_t = net$time1 == time & net$time2==time+1

	}
	net$flow  = mcf@edge_flows[net$ID]
	net_t = net[f_t,] 
    
   	flow = as.data.frame(dplyr::summarize(dplyr::group_by(net_t, mc1, mc2),
													tot_flow = sum(flow)))

    df_mc_index = data.frame(mc = mctnet@metacell_names,
                             mc_ind = seq_along(mctnet@metacell_names))

    flow = dplyr::left_join(flow, rename(df_mc_index,mc1 = mc,mc1_index = mc_ind))
    flow = dplyr::left_join(flow, rename(df_mc_index,mc2 = mc,mc2_index = mc_ind))

    mc_mat = sparseMatrix(  i = flow$mc1_index,
							j = flow$mc2_index,
							x = flow$tot_flow,
								dims = c(length(mctnet@metacell_names),length(mctnet@metacell_names)),
								dimnames = list(mctnet@metacell_names,mctnet@metacell_names))


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

mctnetflow_propagate_from_t = function(mcf, t, mc_p)
{
	mctnet = mcf@mct    
	max_t = ncol(mctnet@mc_t)
	step_m = list()
	probs = matrix(0, nrow = nrow(mctnet@mc_t), ncol=ncol(mctnet@mc_t),dimnames = list(mctnet@metacell_names,colnames(mctnet@mc_t)))
	
    probs[names(mc_p),t] = mc_p


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


#' Compute matrix of flows over cell types
#'
#' @param mcf mcf network object
#' @param min_time minimum time point
#' @param max_time maximum time point
#' @param df_mc_type data.frame with column mc and column type
#'
#' @return a list of matrices show flows from type t to t'
#' @export

mctnetflow_get_type_flows = function(mcf,df_mc_type, time = NULL, max_time = NULL,included_types = NULL) {
	if(is.null(time)) {
		time = 1
	}
	if(is.null(max_time)) {
		max_time = ncol(mcf@mct@mc_t)
	}
	mctnet = mcf@mct
	net = mctnet@network
	if(is.null(mcf@edge_flows)) {
		stop("trying to query uninitialized flows")
	}
	net$flow = mcf@edge_flows[net$ID]

	if(is.factor(df_mc_type$type)) {
		all_types = levels(df_mc_type$type)
	} else {
		all_types = as.character(unique(df_mc_type$type))
	}
	
	df_type_index = data.frame(type = all_types,type_index = seq_along(all_types))
	df_mc_type = df_mc_type %>% left_join(df_type_index)

	mct_mats = list()
	for(t in time:(max_time-1)) {

		f_t = net$time1 == t & net$time2 == t+1
		net_t = net[f_t,] 

		net_t = net_t %>% left_join(select(rename(df_mc_type,mc1 = mc,mc1_type_index = type_index),mc1,mc1_type_index)) %>% left_join(select(rename(df_mc_type,mc2 = mc,mc2_type_index = type_index),mc2,mc2_type_index))

		flow = as.data.frame(summarize(group_by(net_t, mc1_type_index, mc2_type_index),
													tot_flow = sum(flow)))

		mct_mat = sparseMatrix( i = flow$mc1_type_index,
								j = flow$mc2_type_index,
								x = flow$tot_flow,
								dims = c(length(all_types),length(all_types)),
								dimnames = list(all_types,all_types))
		
		if(!is.null(included_types)) {
			mct_mat = mct_mat[included_types,included_types]
		}

		mct_mat = as.matrix(mct_mat)
		mct_mats[[t]] = mct_mat
	}

	return(mct_mats)
}




#' Orders metacells by their type and within each type by color
#'
#' @param mc_metadata data.frame with metacell_names as rownames and type and time metadata fields
#' @param metadata_field_type column name of type metadata
#' @param metadata_field_time column name of time metadata
#' @param df_type_rank data.frame with column rank and column metadata_field_type
#' @param included_metacells specify this if only a subset of the metacells should be used
#'
#' @return a named vector containing the rank of each metacell
#' @export
#' 
#' 
mctnetflow_order_metacells_by_type_and_time = function(mc_metadata,metadata_field_type,metadata_field_time,df_type_rank,included_metacells = NULL) {

    rownames(mc_metadata) = mc_metadata$metacell
	if(is.null(included_metacells)) {
		incuded_metacells = rownames(mc_metadata)
	}

	mc_metadata = left_join(mc_metadata[included_metacells,],df_type_rank,by = metadata_field_type)

	mc_rank = rank(1000*mc_metadata[,"rank"] + mc_metadata[,metadata_field_time])
	names(mc_rank) = included_metacells

	return(mc_rank)
}





mctnetflow_solve_and_plot_network = function(mgraph,
                                    mc_t,
                                    mc_t_raw,
                                    scdb_dir,
                                    net_id,
                                    flow_id,
                                    fig_dir,
                                    cell_type_colors,
                                    df_mc_annotation,
                                    mu = 1000,
                                    mc_t_sd = NULL,
                                    mc_proliferation_rate = NULL,
                                    temporal_bin_time = NULL,
                                    T_cost = 1e+5,
                                    metacell_names = NULL,
                                    t_exp = 1
                                    ) {

    mct = mctnetwork_from_mgraph_and_mc_proliferation_rate( mgraph,
								    					        mc_t,
														        mc_t_raw,
														        mc_proliferation_rate = mc_proliferation_rate,
														        temporal_bin_time = temporal_bin_time,
														        t_exp = t_exp,
														        T_cost = T_cost,
														        metacell_names = metacell_names)
    
    cmp_out = solve_piecewise_linear_mosek_network( mct = mct,
													mc_t_sd = mc_t_sd,
													mu = mu)

    
    mosek_status = cmp_out$mosek_solution_status

    if(mosek_status == 'OPTIMAL') {

        mcf = mctnetflow_comp_propagation(cmp_out$mcf)
        cmp_out$mcf = mcf

        scdb_add_mc2tnetwork(mct = cmp_out$mct,net_id = net_id,scdb_dir = scdb_dir)
        scdb_add_mc2tnetflow(mcf = cmp_out$mcf,flow_id = net_id,scdb_dir = scdb_dir)

        png_filename = sprintf("%s/%s_barplot_original_and_inferred_cell_type_distribution.png",fig_dir,flow_id)
        mctnetflow_barplot_cell_type_frequencies_network(mct = cmp_out$mct,mcf = mcf,cell_type_colors = cell_type_colors,df_mc_annotation = df_mc_annotation,png_filename = png_filename)

        cell_type_colors$rank = c(1:nrow(cell_type_colors))

        mc_rank = mctnetflow_order_metacells_by_type_and_time(mc_metadata = df_mc_annotation,
                                                      metadata_field_type = "cell_type",
                                                      metadata_field_time = "metacell_time",
                                                      df_type_rank = select(cell_type_colors,cell_type,rank),
                                                     included_metacells = metacell_names)

        df_mc_annotation = dplyr::left_join(select(df_mc_annotation,metacell,cell_type),cell_type_colors,by = 'cell_type')

        mc_color = df_mc_annotation$color
        names(mc_color) = df_mc_annotation$metacell

        flow_plot_filename = sprintf("%s/%s_network_flow_plot.png",fig_dir,flow_id)
        mctnetwork_plot_net_new(mct = cmp_out$mct,mcf = mcf,flow_thresh = 1e-6,
                        filename = flow_plot_filename,
                        mc_rank = mc_rank,
                        mc_color = mc_color,plot_over_flow = F,edge_w_scale = 1e-3)

        flow_plot_filename = sprintf("%s/%s_network_flow_plot_with_overflow.png",fig_dir,flow_id)
        mctnetwork_plot_net_new(mct = cmp_out$mct,mcf = mcf,flow_thresh = 1e-6,
                        filename = flow_plot_filename,
                        mc_rank = mc_rank,
                        mc_color = mc_color,plot_over_flow = T,edge_w_scale = 1e-3)

    }

    return(list(cmp_out = cmp_out,mosek_status = mosek_status))
}


#' Creates two barplots of original and inferred cell type frequencies
#'
#' @param mct mctnetwork object
#' @param mcf mctnetflow object
#' @param df_mc_annotation data.frame with columns metacell and cell_type
#' @param cell_type_colors data.frame with columns cell_type and color
#' @param png_filename filename for png image
#' @param with_original_frequency boolean, if original frequency should also be plotted
#' @return a named vector containing the rank of each metacell
#' @export
#' 
#' 
mctnetflow_barplot_cell_type_frequencies_network = function(mct,mcf,df_mc_annotation,cell_type_colors,png_filename = NULL,with_original_frequency = T) {

    df_mc_annotation = as.data.frame(df_mc_annotation)
    rownames(df_mc_annotation) = df_mc_annotation$metacell

    mc_t_raw = mct@mc_t_raw
    
    mc_t_raw_n = t(t(mc_t_raw)/colSums(mc_t_raw))
    ct_ag = tgstat::tgs_matrix_tapply(x = t(mc_t_raw_n),index = df_mc_annotation[rownames(mc_t_raw_n),"cell_type"],sum)

    ct_ag_inf = tgstat::tgs_matrix_tapply(x = t(mcf@mc_t_infer),index = df_mc_annotation[rownames(mcf@mc_t_infer),"cell_type"],sum)
    
    ct_to_col = cell_type_colors$color
    names(ct_to_col) = cell_type_colors$cell_type

    ct_ag = ct_ag[cell_type_colors$cell_type[cell_type_colors$cell_type %in% rownames(ct_ag)],]
    ct_ag_inf = ct_ag_inf[cell_type_colors$cell_type[cell_type_colors$cell_type %in% rownames(ct_ag_inf)],]

    if(!is.null(png_filename)) {
        w = ifelse(with_original_frequency,1000,500)
        png(png_filename,w = w,h = 500)
    }
    if(with_original_frequency) {
            layout(matrix(1:2,ncol = 2,nrow = 1),widths = c(400,400))
        barplot(ct_ag,col = ct_to_col[rownames(ct_ag)],main = "Original frequency")
        barplot(ct_ag_inf,col = ct_to_col[rownames(ct_ag_inf)],main = "Inferred frequency")
    } else {
        barplot(ct_ag_inf,col = ct_to_col[rownames(ct_ag_inf)],main = "Inferred frequency")
    }
    if(!is.null(png_filename)) {
        dev.off()
    }

}