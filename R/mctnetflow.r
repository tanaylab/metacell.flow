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
			if(mct@net$ID != names(flows)) {
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
mcell_new_mctnetflow_old = function(flow_id, net_id,scdb_dir,
							init_mincost = F, flow_tolerance=0.01,
							max_flow_tolerance = 0.05,
							flows=NULL)
{
	mcf = tgMCTNetFlow(net_id, flows = flows)
	if(is.null(flows) & init_mincost) {
		mcf = mctnetflow_gen_mincost(mcf, flow_tolerance,max_flow_tolerance)
	}
	scdb_add_mc2tnetflow(mcf = mcf,flow_id = flow_id,scdb_dir = scdb_dir)
}



generate_piecewise_linear_mosek_problem = function(mct,
                                                   mc_t_sd = NULL,
                                                   total_flow = 0.99,
                                                   min_total_flow = 0.95,
                                                   lambda = 1,
                                                   mu = 1,
                                                   total_flow_scaling = 1,
                                                   sd_over_mean_upper_threshold = 4,
                                                   sd_over_mean_lower_threshold = 0.0125,
                                                   time_bin_capacity_cost_to_zero = NULL) {

    k_inf = 1e+5
    max_cost = 1e+5
    n_edges = nrow(mct@network)

	mc_t_mean = mct@mc_t



    # linear constraints of the network

    cmp_constraints = mctnetwork_linear_constraints(mct = mct,
                                                       total_flow = total_flow*total_flow_scaling,
                                                       min_total_flow = min_total_flow*total_flow_scaling)

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



    n_edges_discrete = 50
    
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

#' Compute network flow solutions for different values of mu and where one time bin is left out
#' The output is a list where each element is the return value of the function solve_piecewise_linear_mosek_network for the different input parameter values of mu and time_bin_out
#'
#'
#' @param df_cv_parameter data.frame with column mu and column time_bin_out
#' @param mct mct object that is the output 
#' @param mc_t_sd optional 
#' 
#' @return list where each element is the return value of the function solve_piecewise_linear_mosek_network for the different input parameter values of mu and time_bin_out
#' @export
#' 
cmp_leave_one_time_bin_out_cross_validation = function(df_cv_parameter,mct,mc_t_sd = NULL) {

	input_parameter_list = lapply(1:nrow(df_cv_parameter),function(i) {
		v = c("time_bin_out" = df_cv_parameter$time_bin_out[i],"mu" = df_cv_parameter$mu[i])
		return(v)
	})
	names(input_parameter_list) = c(1:nrow(df_cv_parameter))

	cmp_list = parallel::mclapply(input_parameter_list,function(input_parameter) {

		    cmp_out = solve_piecewise_linear_mosek_network(mct = mct,
														   mc_t_sd = mc_t_sd,
														   mu = as.numeric(input_parameter["mu"]),
														   time_bin_capacity_cost_to_zero = as.integer(input_parameter["time_bin_out"]))
														   

			return(cmp_out)
	})
	return(list(df_cv_parameter = df_cv_parameter,cmp_list = cmp_list))
}



solve_piecewise_linear_mosek_network = function(mct,
                                                mc_t_sd = NULL,
                                                mu = 1000,
                                                lambda = 1,
                                                total_flow = 0.99,
                                                min_total_flow = 0.95,
                                                total_flow_scaling = 1,
                                                sd_over_mean_lower_threshold = 0.0125,
                                                sd_over_mean_upper_threshold = 4,
                                                time_bin_capacity_cost_to_zero = NULL) {
    
    cmp_piecewise_linear_mosek = generate_piecewise_linear_mosek_problem(mct = mct,
                                                                         mc_t_sd = mc_t_sd,
                                                                         mu = mu,
                                                                         lambda = lambda,
                                                                         total_flow = total_flow,
                                                                         min_total_flow = min_total_flow,
                                                                         total_flow_scaling = total_flow_scaling,
                                                                         time_bin_capacity_cost_to_zero = time_bin_capacity_cost_to_zero,
                                                                         sd_over_mean_lower_threshold = sd_over_mean_lower_threshold,
                                                                         sd_over_mean_upper_threshold = sd_over_mean_upper_threshold)
    
    mct = cmp_piecewise_linear_mosek$mct
    mosek_sol <- mosek(cmp_piecewise_linear_mosek$mosek_problem)
    
    if( mosek_sol$sol$itr$solsta != 'OPTIMAL') {
      warning('Solution status not optimal: ',mosek_sol$sol$itr$solsta)  
    }

    solved_edge_flows = mosek_sol$sol$itr$xx[1:nrow(mct@network)]

    total_flow = total_flow_scaling*0.99
    
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
    mcf@edge_flows = solved_edge_flows*net$flow_scaling[order(net$ID)]/total_flow_scaling

    net$flow = solved_edge_flows[net$ID]*net$flow_scaling/total_flow_scaling
    f_cap = net$mc1 == net$mc2 & net$time1 == net$time2
    net_cap = net[f_cap,]

    mc_time_weight = summarize(group_by(net_cap,mc1,time1),tot_weight = sum(flow))

	mc_index = seq_along(mct@metacell_names)
	names(mc_index) = mct@metacell_names
    mc_t_post = sparseMatrix(i = mc_index[mc_time_weight$mc1],
                         j = as.numeric(mc_time_weight$time1), 
                         x = mc_time_weight$tot_weight)

    mc_t_post = as.matrix(mc_t_post)
    rownames(mc_t_post) = mct@metacell_names
    colnames(mc_t_post) = c(1:ncol(mc_t_post))
    
    mcf@mc_t_infer = mc_t_post

    return(list(mcf = mcf,mct = mct,mosek_solution_status = mosek_sol$sol$itr$solsta))
}











#' Old - not updated to mc2 and convex flow version. Compute mincost flow and update the mct objects with flows and 
#' inferred mc capacity per time 
#'
#' @param mcf flow object
#' @param flow_tolerance how much flow we should miss per time
#'
#' @return an updated flow object
#' @export
mctnetflow_gen_mincost_old = function(mcf, flow_tolerance=0.01,max_flow_tolerance = 0.05)
{
	if (max_flow_tolerance < flow_tolerance | max_flow_tolerance >= 1) {
	  stop("max_flow_tolerance must be larger than flow_tolerance (and smaller than 1)")
	}
  mctnet = scdb_mc2tnetwork(mcf@net_id)
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

mctnetflow_get_flow_mat = function(mcf, mctnet, time)
{
	net = mctnet@network
	if(time == -1) {
		f_t = net$type1 == "norm_f" & net$type2 == "norm_b"
	} else {
		f_t = net$time1 == time & net$time2==time+1

	}

	net$flow  = mcf@edge_flows[net$ID]
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

   	mc_mat = mc_mat[mct@metacell_names,mct@metacell_names]	
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
	mctnet = scdb_mc2tnetwork(mcf@net_id)
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
mctnetflow_gen_gmod_mc_t_old = function(mcf_id, mat_id,  genes, min_percentile=0.05)
{
	mcf = scdb_mc2tnetflow(mcf_id)
	if(is.null(mcf)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	mctnet = scdb_mc2tnetwork(mcf@mct_id)
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
mctnetflow_gen_e_gmt_p_old = function(mcf_id, mat_id,  genes)
{
	mcf = scdb_mc2tnetflow(mcf_id)
	if(is.null(mcf)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	mctnet = scdb_mc2tnetwork(mcf@mct_id)
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

mctnetflow_get_egc_on_cluster_transition_old = function(mcf, min_time, max_time, type1, type2, mc_type=NULL)
{
	mctnet = scdb_mc2tnetwork(mcf@mct_id)
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

	if(is.null(included_metacells)) {
		incuded_metacells = rownames(mc_metadata)
	}
	mc_metadata = left_join(mc_metadata[included_metacells,],df_type_rank,by = metadata_field_type)

	mc_rank = rank(1000*mc_metadata[,"rank"] + mc_metadata[,metadata_field_time])
	names(mc_rank) = included_metacells

	return(mc_rank)
}