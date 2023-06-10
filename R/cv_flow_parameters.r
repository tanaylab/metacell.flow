



crossvalidation_hyperparameter_capacity_cost = function(mct,
                                                        mu_vector = c(10,100,1000,10000),
                                                        mc_t_sd = NULL,
                                                        regularization_log_likelihood = 1e-4) {

    time_points = c(2:(ncol(mct@mc_t) - 1))

    df_cv_parameter = data.frame(   time_bin_out = rep(time_points,length(mu_vector)),
                                    mu_vector = rep(mu_vector,each  = length(time_points)))

    cmp_cv = cmp_leave_one_time_bin_out_cross_validation(df_cv_parameter = df_cv_parameter,
                                                mct = mct,
                                                mc_t_sd = mc_t_sd,
                                                regularization_log_likelihood = regularization_log_likelihood)
    
    cmp_cv$df_likelihood = as.data.frame(cmp_cv$df_likelihood)
    #df_likelihood = cmp_cv$df_likelihood
    #average_l1_distance = df_likelihood %>% group_by()

    #layout(matrix(1:2,ncol = 2,nrow = 1))
    #plot()
    return(cmp_cv)
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
cmp_leave_one_time_bin_out_cross_validation = function(df_cv_parameter,
                                                       mct,
                                                       mc_t_sd = NULL,
                                                       regularization_log_likelihood = 1e-4) {
    
    mc_t_raw_n = t(t(mct@mc_t_raw)/colSums(mct@mc_t_raw))
	input_parameter_list = lapply(1:nrow(df_cv_parameter),function(i) {
		v = c("time_bin_out" = df_cv_parameter$time_bin_out[i],"mu" = df_cv_parameter$mu[i])
		return(v)
	})
	names(input_parameter_list) = c(1:nrow(df_cv_parameter))

	cmp_list = parallel::mclapply(input_parameter_list,function(input_parameter) {

            t_out = as.integer(input_parameter["time_bin_out"])
		    cmp_out = solve_piecewise_linear_mosek_network(mct = mct,
														   mc_t_sd = mc_t_sd,
														   mu = as.numeric(input_parameter["mu"]),
														   time_bin_capacity_cost_to_zero = t_out)

            mc_t_out_infer = cmp_out$mcf@mc_t_infer[,t_out]   

			return(list(mc_t_out_infer = mc_t_out_infer,t_out = t_out,mu = as.numeric(input_parameter["mu"])))
	})


    df_out = do.call(rbind,lapply(cmp_list,function(cv_list) {
        t_out = cv_list$t_out
        mu = cv_list$mu
        mc_t_out_infer = cv_list$mc_t_out_infer
        log_likelihood = sum(log2(mc_t_raw_n[,t_out] + regularization_log_likelihood) * mc_t_raw_n[,t_out] -log2(mc_t_out_infer + regularization_log_likelihood) * mc_t_raw_n[,t_out])
        l1_distance = sum(abs(mc_t_out_infer - mc_t_raw_n[,t_out])) 

        return(c("mu" = mu,"time_bin_out" = as.integer(t_out),"log_likelihood" = log_likelihood,"l1_distance" = l1_distance))
    }))
    

	return(list(df_likelihood = df_out,cmp_list = cmp_list))
}



cv_t_exp_time_bin_out = function(mgraph,
                                    mc_t,
                                    mc_t_raw,
                                    mu,
                                    mc_t_sd = NULL,
                                    mc_proliferation_rate = NULL,
                                    temporal_bin_time = NULL,
                                    T_cost = 1e+5,
                                    metacell_names = NULL,
                                    t_exp_vector = 2^c(-6:6),
                                    regularization_log_likelihood = 1e-4) {
    
    time_points = c(2:(ncol(mc_t) - 1))
    mc_t_raw_n = t(t(mc_t_raw[metacell_names,])/colSums(mc_t_raw[metacell_names,]))
    df_cv_parameter = data.frame(   time_bin_out = rep(time_points,length(t_exp_vector)),
                                    t_exp = rep(t_exp_vector,each  = length(time_points)))

    parameter_index = seq_len(nrow(df_cv_parameter))
    names(parameter_index) = seq_len(nrow(df_cv_parameter))
    cmp_list = parallel::mclapply(parameter_index,function(i) {
        
        t_exp = df_cv_parameter$t_exp[i]
        t_out = df_cv_parameter$time_bin_out[i]

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
														mu = mu,
                                                        time_bin_capacity_cost_to_zero = t_out)

        mc_t_out_infer = cmp_out$mcf@mc_t_infer[,t_out]

		return(list(mc_t_out_infer = mc_t_out_infer,t_out = t_out,t_exp = t_exp,mosek_solution_status = cmp_out$mosek_solution_status,mosek_output = cmp_out$mosek_output))
	})

    df_out = do.call(rbind,lapply(cmp_list,function(cv_list) {
        t_out = cv_list$t_out
        t_exp = cv_list$t_exp
        mc_t_out_infer = cv_list$mc_t_out_infer
        log_likelihood = sum(log2(mc_t_raw_n[,t_out] + regularization_log_likelihood) * mc_t_raw_n[,t_out] -log2(mc_t_out_infer + regularization_log_likelihood) * mc_t_raw_n[,t_out])
        l1_distance = sum(abs(mc_t_out_infer - mc_t_raw_n[,t_out])) 

        return(c("t_exp" = t_exp,"time_bin_out" = as.integer(t_out),"log_likelihood" = log_likelihood,"l1_distance" = l1_distance))
    }))

    return(list(df_likelihood = df_out,cmp_list = cmp_list))
}

cv_cmp_flow_for_t_exp_vector = function(mgraph,
                                        mc_t,
                                        mc_t_raw,
                                        mu,
                                        data_dir,
                                        df_mc_annotation,
                                        cell_type_colors,
                                        net_id,
                                        flow_id,
                                        fig_dir,
                                        mc_t_sd = NULL,
                                        mc_proliferation_rate = NULL,
                                        temporal_bin_time = NULL,
                                        T_cost = 1e+5,
                                        metacell_names = NULL,
                                        t_exp_vector = 2^c(-6:6)) {
    


    if(!dir.exists(data_dir)) {
        stop("data_dir does not exist")
    }
    
    df_parameter = data.frame(  index = seq_len(length(t_exp_vector)),
                                mu = rep(mu,length(t_exp_vector)),
                                t_exp = t_exp_vector)

    parameter_index = seq_len(nrow(df_parameter))
    names(parameter_index) = seq_len(nrow(df_parameter))

    write.table(x = df_parameter,file = paste0(data_dir,'/df_parameter.tsv'),sep = '\t')

    cmp_list = parallel::mclapply(parameter_index,function(i) {
        
        t_exp = df_parameter$t_exp[i]

        net_id = paste('mm_emb_',i)
        flow_id = net_id

        cmp = mctnetflow_solve_and_plot_network(mgraph = mgraph,
                                  mc_t = mc_t,
                                  mc_t_raw = mc_t_raw,
                                  scdb_dir = data_dir,
                                  net_id = net_id,
                                  flow_id = flow_id,
                                  fig_dir = fig_dir,
                                  cell_type_colors = cell_type_colors,
                                  df_mc_annotation = df_mc_annotation,
                                  mu = 1000,
                                  t_exp = t_exp,
                                  mc_proliferation_rate = mc_proliferation_rate,
                                  temporal_bin_time = temporal_bin_time,
                                  metacell_names = metacell_names,
                                  mc_t_sd = mc_t_sd)

        saveRDS(object = cmp,file = sprintf("%s/cmp_out_t_exp_%d.rds",data_dir,i))

		return(cmp$mosek_status)
	})
    df_parameter$mosek_status = do.call(c,cmp_list)
    return(df_parameter)
}

mc_proliferation_from_cell_cycle_expression = function(mc_cc,alpha_exponent,max_proliferation_rate,min_expression = NULL,max_expression = NULL) {

    if(is.null(max_expression)) {
        max_expression = quantile(mc_cc,0.95)
    }
    if(is.null(min_expression)) {
        min_expression = min(mc_cc)
    }

    mc_cc_n = pmax(pmin(mc_cc,max_expression) - min_expression,0)/(max_expression - min_expression)

    mc_proliferation_rate = mc_cc_n^alpha_exponent*max_proliferation_rate

    return(mc_proliferation_rate)
}






compute_mc_t_and_mc_t_sd_from_mc_embryo_mat = function(mc_embryo_mat,
                                                       embryo_to_time_bin,
                                                        included_embryos = NULL,
                                                        included_metacells = NULL,
                                                        zero_threshold = 2e-4) {

    if(is.null(included_metacells)) {
        included_metacells = rownames(mc_embryo_mat)
    }

    if(!is.null(included_embryos)) {
        embryo_to_time_bin = embryo_to_time_bin[included_embryos]
    }

    mc_t = t(tgstat::tgs_matrix_tapply(x = mc_embryo_mat[included_metacells,names(embryo_to_time_bin)],embryo_to_time_bin,mean))
    mc_t_sd = t(tgstat::tgs_matrix_tapply(x = mc_embryo_mat[included_metacells,names(embryo_to_time_bin)],embryo_to_time_bin,sd))

    mc_t[mc_t < zero_threshold] = 0
    mc_t_sd[mc_t < zero_threshold] = 0

    n_t = colSums(mc_t)
    mc_t = t(t(mc_t)/n_t)
    mc_t_sd = t(t(mc_t_sd)/n_t)
    
    return(list(mc_t = mc_t,
                mc_t_sd = mc_t_sd))

}




cv_leave_embryos_out_wrapper = function(mgraph,
                                mc_embryo_mat,
                                mc_t_raw,
                                embryo_to_time_bin,
                                mu,
                                data_dir,
                                mc_proliferation_rate = NULL,
                                temporal_bin_time = NULL,
                                T_cost = 1e+5,
                                metacell_names = NULL,
                                t_exp = 1,
                                included_embryos = NULL,
                                regularization_log_likelihood = 1e-4) {

    n_cv = 10

    cmp_list = lapply(0:(n_cv - 1),function(i) {
        
        if(is.null(included_embryos)) {
            included_embryos = names(embryo_to_time_bin)
        }
        n_embryos  = length(included_embryos)
        embryos_out = included_embryos[c(1:length(included_embryos)) %% n_cv == i]

        embryos_f = setdiff(included_embryos,embryos_out)

        cmp_mc_t = compute_mc_t_and_mc_t_sd_from_mc_embryo_mat(mc_embryo_mat,
                                                        embryo_to_time_bin,
                                                        included_embryos = embryos_f,
                                                        included_metacells = metacell_names)
        
        mc_t = cmp_mc_t$mc_t
        mc_t_sd = cmp_mc_t$mc_t_sd

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

        mc_t_infer = cmp_out$mcf@mc_t_infer

        mc_t_infer_matched = mc_t_infer[,embryo_to_time_bin[embryos_out]]
        colnames(mc_t_infer_matched) = embryos_out
        
        cmp_all = list(mc_t_infer_matched = mc_t_infer_matched,mosek_solution_status = cmp_out$mosek_solution_status)
		return(cmp_all)
	})

    mc_t_infer_matched = do.call(what = cbind,arges = lapply(cmp_list,function(cmp) {
        return(cmp$mc_t_infer_matched)
    }))

    mc_embryo_mat = mc_embryo_mat[,colnames(mc_t_infer_matched)]

    mosek_solution_status = do.call(what = c,arges = lapply(cmp_list,function(cmp) {
        return(cmp$mosek_solution_status)
    }))
    df_mosek_solution_status = data.frame(mosek_solution_status = mosek_solution_status,
                                          cv_index = c(1:n_cv))

    return(list(mc_t_infer_matched = mc_t_infer_matched,mc_embryo_mat = mc_embryo_mat,df_cv_mosek_solution_status = df_mosek_solution_status))    
}


cv_time_bin_out = function(   mgraph,
                                    mc_t,
                                    mc_t_raw,
                                    mu,
                                    time_bin_out,
                                    mc_t_sd = NULL,
                                    mc_proliferation_rate = NULL,
                                    temporal_bin_time = NULL,
                                    T_cost = 1e+5,
                                    metacell_names = NULL,
                                    mct = NULL,
                                    t_exp = 1) {

    if(is.null(mct)) {
        mct = mctnetwork_from_mgraph_and_mc_proliferation_rate( mgraph,
                                                                mc_t,
                                                                mc_t_raw,
                                                                mc_proliferation_rate = mc_proliferation_rate,
                                                                temporal_bin_time = temporal_bin_time,
                                                                t_exp = t_exp,
                                                                T_cost = T_cost,
                                                                metacell_names = metacell_names)       
    }

    cmp_out = solve_piecewise_linear_mosek_network(mct = mct,
                                                    mc_t_sd = mc_t_sd,
                                                    mu = mu,
                                                    time_bin_capacity_cost_to_zero = time_bin_out)

    mc_t_out_infer = cmp_out$mcf@mc_t_infer[,time_bin_out]   

    return(list(mc_t_out_infer = mc_t_out_infer,time_bin_out = time_bin_out,mosek_solution_status = cmp_out$mosek_solution_status))
}



cv_leave_embryos_out = function(mgraph,
                                mc_embryo_mat,
                                mc_t_raw,
                                embryo_to_time_bin,
                                mu,
                                embryos_out,
                                mc_proliferation_rate = NULL,
                                temporal_bin_time = NULL,
                                T_cost = 1e+5,
                                metacell_names = NULL,
                                t_exp = 1,
                                included_embryos = NULL) {

    n_cv = 10

    if(is.null(included_embryos)) {
        included_embryos = names(embryo_to_time_bin)
    }
    n_embryos  = length(included_embryos)

    embryos_f = setdiff(included_embryos,embryos_out)

    cmp_mc_t = compute_mc_t_and_mc_t_sd_from_mc_embryo_mat(mc_embryo_mat,
                                                    embryo_to_time_bin,
                                                    included_embryos = embryos_f,
                                                    included_metacells = metacell_names)
    
    mc_t = cmp_mc_t$mc_t
    mc_t_sd = cmp_mc_t$mc_t_sd

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

    mc_t_infer = cmp_out$mcf@mc_t_infer

    mc_t_infer_matched = mc_t_infer[,embryo_to_time_bin[embryos_out]]
    colnames(mc_t_infer_matched) = embryos_out
    
    cmp_all = list(mc_t_infer_matched = mc_t_infer_matched,mc_embryo_mat_out = mc_embryo_mat[,embryos_out],mosek_solution_status = cmp_out$mosek_solution_status)
    return(cmp_all)
}