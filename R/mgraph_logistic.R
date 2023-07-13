
#' Compute metacell manifold graph using logistic distances and balanced K-nn
#'
#' @param mgraph_id id of new object
#' @param mc_id meta cell id to work with
#' @param ignore_edges provide a data frame with mc1,mc2 pairs of edges to delete manually
#' @param feats_gset gene set name for use for computing distances
#' @param features_exclude list of genes to exclude from the features gene set
#' @param logist_loc the "location" parametr of the logistic function used to determine parametric distances between metacelles
#' @param logist_scale the "lscale" parametr of the logistic function used to determine parametric distances between metacelles
#' 
#' @return data.frame with mc1, mc2 and dist column 
#'
#' @export
metacell2_mgraph_logistic = function(mc2_anndata,
                                 feature_genes,
			                     ignore_edges = NULL, 
			                     features_exclude =NULL,
                                 included_metacells = NULL,
			                     logist_loc = 1, logist_scale = 0.2, 
			                     logist_eps = 4e-5, max_d_fold = 3,
                                 max_deg = 4)
{

	if(!is.null(features_exclude)) {
		feature_genes = setdiff(feature_genes, features_exclude)
	}
    feature_genes = intersect(feature_genes,mc2_anndata$var_names)
	message("got ", length(feature_genes), " feat genes for mc graph construction")

    f_gene_incl = mc2_anndata$var_names %in% feature_genes
    if(!is.null(included_metacells)) {
        f_incl = mc2_anndata$obs_names %in% included_metacells
        mc_gene_mat = mc2_anndata$X[f_incl,f_gene_incl]
    } else {
        mc_gene_mat = mc2_anndata$X[,f_gene_incl]
    }

	mgraph = mgraph_comp_logist(mc_gene_mat,
                                feature_genes,
                                logist_loc,
                                logist_scale, 
                                logist_eps, 
                                max_d_fold)

	if(!is.null(ignore_edges)) {
		all_e = paste(mgraph$mc1, mgraph$mc2, sep="-")
		ig_e = paste(ignore_edges$mc1, ignore_edges$mc2, sep="-")
		ig_re = paste(ignore_edges$mc2, ignore_edges$mc1, sep="-")
		f = all_e %in% c(ig_e,ig_re)
		mgraph= mgraph[!f,]
		message("igoring ", sum(f), " edges")
	}
	return(mgraph)
}

#' Compute logistic distances between metacells
#'
#' @param mc_gene_mat expression matrix metacells (rows) over genes (columns). mc2$X entry of anndata metacell2 object.
#' @param feature_genes feature genes used for logistic distance computation
#' @param loc location parameter of plogis function
#' @param scale scale parameter of plogis function
#' @param eps regularization used for logarithmic expression matrix log2(t(mc_gene_mat[,feature_genes] + eps)
#' @param max_d_fold threshold for cutoff of logistic distance compared to second-nearest neighbor
#' @param max_deg maximal degree of resulting manifold graph
#' 
#' @return data.frame with mc1, mc2 and dist column
#' @export
mgraph_comp_logist_mc2 = function(mc_gene_mat, feature_genes, loc = 1, scale = 0.2, eps = 4e-5, max_d_fold = 3, max_deg = 4)
{

    metacell_names = rownames(mc_gene_mat)
    if(is.null(metacell_names)) {
        warning("no rownames specified for mc_gene_mat")
    }
	
    legc = log2(t(as.matrix(mc_gene_mat[,feature_genes]) + eps))
    legc = as.matrix(legc)

	logist_d = function(x) {
		d = abs(legc - x)
		d = plogis(d, loc, scale)
        d = pmax(d  - plogis(0,loc,scale),0)
		return(colSums(d))
	}
	a = apply(legc, 2, logist_d) 

#connect - d-best outgoing. filter by d_best_ratio < 2
	diag(a) = 1000;
	d_T = apply(a, 1, function(x) sort(x,partial=2)[2])
	a_n = a/d_T
	diag(a) = 0;
	diag(a_n) = 0;

   rank_fr = t(apply(a, 1, rank))
   rank_fr_m = rank_fr
   rank_fr_m[a_n > max_d_fold] = 1000
   rank_fr_m2 = rank_fr_m
   rank_fr_m2[t(a_n) > max_d_fold] = 1000

   edges = apply(rank_fr_m2, 1, function(x) {
                        mc2 = metacell_names[which(x > 0 & x <= max_deg+1)];
                        mc1 = metacell_names[rep(which.min(x), length(mc2))];
                        return(data.frame(mc1 = mc1, mc2=mc2)) })
   ed_df = as.data.frame(do.call(rbind, edges))
	ed_df$dist = apply(ed_df, 1, function(x) 1+a[x[1],x[2]])
	return(list(mgraph = ed_df,logist_dist_raw = a))
}


