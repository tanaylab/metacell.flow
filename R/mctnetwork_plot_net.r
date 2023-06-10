
deform_snail = function(x,y)
{
 	theta = 3.14159 * (1-x)
	expand = 0.4 + x*1.8
	yoff = y - 0.5

	xtag = cos(theta) + yoff*expand*cos(theta)
	ytag = sin(theta) + yoff*expand*sin(theta)
	return(list(xtag, ytag))
}

deform_boom = function(x,y,rmin=0)
{
	r = rmin+(x)*(1-rmin)
	xc = as.character(x)
	normy = rank(rank(x)*max(y,na.rm=T)+y)
#	normy = y
	min_y = tapply(normy,xc,min, na.rm=T)
	max_y = tapply(normy,xc,max, na.rm=T)
	range_y = max_y-min_y+8
	range_y[range_y==0] = 1
	theta = (normy-min_y[xc])/range_y[xc]*2*3.14159
	
	xtag = r*cos(theta)
	ytag = r*sin(theta)
	return(list(xtag, ytag))
}

mctnetwork_mc_rank_from_color_ord = function(mct_id, flow_id, colors_ordered = NULL)
{
  mct = scdb_mctnetwork(mct_id)
  if(is.null(mct)) {
    stop("cannot find mctnet object ", mct_id, " when trying to plot net flows")
  }
  mcf = scdb_mctnetflow(flow_id)
  if(is.null(mcf)) {
    stop("cannot find mcfnet object ", flow_id, " when trying to plot net flows")
  }
  net = mct@network
  net$flow = mcf@edge_flows
  mct@network = net
  
  mc = scdb_mc(mct@mc_id)
  mc_mean_age = apply(mct@mc_t, 1, function(x) {return(mean(x*(1:length(x))/sum(x))) })
  #   mc_mean_age = 1:length(mc@colors)
  col_to_rank = c(1:length(colors_ordered))
  names(col_to_rank) = colors_ordered
  
  mc_rank = 1000*col_to_rank[mc@colors] + mc_mean_age
  mc_rank = rank(mc_rank)
  names(mc_rank) = as.character(1:length(mc_rank))
  
  mc_rank1 = mc_rank
  
  #cluster flow graph
  #rank clusters by mean rank
  #rank mc by membership + time
  #plot(sign(dfore)*log2(abs(dfore)), sign(dback)*log2(abs(dback)), pch=19, col=mc@colors, cex=1.5)
  return(mc_rank1)
}


mctnetwork_plot_net_new = function(mct,
                                   mcf,
                                   filename,
                                   mc_color = NULL,
						           mc_rank = NULL, 
						           edge_w_scale=5e-4,
                                   max_lwd = 10,
						           flow_thresh = 1e-4,
						           w = 2000,h = 2000,
						           mc_cex = 0.5,
						           dx_back = 0.15,
                                   plot_mc_ids = F,
                                   plot_over_flow = F,
                                   mc_t_infer = NULL) {
    

  
    if(is.null(mct)) {
        stop("mct object is NULL")
    }
    net = mct@network
    net$flow = mcf@edge_flows[net$ID]

	if(is.null(mcf@edge_flows) | sum(mcf@edge_flows)==0) {
		stop("flows seems not to be initialized in mcf maybe rerun the mincost algorithm?")
	}

    if(is.null(mc_color)) {
        mc_color = rep('gray',length(mct@metacell_names))
        names(mc_color) = mct@metacell_names
    }
    if(is.null(names(mc_color))) {
        stop("mc_color needs to have metacell_names as names")
    }
	
    if(is.null(mc_rank)) {
        mc_rank = c(1:length(mct@metacell_names))
        names(mc_rank) = mct@metacell_names
    }
    mc_rank["src"] = length(mc_rank)/2
    mc_rank["sink"] = length(mc_rank)/2
    
  
    #add growth mc
    
    net$edge_color = ifelse(net$type1 == "sink",mc_color[as.character(net$mc2)],mc_color[as.character(net$mc1)])
    
    if(plot_over_flow) {
        
        if(is.null(mc_t_infer)) {
            mc_t_infer = mcf@mc_t_infer
        }

        max_deviation = 0.5
        reg = 1e-3
        mc_relative_deviation = (mc_t_infer - mct@mc_t)/(mct@mc_t + reg)

        df_dev = mc_relative_deviation %>% as.data.frame %>% tibble::rownames_to_column(var = 'mc1') %>% tidyr::pivot_longer(cols = !mc1,names_to = "time1",values_to = "mc_dev")
        df_dev$time1 = as.double(df_dev$time1)
        
        f_cap = net$mc1 == net$mc2 & net$time1 == net$time2
        net_cap = net[f_cap,]
        net_cap = left_join(net_cap,df_dev,by = c('mc1',"time1"))

        shades = rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(101))
        capacity_edge_color = shades[floor(100*pmax(pmin((net_cap$mc_dev + max_deviation)/2/max_deviation,1),0)) + 1]

        net$edge_color[f_cap] = capacity_edge_color

    }
    f= net$flow > flow_thresh    
	nn = net[f,]
    
	x1 = nn$time1
	x2 = nn$time2
	y1 = as.numeric(mc_rank[as.character(nn$mc1)])
	y2 = as.numeric(mc_rank[as.character(nn$mc2)])

	x1 = ifelse(nn$type1 == "norm_b",x1-dx_back,x1)
	x2 = ifelse(nn$type2 == "norm_b",x2-dx_back,x2)
	y1 = ifelse(nn$type1 == "src", max(y1)/2, y1)
	y2 = ifelse(nn$type2 == "sink", max(y2)/2, y2)
	y2 = ifelse(nn$type2 == "sink", NA, y2)

	nn$mc1 = ifelse(nn$type1 == "src", nn$mc2, nn$mc1)

	
	png(filename, width = w,height = h)

    f_cap = nn$mc1 == nn$mc2 & nn$time1 == nn$time2

    plot(x = c(x1[f_cap],x2[f_cap],c(-1)), y = c(y1[f_cap],y2[f_cap],max(y2)/2), pch=19, 
         col=c(mc_color[c(nn$mc1[f_cap],nn$mc2[f_cap])],"gray"),cex=c(pmin(nn$flow[f_cap]/edge_w_scale, max_lwd),1)/10)

    #plot(c(x1,x2), c(y1,y2), pch=19, col=mc_color[c(nn$mc1,nn$mc2)],cex=mc_cex)


    segments(x0 = x1,x1 = x2,y0 = y1,y1 = y2,col  =  nn$edge_color,lwd = pmin(nn$flow/edge_w_scale, max_lwd))
    #points(c(x1,x2), c(y1,y2), pch=19, col=mc_color[c(nn$mc1,nn$mc2)],cex=1)
    points(c(x1[f_cap],x2[f_cap],c(-1)), c(y1[f_cap],y2[f_cap],max(y2)/2), pch=19, col=c(mc_color[c(nn$mc1[f_cap],nn$mc2[f_cap])],"gray"),cex=c(pmin(nn$flow[f_cap]/edge_w_scale, max_lwd),1)/10)

	if(plot_mc_ids) {
				text(x1[f1]-0.2,y1[f1], labels = nn$mc1[f1], cex=1)
	}	
	dev.off()
}

write_flow = function(file_name,mct_id, flow_id, mcnames, flow_thresh = 1e-4){
    mct = scdb_mctnetwork(mct_id)
    mcf = scdb_mctnetflow(flow_id)
    net = mct@network
    net$flow = mcf@edge_flows
    mc = scdb_mc(mct@mc_id)
    f = net$flow > flow_thresh
    nn = net[f,]
    named_nn = nn
    named_nn$mc1[named_nn$type1=="src"] = "src"
    named_nn$mc2[named_nn$type2=="sink"] = "sink"
    f_not_source = named_nn$mc1!="src"
    named_nn$mc1[f_not_source] = mcnames[named_nn$mc1[f_not_source]]
    f_not_sink = named_nn$mc2!="sink"
    named_nn$mc2[f_not_sink] = mcnames[named_nn$mc2[f_not_sink]]
    
    data.table::fwrite(named_nn[,c("mc1", "mc2", "time1", "time2", "type1", "type2", "flow")], file_name)
}

mctnetwork_plot_propagation = function(mct_id, time, p_anchor, fn, 
						mc_ord = NULL, colors_ordered = NULL,
						edge_w_scale=5e-4, 
						w = 2000,h = 2000, mc_cex=1)
{
	mct = scdb_mctnetwork(mct_id)
	if(is.null(mct)) {
		stop("cannot find mctnet object ", mct_id, " when trying to plot net flows anchors")
	}
	probs = mctnetwork_propogate_from_t(mct, time, p_anchor)
	mctnetwork_plot_net(mct_id, fn=fn, 
					propogate = probs$step_m,
					mc_ord = mc_ord, colors_ordered = colors_ordered,
					edge_w_scale = edge_w_scale, 
					w = w, h = h, mc_cex = mc_cex)
}
