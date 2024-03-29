
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

mctnetwork_plot_net = function(mct_id,flow_id, fn, 
						mc_ord = NULL, colors_ordered = NULL,
						propogate=NULL,
						mc_t_score = NULL,
						edge_w_scale=5e-4, 
						flow_thresh = 1e-4,
						w = 2000,h = 2000,
						mc_cex = 0.5,
						dx_back = 0.15, dy_ext = 0.4,
						sigmoid_edge = F, grad_col_edge = F,
						plot_mc_ids = F,miss_color_thresh = 0.5,
						func_deform=NULL,
						score_shades = colorRampPalette(c("lightgray", "gray", "darkgray", "lightpink", "pink", "red", "darkred"))(1000))
{
	if(!is.null(propogate) | !is.null(mc_t_score)) {
		dx_back = 0
		dy_ext = 0
	}
	
  mct = scdb_mctnetwork(mct_id)
  mcf = scdb_mctnetflow(flow_id)
  if(is.null(mct)) {
    stop("cannot find mctnet object ", mct_id, " when trying to plot net flows")
  }
  net = mct@network
  net$flow = mcf@edge_flows
	mc = scdb_mc(mct@mc_id)
	if(is.null(mc)) {
		stop("cannot find mc object ", mct@mc_id, " matching the mc id in the mctnetwork object! db mismatch? recreate objects?")
	}
	if(is.null(mcf@edge_flows) | sum(mcf@edge_flows)==0) {
		stop("flows seems not to be initialized in mct id ", mct_id, " maybe rerun the mincost algorithm?")
	}
	names(mc@colors) = as.character(1:length(mc@colors))
	
#color_ord = read.table("config/atlas_type_order.txt", h=T, sep="\t")
#order MCs by type, mean age
  if(is.null(mc_ord)) {
	 if(is.null(colors_ordered)) {
			stop("specify either mc_ord or color ord when plotting mctnet network")
	 }
	 mc_rank = mctnetwork_mc_rank_from_color_ord(mct_id,flow_id, colors_ordered)
  } else {
    mc_rank = rep(-1,length(mc_ord))
    mc_rank[mc_ord] = c(1:length(mc_ord))
    names(mc_rank) = as.character(1:length(mc_rank))
  }
  
  mc_rank["-2"] = 0
  mc_rank["-1"] = length(mc_rank)/2
  
  #add growth mc

	f= net$flow > flow_thresh
	nn = net[f,]
	x1 = nn$time1
	x2 = nn$time2
	y1 = as.numeric(mc_rank[as.character(nn$mc1)])
	y2 = as.numeric(mc_rank[as.character(nn$mc2)])

	x1 = ifelse(nn$type1 == "growth", x1 + 0.3, x1)
	x2 = ifelse(nn$type2 == "growth", x2 + 0.3, x2)
	x1 = ifelse(nn$type1 == "norm_b" | nn$type1 == "extend_b",x1-dx_back,x1)
	x2 = ifelse(nn$type2 == "norm_b" | nn$type2 == "extend_b",x2-dx_back,x2)
	y1 = ifelse(nn$type1 == "src", max(y1)/2, y1)
	y2 = ifelse(nn$type2 == "sink", max(y2)/2, y2)
	y2 = ifelse(nn$type2 == "sink", NA, y2)
	y1 = ifelse(nn$type1 == "growth", y2+2.5, y1)
	y2 = ifelse(nn$type2 == "growth", y1+2.5, y2)
	y1 = ifelse(nn$type1 == "extend_b" | nn$type1 == "extend_f",y1+dy_ext, y1)
	y2 = ifelse(nn$type2 == "extend_f" | nn$type2 == "extend_b",y2+dy_ext, y2)

	if(!is.null(func_deform)) {
		x1 = ifelse(nn$type1 == "growth", NA, x1)
		x2 = ifelse(nn$type2 == "growth", NA, x2)
		y1 = ifelse(nn$type1 == "growth", NA, y1)
		y2 = ifelse(nn$type2 == "growth", NA, y2)
		min_x = min(c(x1,x2),na.rm=T)
		min_y = min(c(y1,y2),na.rm=T)
		range_x = (max(c(x1,x2),na.rm=T)-min_x)
		range_y = (max(c(y1,y2),na.rm=T)-min_y)
		ax = c((x1-min_x)/range_x, (x2-min_x)/range_x)
		ay = c((y1-min_y)/range_y, (y2-min_y)/range_y)
		atag = func_deform(ax, ay)
		n1 = length(x1)
		n = length(ax)
		x1 = atag[[1]][1:n1]
		x2 = atag[[1]][(n1+1):n]
		y1 = atag[[2]][1:n1]
		y2 = atag[[2]][(n1+1):n]
		
if(0) {
		xy1 = func_deform((x1-min_x)/range_x, (y1-min_y)/range_y)
		xy2 = func_deform((x2-min_x)/range_x, (y2-min_y)/range_y)
		x1 = xy1[[1]]
		x2 = xy2[[1]]
		y1 = xy1[[2]]
		y2 = xy2[[2]]
}
	}

	nn$mc1 = ifelse(nn$type1 == "src", nn$mc2, nn$mc1)
	
	png(fn, width = w,height = h)
	f_overflow = nn$type2=="extend_f" & nn$cost > 100
	f_underflow = nn$type2 == "norm_f" & nn$cost < -100
#  nn$flow/(1e-8 + nn$capacity) < miss_color_thresh
	if(is.null(propogate) & is.null(mc_t_score)) {
		plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
		mc_rgb = col2rgb(mc@colors)/256
		f = nn$mc1>0 & nn$mc2 > 0
		m1 = as.numeric(nn$mc1[f])
		m2 = as.numeric(nn$mc2[f])
		seg_df = data.frame(x1 = x1[f], y1=y1[f], dx=x2[f]-x1[f], dy=y2[f]-y1[f], 
									r1 = mc_rgb["red",m1],
									r2 = mc_rgb["red",m2],
									g1 = mc_rgb["green",m1],
									g2 = mc_rgb["green",m2],
									b1 = mc_rgb["blue",m1],
									b2 = mc_rgb["blue",m2])

		for(alpha in seq(0,0.98,0.02)) {
			beta = alpha
			beta5 = alpha+0.02
			if(sigmoid_edge) {
			beta = plogis(alpha,loc=0.5,scale=0.1)
			beta5 = plogis(alpha+0.02,loc=0.5,scale=0.1)
			}
			sx1 = seg_df$x1+alpha*seg_df$dx
			sx2 = seg_df$x1+(alpha+0.02)*seg_df$dx
			sy1 = seg_df$y1+beta*seg_df$dy
			sy2 = seg_df$y1+beta5*seg_df$dy
			alpha_col = ifelse(grad_col_edge, alpha,0)
			rgb_r = seg_df$r2*alpha_col+seg_df$r1*(1-alpha_col)
			rgb_g = seg_df$g2*alpha_col+seg_df$g1*(1-alpha_col)
			rgb_b = seg_df$b2*alpha_col+seg_df$b1*(1-alpha_col)
			cols = rgb(rgb_r, rgb_g, rgb_b)
			segments(sx1, sy1, sx2, sy2, 
				col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
				lwd=pmin(nn$flow/edge_w_scale, 10))
		}
#		segments(x1,y1,x2,y2, 
#				col=ifelse(nn$type2=="growth", "black", mc@colors[nn$mc1]), 
#				lwd=pmin(nn$flow/edge_w_scale, 10))
		f = f_overflow; segments(x1[f],y1[f],x2[f],y2[f], col="red", 
									lwd=pmin(nn$flow[f]/edge_w_scale, 10))
		f = f_underflow; segments(x1[f],y1[f],x2[f],y2[f], col="blue", 
	         					lwd=pmin((nn$capacity[f] - nn$flow[f])/edge_w_scale,10))
#		points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=1)
	} else if(!is.null(mc_t_score)) {
		plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
		mc_t_score = pmax(mc_t_score, quantile(mc_t_score,0.03))
		mc_t_score = pmin(mc_t_score, quantile(mc_t_score,0.97))
		mc_t_score = mc_t_score-min(mc_t_score)
		mc_t_score = mc_t_score/max(mc_t_score)
		mc_t_score = floor(1+999*mc_t_score)
		f = nn$mc1>0 & nn$mc2>0 & nn$time1>0
		max_mc = nrow(mc_t_score)
		m1 = as.numeric(nn$mc1[f]) 
		m2 = as.numeric(nn$mc2[f]) 
		score1 = rep(1, nrow(nn))
		score2 = rep(1, nrow(nn))
		score1[f] = mc_t_score[m1+(nn[f,"time1"]-1)*max_mc]
		score2[f] = mc_t_score[m2+(nn[f,"time2"]-1)*max_mc]

		seg_df = data.frame(x1 = x1, y1=y1, dx=x2-x1, dy=y2-y1, 
														score1= score1, dscore=score2-score1)

		for(alpha in seq(0,0.95,0.05)) {
			x1 = seg_df$x1+alpha*seg_df$dx
			x2 = seg_df$x1+(alpha+0.05)*seg_df$dx
			y1 = seg_df$y1+alpha*seg_df$dy
			y2 = seg_df$y1+(alpha+0.05)*seg_df$dy
			cols = score_shades[floor(seg_df$score1 + (alpha+0.025)*seg_df$dscore)]
			segments(x1, y1, x2, y2, 
				col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
				lwd=pmin(nn$flow/edge_w_scale, 10))
		}
		rect(seg_df$x1-0.12,seg_df$y1-0.5, seg_df$x1+0.12, seg_df$y1+0.5,
					col = mc@colors[nn$mc1], border=NA)
		rect(seg_df$x2-0.12,seg_df$y2-0.5, seg_df$x2+0.12, seg_df$y2+0.5,
					col = mc@colors[nn$mc2], border=NA)
	} else {
		plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
		max_time = length(propogate)
		m1 = as.numeric(nn$mc1) 
		m2 = as.numeric(nn$mc2) 
		max_m = ncol(propogate[[1]])
		prop_flow = rep(0, nrow(nn))
		for(t in 1:max_time) {
			f = (nn$time1 == t) & nn$mc1>0 & nn$mc2>0
			prop_flow[f] = propogate[[t]][m1[f]+max_m*(m2[f]-1)]
		}
		segments(x1,y1,x2,y2, 
				col=ifelse(nn$type2=="growth", "black", mc@colors[nn$mc1]), 
				lwd=pmin(prop_flow/edge_w_scale, 10))
		points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=m_cex)
	}

	if(plot_mc_ids) {
		f1 = nn$type1!="growth" 
		text(x1[f1]-0.2,y1[f1], labels = nn$mc1[f1], cex=1)
#	  text(c(x1[f1],x2[f2]),c(y1[f1],y2[f2]),labels = c(nn$mc1[f1],nn$mc2[f2]), cex=1)
	}
	
	dev.off()
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

mm_mctnetwork_plot_net = function(mct_id, flow_id, fn, 
                                  mc_ord = NULL, colors_ordered = NULL,
                                  propogate=NULL,
                                  mc_t_score = NULL,
                                  edge_w_scale=5e-4, 
                                  w = 2000,h = 2000,
                                  mc_cex = 0.5,
                                  dx_back = 0.15, dy_ext = 0.4,
                                  sigmoid_edge = F, grad_col_edge = F,
                                  plot_mc_ids = F,miss_color_thresh = 0.5,
                                  func_deform=NULL,
                                  plot_background_as_grey = F,
                                  bg_col = "gray90",
                                  score_shades = colorRampPalette(c("lightgray", "gray", "darkgray", "lightpink", "pink", "red", "darkred"))(1000),
                                  bg_scale = 1,
                                  fr_scale = 2,
                                  max_lwd = 10,
                                  plot_pdf = FALSE,
                                  show_over_under_flow = T,
                                  show_axes = T,
                                  bg = "white")
{
  if(!is.null(propogate) | !is.null(mc_t_score)) {
    dx_back = 0
    dy_ext = 0
  }
  mct = scdb_mctnetwork(mct_id)
  mcf = scdb_mctnetflow(flow_id)
  if(is.null(mct)) {
    stop("cannot find mctnet object ", mct_id, " when trying to plot net flows")
  }
  net = mct@network
  net$flow = mcf@edge_flows
  mc = scdb_mc(mct@mc_id)
  if(is.null(mc)) {
    stop("cannot find mc object ", mct@mc_id, " matching the mc id in the mctnetwork object! db mismatch? recreate objects?")
  }
  if(is.null(mcf@edge_flows) | sum(mcf@edge_flows)==0) {
    stop("flows seems not to be initialized in mct id ", mct_id, " maybe rerun the mincost algorithm?")
  }
  names(mc@colors) = as.character(1:length(mc@colors))
  
  #color_ord = read.table("config/atlas_type_order.txt", h=T, sep="\t")
  #order MCs by type, mean age
  if(is.null(mc_ord)) {
    if(is.null(colors_ordered)) {
      stop("specify either mc_ord or color ord when plotting mctnet network")
    }
    mc_rank = mctnetwork_mc_rank_from_color_ord(mct_id, flow_id, colors_ordered)
  } else {
    mc_rank = rep(-1,length(mc_ord))
    mc_rank[mc_ord] = c(1:length(mc_ord))
    names(mc_rank) = as.character(1:length(mc_rank))
  }
  
  mc_rank["-2"] = 0
  mc_rank["-1"] = length(mc_rank)/2
  
  #add growth mc
  
  f= net$flow > 1e-4
  nn = net[f,]
  x1 = nn$time1
  x2 = nn$time2
  y1 = as.numeric(mc_rank[as.character(nn$mc1)])
  y2 = as.numeric(mc_rank[as.character(nn$mc2)])
  
  x1 = ifelse(nn$type1 == "growth", x1 + 0.3, x1)
  x2 = ifelse(nn$type2 == "growth", x2 + 0.3, x2)
  x1 = ifelse(nn$type1 == "norm_b" | nn$type1 == "extend_b",x1-dx_back,x1)
  x2 = ifelse(nn$type2 == "norm_b" | nn$type2 == "extend_b",x2-dx_back,x2)
  y1 = ifelse(nn$type1 == "src", max(y1)/2, y1)
  y2 = ifelse(nn$type2 == "sink", max(y2)/2, y2)
  y2 = ifelse(nn$type2 == "sink", NA, y2)
  y1 = ifelse(nn$type1 == "growth", y2+2.5, y1)
  y2 = ifelse(nn$type2 == "growth", y1+2.5, y2)
  y1 = ifelse(nn$type1 == "extend_b" | nn$type1 == "extend_f",y1+dy_ext, y1)
  y2 = ifelse(nn$type2 == "extend_f" | nn$type2 == "extend_b",y2+dy_ext, y2)
  
  if(!is.null(func_deform)) {
    min_x = min(c(x1,x2),na.rm=T)
    min_y = min(c(y1,y2),na.rm=T)
    range_x = (max(c(x1,x2),na.rm=T)-min_x)
    range_y = (max(c(y1,y2),na.rm=T)-min_y)
    xy1 = func_deform((x1-min_x)/range_x, (y1-min_y)/range_y)
    xy2 = func_deform((x2-min_x)/range_x, (y2-min_y)/range_y)
    x1 = xy1[[1]]
    x2 = xy2[[1]]
    y1 = xy1[[2]]
    y2 = xy2[[2]]
  }
  
  nn$mc1 = ifelse(nn$type1 == "src", nn$mc2, nn$mc1)
  
  if(plot_pdf) {
    pdf(fn,width = w,height =h,useDingbats = F)
  } else {
    png(fn, width = w,height = h,bg = bg)
  }
  
  f_overflow = nn$type2=="norm_f" & nn$cost > 100
  f_underflow = nn$type2 == "norm_f" & nn$cost < -100
  #  nn$flow/(1e-8 + nn$capacity) < miss_color_thresh
  if(is.null(propogate) & is.null(mc_t_score)) {
    
    if(show_axes) {
      plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
    } else {
      plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex,axes = F,xlab = "",ylab = "")
    }
    
    mc_rgb = col2rgb(mc@colors)/256
    f = nn$mc1>0 & nn$mc2 > 0
    m1 = as.numeric(nn$mc1[f])
    m2 = as.numeric(nn$mc2[f])
    seg_df = data.frame(x1 = x1[f], y1=y1[f], dx=x2[f]-x1[f], dy=y2[f]-y1[f], 
                        r1 = mc_rgb["red",m1],
                        r2 = mc_rgb["red",m2],
                        g1 = mc_rgb["green",m1],
                        g2 = mc_rgb["green",m2],
                        b1 = mc_rgb["blue",m1],
                        b2 = mc_rgb["blue",m2])
    
    for(alpha in seq(0,0.98,0.02)) {
      beta = alpha
      beta5 = alpha+0.02
      if(sigmoid_edge) {
        beta = plogis(alpha,loc=0.5,scale=0.1)
        beta5 = plogis(alpha+0.02,loc=0.5,scale=0.1)
      }
      sx1 = seg_df$x1+alpha*seg_df$dx
      sx2 = seg_df$x1+(alpha+0.02)*seg_df$dx
      sy1 = seg_df$y1+beta*seg_df$dy
      sy2 = seg_df$y1+beta5*seg_df$dy
      alpha_col = ifelse(grad_col_edge, alpha,0)
      rgb_r = seg_df$r2*alpha_col+seg_df$r1*(1)
      rgb_g = seg_df$g2*alpha_col+seg_df$g1*(1)
      rgb_b = seg_df$b2*alpha_col+seg_df$b1*(1)
      cols = rgb(rgb_r, rgb_g, rgb_b)
      segments(sx1, sy1, sx2, sy2, 
               col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
               lwd=pmin(nn$flow/edge_w_scale, max_lwd))
    }
    if(show_over_under_flow) {
      f = f_overflow; segments(x1[f],y1[f],x2[f],y2[f], col="red", 
                               lwd=pmin(nn$flow[f]/edge_w_scale, max_lwd))
      f = f_underflow; segments(x1[f],y1[f],x2[f],y2[f], col="blue", 
                                lwd=pmin((nn$capacity[f] - nn$flow[f])/edge_w_scale,max_lwd))
    }
    
    points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=1)
    if(show_over_under_flow) {
      f = f_overflow; segments(x1[f],y1[f],x2[f],y2[f], col="red", 
                               lwd=pmin(nn$flow[f]/edge_w_scale, max_lwd))
      f = f_underflow; segments(x1[f],y1[f],x2[f],y2[f], col="blue", 
                                lwd=pmin((nn$capacity[f] - nn$flow[f])/edge_w_scale,max_lwd))
    }
    
  } else if(!is.null(mc_t_score)) {
    plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
    mc_t_score = pmax(mc_t_score, quantile(mc_t_score,0.03))
    mc_t_score = pmin(mc_t_score, quantile(mc_t_score,0.97))
    mc_t_score = mc_t_score-min(mc_t_score)
    mc_t_score = mc_t_score/max(mc_t_score)
    mc_t_score = floor(1+999*mc_t_score)
    f = nn$mc1>0 & nn$mc2>0 & nn$time1>0
    max_mc = nrow(mc_t_score)
    m1 = as.numeric(nn$mc1[f]) 
    m2 = as.numeric(nn$mc2[f]) 
    score1 = rep(1, nrow(nn))
    score2 = rep(1, nrow(nn))
    score1[f] = mc_t_score[m1+(nn[f,"time1"]-1)*max_mc]
    score2[f] = mc_t_score[m2+(nn[f,"time2"]-1)*max_mc]
    
    seg_df = data.frame(x1 = x1, y1=y1, dx=x2-x1, dy=y2-y1, 
                        score1= score1, dscore=score2-score1)
    
    for(alpha in seq(0,0.95,0.05)) {
      x1 = seg_df$x1+alpha*seg_df$dx
      x2 = seg_df$x1+(alpha+0.05)*seg_df$dx
      y1 = seg_df$y1+alpha*seg_df$dy
      y2 = seg_df$y1+(alpha+0.05)*seg_df$dy
      cols = score_shades[floor(seg_df$score1 + (alpha+0.025)*seg_df$dscore)]
      segments(x1, y1, x2, y2, 
               col=ifelse(nn$type2=="growth" | nn$type1=="source" | nn$type2=="sink", "gray", cols), 
               lwd=pmin(nn$flow/edge_w_scale, max_lwd))
    }
    rect(seg_df$x1-0.12,seg_df$y1-0.5, seg_df$x1+0.12, seg_df$y1+0.5,
         col = mc@colors[nn$mc1], border=NA)
    rect(seg_df$x2-0.12,seg_df$y2-0.5, seg_df$x2+0.12, seg_df$y2+0.5,
         col = mc@colors[nn$mc2], border=NA)
  } else {
    
    #mc_col_bg = alpha(c(c("gray","gray"),mc@colors),0.01)
    #names(mc_col_bg) = c(c("-2","-1"),as.character(c(1:length(mc@colors))))
    
    
    plot(c(x1,x2), c(y1,y2), pch=19, col=bg_col,cex=mc_cex*bg_scale*0.2,axes = F,xlab = "",ylab = "")
    
    segments(x1,y1,x2,y2,
             col=bg_col,
             lwd=pmin(nn$flow/edge_w_scale, max_lwd)*bg_scale)
    
    #plot(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex)
    max_time = length(propogate)
    m1 = as.numeric(nn$mc1) 
    m2 = as.numeric(nn$mc2) 
    max_m = ncol(propogate[[1]])
    prop_flow = rep(0, nrow(nn))
    for(t in 1:max_time) {
      f = (nn$time1 == t) & nn$mc1>0 & nn$mc2>0
      prop_flow[f] = propogate[[t]][m1[f]+max_m*(m2[f]-1)]
    }
    
    prop_flow[prop_flow/edge_w_scale < 0.1] = 0
    
    segments(x1,y1,x2,y2, 
             col=ifelse(nn$type2=="growth", "black", mc@colors[nn$mc1]), 
             lwd=pmin(prop_flow/edge_w_scale, max_lwd)*fr_scale)
    points(c(x1,x2), c(y1,y2), pch=19, col=mc@colors[c(nn$mc1,nn$mc2)],cex=mc_cex*rep(pmin(prop_flow/edge_w_scale, max_lwd),2)*0.1*fr_scale)
  }
  
  if(plot_mc_ids) {
    f1 = nn$type1!="growth" 
    text(x1[f1]-0.2,y1[f1], labels = nn$mc1[f1], cex=1)
    #	  text(c(x1[f1],x2[f2]),c(y1[f1],y2[f2]),labels = c(nn$mc1[f1],nn$mc2[f2]), cex=1)
  }
  
  dev.off()
}
