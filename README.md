
<!-- README.md is generated from README.Rmd. Please edit that file -->
# metacell.flow

<!-- badges: start -->
<!-- badges: end -->
The goal of metacell.flow is to infer differentiation flow models from single cell RNA-seq data and metacell models representing it, with time tags per cell. Potential techniques for inferring or sampling time tags will be under invesigation. This package is still in research mode.

## Installation

You can install the released version of metacell.flow from github with:

\`\`\` r remotes::install\_github("tanaylab/metacell.flow")

## Example

``` r
library("metacell")
library("metacell.flow")
scdb_init("scrna_db/", force_reinit=T)
scdb_flow_init()
tgconfig::override_params("config/netflow.yaml","metacell")

example_leak_table_sing_emb = function(mc_id, leak_emb_endo, leak_exe_endo) 
{
  mc = scdb_mc(mc_id)
  legc = log2(mc@e_gc + 1e-5)
  mc_leak = rep(0,ncol(legc))
  
  # first separation embryonic endoderm (including node/notochord) from meso/-ectoderm 
  x1 = -16
  y1 = -12
  x2 = -12
  y2 = -16
  
  b_emb_endo = (y2 - y1)/(x2 - x1)
  a_emb_endo = (y1*x2 - y2*x1)/(x2 - x1)
  
  f_endo = (legc["Foxa1",]  > a_emb_endo + b_emb_endo*legc["Foxa2",])
  mc_leak[f_endo] = leak_emb_endo
  
  # second separation extraembryonic from embryonic endoderm uses Apoe
  x1 = -8.4
  y1 = -14
  x2 = -11
  y2 = -8.4
  
  b_exe_endo = (y2 - y1)/(x2 - x1)
  a_exe_endo = (y1*x2 - y2*x1)/(x2 - x1)
  f_exe = (legc["Ttr",]  > a_exe_endo + b_exe_endo*legc["Apoe",])
  mc_leak[f_exe] = leak_exe_endo
  return(mc_leak)
}

mc_leak = example_leak_table_sing_emb(
            mc_id = "sing_emb_wt10_recolored",
            leak_emb_endo = 0.12,
            leak_exe_endo = 0.17)

write.table(as.data.frame(mc_leak), sep="\t", quote=F, file="scrna_db/leak.sing_emb_wt10_recolored.txt")

mat_id = "sing_emb_wt10"
mc_id = "sing_emb_wt10_recolored"
mgraph_id = "sing_emb_wt10_logist"
net_id = "sing_emb_wt10_logist"

mat = scdb_mat(mat_id)
md = mat@cell_metadata
cell_time = md[,"age_group"]

mgraph = scdb_mgraph(mgraph_id)
names(cell_time) = rownames(md)

leak_fn = "scrna_db/leak.sing_emb_wt10_recolored.txt"

mcell_mctnet_from_mgraph(net_id, 
                                mgraph_id,
                                cell_time,
                                leak_fn,
                            t_exp = 1,
                                T_cost = 1e+5,
                                capacity_var_factor = NULL,
                                off_capacity_cost1 = 1,
                                off_capacity_cost2 = 1000,
                                k_norm_ext_cost = 1,
                                k_ext_norm_cost = 1,
                                k_ext_ext_cost = 1)

##now generate the flow itself

flow_id = "sing_emb_wt10_logist"
fig_dir = "figs/"

flow_tolerance = 0.01

network_color_ord = NULL
  
message("generate flows")

mcell_new_mctnetflow(flow_id, net_id, 
                            init_mincost = T, flow_tolerance=0.01)

message("solved network flow problem")

mcf = scdb_mctnetflow(flow_id)
  
#compute propagatation forward and background
mcf = mctnetflow_comp_propagation(mcf)
  
#adding back the object with the network and flows
scdb_add_mctnetflow(flow_id, mcf)
```
