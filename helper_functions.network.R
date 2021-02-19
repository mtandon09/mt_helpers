plot_interaction_network <- function(diff_exp_results=NULL, query_genes=NULL, pval_cutoff=0.05, n_genes_cutoff=10,
                             get_neighbors=FALSE,
                             gene_column="gene",fc_column="logFC",pval_column="adj.P.Val",
                             maf=NULL, min_mutated_samples = 1,
                             savename=NULL,
                             gene_interaction_saved_data="data/gene_interactions.Rdata") {
  
  require(igraph)
  require(qgraph)
  require(ggnetwork)
  require(network)
  require(intergraph)
  require(RColorBrewer)
  require(ggnewscale)
  require(maftools)
  
  if (file.exists(gene_interaction_saved_data)) {
    load(gene_interaction_saved_data)
  } else {
    stop("Need saved network interaction data.")
  }
  
  full_interaction_network <- graph_from_data_frame(gene_int_mat.uniq)
  
  if (is.null(diff_exp_results)) {
    diff_exp_results <- as.data.frame(matrix(0, nrow=1, ncol=3))
    colnames(diff_exp_results) <- c(gene_column, fc_column, pval_column)
  }
  diff_exp_results <- diff_exp_results[order(diff_exp_results[,pval_column], decreasing = F),]
  
  # browser()
  if (is.null(query_genes)) {
    query_genes <- diff_exp_results[diff_exp_results[,pval_column] < pval_cutoff,gene_column]
    query_genes <- query_genes[1:n_genes_cutoff]
    title_text <- paste0("top ",length(query_genes), " genes with p-value < ", pval_cutoff)
  } else {
    title_text <- paste0(length(query_genes), " specified genes")
  }
  title_text <- paste0(title_text, ifelse(get_neighbors, "\nNeighbors included",""))
  
  query_genes <- query_genes[query_genes %in% names(V(full_interaction_network))]
  if (length(query_genes)==0) {
    warning(paste0("Top ", n_genes_cutoff, " not found in interaction data"))
    # stop(paste0("Top ", n_genes_cutoff, " not found in interaction data"))
    return(NA)
  }
  
  if (get_neighbors) {
    neighborhood <- ego(full_interaction_network, nodes=V(full_interaction_network)[query_genes])
    query_nodes <- unlist(neighborhood)
  } else {
    query_nodes <- V(full_interaction_network)[query_genes]
  }
  
  my_interaction_graph <- simplify(induced_subgraph(full_interaction_network,query_nodes), edge.attr.comb="sum")
  match_idx <- match(V(my_interaction_graph)$name,diff_exp_results[,gene_column],nomatch = 0)
  diff_exp <- diff_exp_results[match_idx,c(gene_column, fc_column, pval_column)]
  
  if (! is.null(maf)) {
    if (class(maf)=="character") {
      if (file.exists(maf)) {
        mymaf <- read.maf(maf)
      } else {
        stop(paste0("MAF file does not exist: ",maf))
      }
    } else if (class(maf)=="MAF") {
      mymaf <- maf
    } else {
      stop("Argument 'maf' must be a path to a MAF file or a maftools object.")
    }
    # browser()
    gene_mut_summary <- mymaf@gene.summary
    # flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
    # gene_mut_summary <- gene_mut_summary[!gene_mut_summary$Hugo_Symbol %in% flag_genes,]
    num_samples <- as.integer(mymaf@summary$summary[mymaf@summary$ID=="Samples"])
    gene_mut_summary$AlteredFraction <- gene_mut_summary$AlteredSamples/num_samples
    
    match_idx <- match(V(my_interaction_graph)$name,gene_mut_summary$Hugo_Symbol,nomatch = 0)
    mut_status <- gene_mut_summary[match_idx,c("Hugo_Symbol",	"AlteredSamples",	"AlteredFraction")]
    mut_status$is_altered <- ifelse(mut_status$AlteredSamples >= min_mutated_samples, 
                                    paste0("Mutated in ", min_mutated_samples, " or more samples"),
                                    NA)
  } else {
    mut_status <- data.frame(Hugo_Symbol=V(my_interaction_graph)$name)
    mut_status$is_altered <- NA
  }
  
  curr_graph <- my_interaction_graph
  my_attrs <- data.frame(vname=V(curr_graph)$name)
  
  my_attrs$logFC <- diff_exp[match(my_attrs$vname, diff_exp$gene), fc_column]
  my_attrs$pval <- diff_exp[match(my_attrs$vname, diff_exp$gene), pval_column]
  my_attrs$pval_binary <- ifelse(my_attrs$pval < pval_cutoff, paste0("p < ",pval_cutoff), 
                                 ifelse(my_attrs$pval >= pval_cutoff, "ns",NA))
  my_attrs$pval_binary[is.na(my_attrs$pval_binary)] <- "No data" 
  
  my_attrs$mutated <- mut_status$is_altered[match(my_attrs$vname, mut_status$Hugo_Symbol)]
  my_attrs$mutated[is.na(my_attrs$mutated)] <- "No data"
  # print(table(my_attrs$mutated))
  # browser()
  curr_graph <- set_vertex_attr(curr_graph,"logFC",V(curr_graph)[my_attrs$vname],my_attrs$logFC)
  curr_graph <- set_vertex_attr(curr_graph,"pval",V(curr_graph)[my_attrs$vname],my_attrs$pval)
  curr_graph <- set_vertex_attr(curr_graph,"pval_binary",V(curr_graph)[my_attrs$vname],my_attrs$pval_binary)
  curr_graph <- set_vertex_attr(curr_graph,"mutated",V(curr_graph)[my_attrs$vname],my_attrs$mutated)
  
  nonsig_nodes <- V(curr_graph)$name[V(curr_graph)$pval > pval_cutoff]
  
  full_strengths <- strength(curr_graph, mode="all")
  curr_graph <- set_vertex_attr(curr_graph,"strength",V(curr_graph)[names(full_strengths)],full_strengths)
  
  # isolated_nodes <- V(curr_graph)$name[full_strengths < 1]
  # isolated_nodes <- V(curr_graph)$name[match(names(full_strengths)[full_strengths < 1], V(curr_graph)$name)]
  isolated_nodes <- V(curr_graph)$name[V(curr_graph)$strength < 1]
  
  if (length(isolated_nodes)/length(V(curr_graph)) > 0.8) {
    warning("Too few edges, so skipping...")
    # isolated_nodes <- c()
    return(NA)
  }
  nodes_to_remove <- union(nonsig_nodes, isolated_nodes)
  nodes_to_remove <- nodes_to_remove[!is.na(nodes_to_remove)]
  
  
  curr_graph <- simplify(delete_vertices(curr_graph, nodes_to_remove), edge.attr.comb="sum")
  
  strengths <- strength(curr_graph, mode="out")
  edgelist<-data.frame(get.edgelist(curr_graph))
  weights <- strengths[match(edgelist[,1], names(strengths))]
  # weights <- 10^(1-weights/max(weights))
  weights <- 2^(1-weights/max(weights))
  curr_graph <- set_edge_attr(curr_graph, "weights", value=weights)
  
  change_bin <- ifelse(edge_attr(curr_graph,"change")==0,"Other",ifelse(edge_attr(curr_graph,"change")<0,"Inhibition", "Activation"))
  curr_graph <- set_edge_attr(curr_graph,"interaction_type",value=change_bin)
  
  if (length(V(curr_graph)$name) < 1) {
    browser()
  }
  
  plotgraph <- asNetwork(curr_graph)
  
  shape_vals <- c(21, 22, 23, 24, 25)
  names(shape_vals) <- sort(unique(get.vertex.attribute(plotgraph, "mutated")))
  
  # outline_vals <- c("gold","grey40","grey90")
  # names(outline_vals) <- sort(unique(get.vertex.attribute(plotgraph, "mutated")))
  
  # interaction_colors <- c("Unknown"="grey70", "Inhibition"="blue2","Activation"="red3")
  interaction_colors <- c("Other"="grey70", "Inhibition"="blue2","Activation"="gold")
  
  
  alpha_vals <- c(1, 0.25, 0.1)
  names(alpha_vals) <- rev(sort(unique(get.vertex.attribute(plotgraph, "pval_binary"))))
  
  outline_vals <- c("grey30","grey80", "grey90")
  names(outline_vals) <- names(alpha_vals)
  
  # layout_area_param=length(network.vertex.names(plotgraph))^1
  layout_area_param=NULL
  plotdata <- ggnetwork(plotgraph,weights = "weights")
  edgesize <- 1.5
  nedges <- length(E(curr_graph))
  edgesize <- ifelse(nedges > 1000, 0.25, 
                     ifelse(nedges > 100, 0.5,
                            ifelse(nedges > 10, 1, edgesize)))
  
  nnodes <- length(V(curr_graph))
  nodesize <- ifelse(nnodes > 1000, 2, 
                     ifelse(nnodes > 100, 5,
                            ifelse(nnodes > 10, 8, nnodes)))
  
  suppressWarnings(fc_color_limit <- max(abs(get.vertex.attribute(plotgraph, "logFC")), na.rm = T) * c(-1, 1))
  # fc_colors <- colorRampPalette(rev(brewer.pal(11,"PiYG")))(100)
  # fc_colors.values <- seq(range(plotdata$logFC)[1],range(plotdata$logFC)[2],length.out = 100)
  # browser()
  # fc_range <- range(plotdata$logFC, na.rm=T)
  # fc_colors.values <- seq(max(abs(fc_range))*-1,max(abs(fc_range))*1,length.out = 100)
  if (sum(is.infinite(fc_color_limit)) > 0) {
    fc_color_limit <- c(0,0)
    showFClegend=F
  } else {
    fc_color_limit <- max(abs(get.vertex.attribute(plotgraph, "logFC")), na.rm = T) * c(-1, 1)
    showFClegend=T
  }
  
  # gene_name_colors <- c("#ff3300","#cc0000","grey50")
  gene_name_colors <- c("#ff3300","#cc0000","grey50")
  names(gene_name_colors) <- sort(unique(my_attrs$mutated))
  gene_name_colors["No data"] <- "black"
  
  # browser()
  
  myplot <- ggplot(plotdata, 
                   # layout = "fruchtermanreingold", cell.jitter = 2, niter=1000, area=layout_area_param,
                   # layout = "mds",
                   # layout = "spring", repulse=T, mass=0.1, k=0.001, kfr=0.01, repeqdis=0.5,
                   # layout = "kamadakawai", niter=10000,initemp=1000,cool.exp=0.1,
                   # layout = "mds", niter=10, var="geodist", dist="maximum",
                   layout = "hall", niter=100,
                   aes(x = x, y = y, xend = xend, yend = yend))+#,
                   # layout = "target", niter=1000) +
    # geom_edges(color= "grey50", alpha=0.5,size=0.1,curvature = 0.2,
    geom_edges(data = subset(plotdata, interaction_type %in% "Other"),
               aes(color= interaction_type), alpha=0.5,size=edgesize,curvature = 0.2,
               arrow = arrow(length = unit(0, "pt"), type = "closed", angle=90)) +
    geom_edges(data = subset(plotdata, interaction_type %in% "Activation"),
               aes(color= interaction_type), alpha=0.5,size=edgesize,curvature = 0.2,
               arrow = arrow(length = unit(3, "pt"), type = "closed")) +
    geom_edges(data = subset(plotdata, interaction_type %in% "Inhibition"),
               aes(color= interaction_type), alpha=0.5,size=edgesize,curvature = 0.2,
               arrow = arrow(length = unit(3, "pt"), type = "closed", angle=90)) +
    scale_color_manual(values=interaction_colors)+
    labs(color = "Gene Interaction") +
    new_scale_color()+
    # geom_nodes(aes(fill=logFC, alpha=pval_binary),
    # geom_nodes(data = subset(plotdata, mutated %in% "Inhibition"),
    geom_nodes(aes(fill=logFC, shape=mutated, alpha=pval_binary, color=pval_binary,size = strength),
               show.legend = showFClegend) +
               # size = nodesize, show.legend = showFClegend) +
               # color = "grey50", size = nodesize, lwd=2) +
               # shape = 21, color = "grey50") +#, size = nodesize, lwd=2) +
    scale_color_manual(values=outline_vals)+
    labs(color = "logFC pval\n(outline color)",
         fill = "Expression logFC") +
    scale_shape_manual(values=shape_vals) +
    labs(shape = "Mutations") +
    scale_fill_distiller(palette = "PiYG",na.value = "grey80", limit=fc_color_limit, direction = 1) +
    # scale_fill_gradientn(colors=fc_colors,values=fc_colors.values,na.value = "grey80") +
    scale_alpha_manual(values=alpha_vals, guide=FALSE)+
    new_scale_color()+
    # geom_nodetext(aes( label = vertex.names, color=mutated ),
    geom_nodetext_repel(aes( label = vertex.names, color=mutated ),
    # geom_nodetext(aes( label = vertex.names, color=logFC ),
                  # fontface = "bold",size=1, color="steelblue") +
                  # fontface = "bold",size=1, show.legend = F) +
                  fontface = "bold",size=nodesize/3, show.legend = T) +
    # scale_colour_gradient2(low="grey70",mid="black",high="grey70") +
    scale_colour_manual(values=gene_name_colors) +
    labs(color = "Mutations") +
    ggtitle(paste0("Interactions among ", title_text,"\n(",length(V(curr_graph))," genes plotted)")) +
    theme_blank()
  
  if (! is.null(savename)) {
    ggsave(myplot, filename = savename, width=8, height=8) 
  } #else {
    # print(myplot)
  # }
  
  
  return(myplot)
}



