# source("https://gist.githubusercontent.com/mtandon09/4a870bf4addbe46e784059bce0e5d8d6/raw/dc2927aa3e6a09b34a39f8346b5ebcfd41ce2a6d/install_R_dependencies.R")

fungin_initialize <- function(fungin_data_dir=file.path("data","fungin"), sourcedb="reactomeFI",
                              rds_file=NULL, fromScratch=F, genome="hg38", string_score_threshold=NULL) {
  require(igraph)
  
  sourcedb=tolower(sourcedb)
  if (! sourcedb %in% c("reactomefi","stringdb")) {
    stop("'sourcedb' must be one of 'rectomeFI' or 'stringDB'.")
  }
  
  fungin_data_dir <- file.path(fungin_data_dir,sourcedb)
  if (is.null(rds_file)) {
    rds_file <- ifelse(sourcedb=="stringdb", file.path(fungin_data_dir,paste0("fungin.",sourcedb,".",genome,".Rds")),
                                               file.path(fungin_data_dir,paste0("fungin.",sourcedb,".Rds")))
  }
  if (!file.exists(rds_file) | fromScratch) {
    message(paste0("Building network data from scratch (",sourcedb,")...\n"))
    if (sourcedb=="reactomefi") {
      rds_file <- gather_FI_network_data(fungin_data_dir,rds_file=rds_file, fromScratch=fromScratch)
      
      
    } else if (sourcedb=="stringdb") {
      rds_file <- gather_stringdb_data(fungin_data_dir,rds_file=rds_file, fromScratch=fromScratch, genome=genome)
    }
  }
  load(rds_file)
  full_interaction_network <- delete.edges(full_interaction_network, E(full_interaction_network)[is.na(edge_attr(full_interaction_network,"combined_score"))])
  if (sourcedb=="stringdb" & !is.null(string_score_threshold)) {
    # browser()
    # get.edge.attribute(full_interaction_network,"combined_score")
    # full_interaction_network <- delete.edges(full_interaction_network, which(E(full_interaction_network)$combined_score >= string_score_threshold)-1)
    # string_score_threshold=800
    # load(rds_file)
    # dim(get.edgelist(full_interaction_network))
    # range(edge_attr(full_interaction_network,"combined_score"))
    # # full_interaction_network <- full_interaction_network %>% delete_edges(which(edge_attr(full_interaction_network,"combined_score") >= string_score_threshold)-1)
    full_interaction_network <- delete.edges(full_interaction_network, E(full_interaction_network)[!edge_attr(full_interaction_network,"combined_score") >= string_score_threshold])
    # dim(get.edgelist(full_interaction_network))
    # range(edge_attr(full_interaction_network,"combined_score"))
    # # isolated_nodes <- which(strength(full_interaction_network, mode="all") < 1)
    isolated_nodes <- which(degree(full_interaction_network) < 1)
    nodes_to_remove <- isolated_nodes[!is.na(isolated_nodes)]
    full_interaction_network <- simplify(delete_vertices(full_interaction_network, nodes_to_remove), edge.attr.comb="first")
  }
  full_interaction_network <- fungin_add_attr(full_interaction_network)
  
  
  return(full_interaction_network)
}


fungin_query <- function(fungin_graph, 
                         selection_method=c("both","deg","maf"),
                         deg_table=NULL, pval_cutoff=0.05, n_genes_cutoff=NULL,
                         gene_column="gene",fc_column="logFC",pval_column="adj.P.Val",
                         maf=NULL, min_mutated_samples = 1,
                         fillVar="mutated",
                         get_neighbors = F,
                         pathways=NULL,
                         max_pathways = 7,
                         savename=NULL) {
  
  
  # browser()
  query_genes <- fungin_gene_selector(selection_method=selection_method,
                                      deg_table=deg_table, pval_cutoff=pval_cutoff, n_genes_cutoff=n_genes_cutoff,
                                      gene_column=gene_column,fc_column=fc_column,pval_column=pval_column,
                                      maf=maf, min_mutated_samples = 1)
  if (length(query_genes) < 1) {
    message("No genes left after filtering with selection method.")
    return(invisible(NULL))
  }
  
  trimmed_network <- fungin_trim(fungin_graph, query_genes = query_genes, get_neighbors = get_neighbors)
  
  anno_network <- fungin_annotate(trimmed_network,
                                 diff_exp_results=deg_table, pval_cutoff=pval_cutoff, n_genes_cutoff=n_genes_cutoff,
                                 # gene_column="name",fc_column="logFC",pval_column="pval",
                                 gene_column=gene_column,fc_column=fc_column,pval_column=pval_column,
                                 maf=maf, min_mutated_samples = min_mutated_samples,
                                 pathways=pathways, max_pathways = max_pathways)
  
  myplot <- fungin_plot(anno_network, marker_size=40, fillVar=,savename=savename)
  
  return(invisible(myplot))
}

fungin_gene_selector <- function(selection_method=c("both","deg","maf"),
                                 deg_table=NULL, pval_cutoff=0.05, n_genes_cutoff=NULL,
                                 gene_column="gene",fc_column="logFC",pval_column="adj.P.Val",
                                 maf=NULL, min_mutated_samples = 1, pathogenic_only=F) {
  require(igraph)
  # all_vertex_attr <- data.frame(do.call(cbind,vertex_attr(fungin_graph)), stringsAsFactors = F)
  if (!is.null(deg_table)) {
    if (any(! c(gene_column, fc_column, pval_column) %in% colnames(deg_table))) {
      stop(paste0("Diff exp results are missing these columns: ", paste0(setdiff(c(gene_column, fc_column, pval_column), colnames(deg_table)), collapse=",")))
    }
    deg_table <- deg_table[,c(gene_column, fc_column, pval_column)]
    colnames(deg_table) <- c("name","logFC","pval")
  }
  
  if (is.null(n_genes_cutoff) || all(is.na(deg_table$pval))) {
    n_genes_cutoff <- nrow(deg_table)
  }
  deg_table <- deg_table[order(deg_table$pval, decreasing = F),]
  
  deggenes <- deg_table$name[deg_table$pval < pval_cutoff]
  deggenes <- deggenes[!is.na(deggenes)] 
  deggenes <- deggenes[1:min(c(length(deggenes), n_genes_cutoff))] 
  
  mafgenes <- c()
  if (! is.null(maf)) {
    require(maftools)
    mymaf=NULL
    if ("data.frame" %in% class(maf)) {
      if (nrow(maf) > 0) {
        mymaf <- maftools::read.maf(maf)
      }
    } else if (class(maf) == "MAF") {
      mymaf <- maf
    } 
    
    if (!is.null(mymaf)) {
      nsamples=nrow(mymaf@variants.per.sample)
      if (min_mutated_samples > 0 & min_mutated_samples < 1) {
        min_mutated_samples=min_mutated_samples*nsamples
      }
      mafgenes <- mymaf@gene.summary$Hugo_Symbol[mymaf@gene.summary$AlteredSamples >= min_mutated_samples]
      if (pathogenic_only) {
        if ("CLIN_SIG" %in% colnames(mymaf@data)) {
          mafgenes <- unique(mymaf@data$Hugo_Symbol[grepl("pathogenic",mymaf@data$CLIN_SIG)])
          message(paste0("Using ", length(mafgenes), " genes from MAF with pathogenic annotation..."))
        } else {
          warning("ClinVar annotations not found... ignoring pathogenic_only option")
        }
      }
    } else {
      warning("maf argument must be a MAF object or data.frame-like. No genes selected from maf.")
    }
  }
  
  selection_method <- tolower(selection_method)
  if ("both" %in% selection_method) {selection_method <- c("maf","deg")}
  query_genes <- c()
  if ("maf" %in% selection_method) {
    query_genes <- c(query_genes, mafgenes)
  }
  if ("deg" %in% selection_method) {
    query_genes <- c(query_genes, deggenes)
  }
  query_genes <- unique(query_genes)
  
  return(query_genes)
}

fungin_add_attr <- function(igraph_obj) {
  all_node_attrs <- names(edge_attr(igraph_obj))
  all_edge_attrs <- names(edge_attr(igraph_obj))
  if (! "logFC" %in% all_node_attrs) {
    igraph_obj <- set_vertex_attr(igraph_obj,name="logFC",value=NA)
  }
  if (! "pval" %in% all_node_attrs) {
    igraph_obj <- set_vertex_attr(igraph_obj,name="pval",value=NA)
  }
  if (! "pval_binary" %in% all_node_attrs) {
    igraph_obj <- set_vertex_attr(igraph_obj,name="pval_binary",value=NA)
  }
  if (! "mutated" %in% all_node_attrs) {
    igraph_obj <- set_vertex_attr(igraph_obj,name="mutated",value=NA)
  }
  if (! "strength" %in% all_node_attrs) {
    full_strengths <- strength(igraph_obj, mode="all")
    igraph_obj <- set_vertex_attr(igraph_obj,"strength",V(igraph_obj)[names(full_strengths)],full_strengths)
  }
  if (! "changetype" %in% all_edge_attrs) {
    igraph_obj <- set_edge_attr(igraph_obj,name="changetype",value="Other")
  }
  return(igraph_obj)
}

fungin_trim <- function(full_interaction_network, query_genes=NULL, get_neighbors=FALSE) {
  require(igraph)
  
  net_genes <- query_genes[query_genes %in% names(V(full_interaction_network))]
  if (length(net_genes)==0) {
    if (length(query_genes)>0){
      nquery=min(length(query_genes), 5)
      message(paste0("Query genes: ", paste0(query_genes[1:nquery], collapse=", "),",... (and ", length(query_genes)-nquery, " more)"))
    }
    stop("Query genes not found in interaction data")
    # stop(paste0("Top ", n_genes_cutoff, " not found in interaction data"))
    # return(NA)
  }
  
  # browser()
  
  query_nodes <- V(full_interaction_network)[net_genes]
  my_interaction_graph <- simplify(induced_subgraph(full_interaction_network,query_nodes), edge.attr.comb="first")
  
  if (get_neighbors=="auto"){
    get_neighbors=F
    if (length(E(my_interaction_graph)) < 1) {
      get_neighbors=T
    }
  }
  if (length(net_genes) > 100) {
    if (get_neighbors) {
      warning("Too many genes to get neighbors; setting 'get_neighbors' to FALSE...")
      get_neighbors=F
    }
  }
  
  if (get_neighbors) {
    print("Getting neighbors...")
    neighborhood <- ego(full_interaction_network, nodes=V(full_interaction_network)[net_genes])
    query_nodes <- unlist(neighborhood)
    my_interaction_graph <- simplify(induced_subgraph(full_interaction_network,query_nodes), edge.attr.comb="first")
  }
  
  # my_interaction_graph <- simplify(induced_subgraph(full_interaction_network,query_nodes))#, edge.attr.comb="sum")
  
  full_strengths <- strength(my_interaction_graph, mode="all")
  my_interaction_graph <- set_vertex_attr(my_interaction_graph,"strength",V(my_interaction_graph)[names(full_strengths)],full_strengths)
  
  isolated_nodes <- V(my_interaction_graph)$name[V(my_interaction_graph)$strength < 1]
  if (length(isolated_nodes)/length(V(my_interaction_graph)) > 0.8) {
    # stop(paste0("Not enough interactions left to plot.",ifelse(get_neighbors,""," Try setting 'get_neighbors' to TRUE.")))
    # warning(paste0("Not enough interactions to plot.",ifelse(get_neighbors,""," Try setting 'get_neighbors' to TRUE.")))
    # return(NA)
  }
  # nodes_to_remove <- union(nonsig_nodes, isolated_nodes)
  nodes_to_remove <- isolated_nodes[!is.na(isolated_nodes)]
  
  # browser()
  # curr_graph <- simplify(delete_vertices(curr_graph, nodes_to_remove), edge.attr.comb="sum")
  # my_interaction_graph <- simplify(delete_vertices(my_interaction_graph, nodes_to_remove))
  my_interaction_graph <- simplify(delete_vertices(my_interaction_graph, nodes_to_remove), edge.attr.comb="first")
  
  
  return(my_interaction_graph)
}


fungin_annotate <- function(fungin_graph, 
                            diff_exp_results=NULL, pval_cutoff=0.05, n_genes_cutoff=NULL,
                            gene_column="gene",fc_column="logFC",pval_column="adj.P.Val",
                            maf=NULL, min_mutated_samples = 1,
                            pathways=NULL,
                            max_pathways = 7,
                            fillVar="mutated") {
  require(igraph)
  # browser()
  all_vertex_attr <- data.frame(do.call(cbind,vertex_attr(fungin_graph)), stringsAsFactors = F)
  if (is.null(diff_exp_results)) {
    # diff_exp_results <- as.data.frame(matrix(0, nrow=1, ncol=3))
    # colnames(diff_exp_results) <- c(gene_column, fc_column, pval_column)
    diff_exp_results <- data.frame(name=names(V(fungin_graph)),
                                   all_vertex_attr[,c("logFC","pval")],
                                   stringsAsFactors = F)
  } else {
    # browser()
    if (any(! c(gene_column, fc_column, pval_column) %in% colnames(diff_exp_results))) {
      stop(paste0("Diff exp results are missing these columns: ", paste0(setdiff(c(gene_column, fc_column, pval_column), colnames(diff_exp_results)), collapse=",")))
    }
    diff_exp_results <- diff_exp_results[,c(gene_column, fc_column, pval_column)]
    colnames(diff_exp_results) <- c("name","logFC","pval")
  }
  
  diff_exp <- diff_exp_results[order(diff_exp_results$pval, decreasing = F),]
  if ( ! is.null(n_genes_cutoff) ) {
    diff_exp <- diff_exp[1:min(c(n_genes_cutoff, nrow(diff_exp))),]
  }
  mut_status <- data.frame(name=names(V(fungin_graph)),
                           mutated=all_vertex_attr[,c("mutated")],
                           stringsAsFactors = F)
  if (! is.null(maf)) {
    # browser()
    mut_summary <- make_mut_summary(maf)
    mut_summary$mutated <- ifelse(mut_summary$AlteredSamples >= min_mutated_samples, 
                                   paste0("Mutated in ", min_mutated_samples, " or more samples"),
                                   NA)
    mut_status <- data.frame(name=mut_summary$Hugo_Symbol,
                             mutated=mut_summary$mutated,
                             Nmutated=mut_summary$MutatedSamples,
                             stringsAsFactors = F)
  }
  
  curr_graph <- fungin_graph
  nodetypes <- rep("Gene",length(V(curr_graph)))
  if (! is.null(pathways) ){
    # browser()
    pathways <- data.frame(pathways)
    pathways <- pathways[pathways[,1] %in% names(V(curr_graph)),]
    # min_genes_for_pathway = 3
    # path_to_keep <- names(table(pathways[,2])[table(pathways[,2])>=min_genes_for_pathway])
    # max_pathways = 7
    path_to_keep <- names(sort(table(pathways[,2]), decreasing = T)[1:min(c(max_pathways, length(unique(pathways[,2]))))])
    pathways <- pathways[pathways[,2] %in% path_to_keep,]
    if (nrow(pathways) > 0) {
      pathgraph <- igraph::graph_from_data_frame(pathways, directed = is_directed(curr_graph))
      # pathgraph <- fungin_add_attr(pathgraph)
      curr_graph <- curr_graph %u% pathgraph
      curr_graph <- set_edge_attr(curr_graph, "changetype", which(is.na(get.edge.attribute(curr_graph,"changetype"))),value = "Pathway")
      nodetypes <- ifelse(names(V(curr_graph)) %in% pathways[,2], "Pathway","Gene")
      # curr_graph <- set_vertex_attr(curr_graph,"nodetype",V(curr_graph),nodetypes)
    }
    
  }
  # browser()
  my_attrs <- data.frame(vname=V(curr_graph)$name, stringsAsFactors = F)
  
  my_attrs$logFC <- diff_exp$logFC[match(my_attrs$vname, diff_exp$name)]
  my_attrs$pval <- diff_exp$pval[match(my_attrs$vname, diff_exp$name)]
  my_attrs$pval_binary <- ifelse(my_attrs$pval < pval_cutoff, paste0("p < ",pval_cutoff), 
                                 ifelse(my_attrs$pval >= pval_cutoff, "ns",NA))
  my_attrs$pval_binary[is.na(my_attrs$pval_binary)] <- "No data" 
  
  my_attrs$nodetype <- nodetypes
  my_attrs$mutated <- mut_status$mutated[match(my_attrs$vname, mut_status$name)]
  my_attrs$mutated[is.na(my_attrs$mutated)] <- "No data"
  my_attrs$mutated[my_attrs$nodetype=="Pathway"] <- "Pathway"
  my_attrs$Nmutated <- mut_status$Nmutated[match(my_attrs$vname, mut_status$name)]
  # print(table(my_attrs$mutated))
  # browser()
  curr_graph <- set_vertex_attr(curr_graph,"logFC",V(curr_graph)[my_attrs$vname],my_attrs$logFC)
  curr_graph <- set_vertex_attr(curr_graph,"pval",V(curr_graph)[my_attrs$vname],my_attrs$pval)
  curr_graph <- set_vertex_attr(curr_graph,"pval_binary",V(curr_graph)[my_attrs$vname],my_attrs$pval_binary)
  curr_graph <- set_vertex_attr(curr_graph,"mutated",V(curr_graph)[my_attrs$vname],my_attrs$mutated)
  curr_graph <- set_vertex_attr(curr_graph,"Nmutated",V(curr_graph)[my_attrs$vname],my_attrs$Nmutated)
  curr_graph <- set_vertex_attr(curr_graph,"nodetype",V(curr_graph)[my_attrs$vname],my_attrs$nodetype)
  
  # nonsig_nodes <- V(curr_graph)$name[V(curr_graph)$pval > pval_cutoff]
  
  
  return(curr_graph)
}





fungin_plot <- function(fungin_graph, marker_size=40, fillVar="logFC",savename=NULL) {

  
  # browser()
  require(ggnetwork)
  require(ggnewscale)
  layout_area_param=NULL
  plotdata <- ggnetwork(fungin_graph,scale = F, stringsAsFactors=F)
  if (nrow(plotdata) < 1 ){
    warning("Nothing to plot")
    return(NULL)
  }
  edgesize <- 1.5
  # shape_vals <- c(21, 22, 23, 24, 25)
  # names(shape_vals) <- sort(unique(plotdata$mutated))
  shape_vals <- c(21, 23)
  names(shape_vals) <- c("No data",grep("Mutated",unique(plotdata$mutated), value=T))
  
  interaction_colors <- c("Other"="grey70", "Inhibition"="blue2","Activation"="gold","Pathway"="wheat2")
  
  # alpha_vals <- c(1, 0.25, 0.1)
  # names(alpha_vals) <- rev(sort(unique(plotdata$pval_binary)))
  lty_vals <- c(1, 2, 3)
  names(lty_vals) <- rev(sort(unique(plotdata$pval_binary)))
  
  outline_vals <- c("grey30","grey60", "grey90")
  names(outline_vals) <- names(lty_vals)
  
  nedges <- length(E(fungin_graph))
  edgesize <- ifelse(nedges > 1000, 0.25, 
                     ifelse(nedges > 100, 0.5,
                            ifelse(nedges > 10, 2, edgesize)))
  
  nnodes <- length(V(fungin_graph))
  nodesize <- ifelse(nnodes > 1000, 2, 
                     ifelse(nnodes > 100, 5,
                            ifelse(nnodes > 10, 10, nnodes)))
  
  if (!fillVar %in% colnames(plotdata)) {
    warning("fillVar not found, using logFC for fill colors...")
    fillVar="logFC"
  }
  plotdata$fillvar <- as.numeric(plotdata[,fillVar])
  mypal <- "PiYG"
  # browser()
  maxFillValue<-max(abs(plotdata$fillvar), na.rm = T)
  color_limits <- c(-maxFillValue, maxFillValue)
  showFClegend=T
  if (all(is.na(plotdata$fillvar))) {
    plotdata$logFC <- 0
    showFClegend=F
    # suppressWarnings(maxFC<-max(abs(vis.nodes$logFC), na.rm = T))
  } else {
    if (min(plotdata$fillvar, na.rm=T) >= 0) {
      mypal <- "Reds"
      # color_limits <- range(plotdata$fillvar, na.rm = T)
      color_limits <- c(0, max(plotdata$fillvar, na.rm = T))
    }
  }
  
  
  # # gene_name_colors <- c("#ff3300","#cc0000","grey50")
  # # gene_name_colors <- c("#ff3300","#cc0000","grey50")
  # gene_name_colors <- c("#66dbff","black")
  # # names(gene_name_colors) <- sort(unique(plotdata$changetype))
  # names(gene_name_colors) <- c("Pathway","Gene")
  # # gene_name_colors["No data"] <- "black"
  
  # plotdata$mutated[plotdata$nodetype=="Pathway"] <- "Pathway"
  gene_name_colors <- c("grey20","wheat4","gold")
  names(gene_name_colors) <- c("No data",grep("Mutated","Pathway",unique(plotdata$mutated), value=T))
  
  gene_name_colors <- gene_name_colors[names(gene_name_colors) %in% unique(plotdata$mutated)]
  
  # browser()
  
  edgeplot <- ggplot(plotdata, 
                   # layout = "fruchtermanreingold", cell.jitter = 2, niter=1000, area=layout_area_param,
                   # layout = "mds",
                   # layout = "spring", repulse=T, mass=0.1, k=0.001, kfr=0.01, repeqdis=0.5,
                   layout = "kamadakawai", niter=10000,initemp=1000,cool.exp=0.1, #weights="strength",
                   # layout = "mds", niter=10, var="geodist", dist="maximum",
                   # layout = "hall", niter=100,
                   aes(x = x, y = y, xend = xend, yend = yend))+#,
    # layout = "target", niter=1000) +
    # geom_edges(color= "grey50", alpha=0.5,size=0.1,curvature = 0.2,
    geom_edges(data = subset(plotdata, changetype %in% "Other"),
               aes(color= changetype), alpha=0.3,size=edgesize,curvature = 0.2,
               arrow = arrow(length = unit(0, "pt"), type = "closed", angle=90)) +
    geom_edges(data = subset(plotdata, changetype %in% "Activation"),
               aes(color= changetype), alpha=0.5,size=edgesize,curvature = 0.2,
               arrow = arrow(length = unit(3, "pt"), type = "closed")) +
    geom_edges(data = subset(plotdata, changetype %in% "Inhibition"),
               aes(color= changetype), alpha=0.5,size=edgesize,curvature = 0.2,
               arrow = arrow(length = unit(3, "pt"), type = "closed", angle=90)) +
    geom_edges(data = subset(plotdata, changetype %in% "Pathway"), lty=2,
               aes(color= changetype), alpha=0.5,size=edgesize*0.5,curvature = 0.2,
               arrow = arrow(length = unit(0, "pt"), type = "closed", angle=90)) +
    scale_color_manual(values=interaction_colors)+
    labs(color = "Gene Interaction") 
  
  nodeplot <- edgeplot +
    new_scale_color()+
    # geom_nodes(aes(fill=logFC, alpha=pval_binary),
    # geom_nodes(data = subset(plotdata, mutated %in% "Inhibition"),
    # geom_nodes(aes(fill=fillvar, shape=mutated, alpha=pval_binary, color=pval_binary,size = strength),
    geom_nodes(aes(fill=fillvar, shape=mutated, color=pval_binary,size = strength),
               show.legend = showFClegend) +
    # size = nodesize, show.legend = showFClegend) +
    # color = "grey50", size = nodesize, lwd=2) +
    # shape = 21, color = "grey50") +#, size = nodesize, lwd=2) +
    scale_color_manual(values=outline_vals, guide="none")+
    labs(color = "logFC pval\n(outline color)",
         fill = fillVar) +
    scale_shape_manual(values=shape_vals) +
    labs(shape = "Mutations") +
    scale_fill_distiller(palette = mypal,na.value = "grey80", limit=color_limits, n.breaks=5, direction = 1) #+
    # scale_fill_gradientn(colors=fc_colors,values=fc_colors.values,na.value = "grey80") +
    # scale_alpha_manual(values=alpha_vals, guide="none")
    # scale_linetype_manual(values=lty_vals)
  # nodeplot
  annotatedplot <- nodeplot +
    new_scale_color()+
    geom_nodetext_repel(aes( label = name, color=mutated ),
                        # geom_nodetext(aes( label = vertex.names, color=logFC ),
                        # fontface = "bold",size=1, color="steelblue") +
                        # fontface = "bold",size=1, show.legend = F) +
                        fontface = "bold",#size=nodesize/3, 
                        show.legend = T) +
    # scale_colour_gradient2(low="grey70",mid="black",high="grey70") +
    scale_colour_manual(values=gene_name_colors) +
    labs(color = "Label Color") +
    ggtitle(paste0(length(V(fungin_graph))," genes plotted")) +
    theme_blank()
    
  # annotatedplot
  if (! is.null(savename)) {
    plotheight=scales::rescale(sqrt(nnodes), to=c(8,36), from=c(1,100))
    plotheight=min(c(plotheight, 40))
    ggsave(annotatedplot, filename = savename, width=plotheight, height=plotheight) 
  } #else {
  # print(myplot)
  # }
  
  
  return(annotatedplot)
  
}
fungin_plot_interactive <- function(fungin_graph,fill_var="mutated",marker_size=40) {
  # browser()
  plotdata <- fungin_make_plotdata(fungin_graph)
  nodedata <- plotdata$nodes
  edgedata <- plotdata$edges
  require(plotly)
  
  node_plot <- plot_ly(data=nodedata, 
                       x = ~Xn, y = ~Yn, 
                       mode = "markers", 
                       text = nodedata$name, hoverinfo = "text",
                       marker = list(
                         color = nodedata$color.background,
                         shape = nodedata$shape,
                         size = marker_size,
                         opacity=0.5,
                         line = list(
                           # color = 'red2',
                           color = nodedata$color.background,
                           width = 5
                         )#,
                         # selected = list(maker=list(opacity=0.1, color="red",size=20))
                       )
  )
  
  
  
  network <- node_plot %>% 
    add_annotations(
      data=edgedata,
      x = ~startx,
      y = ~starty,
      ax = ~endx,
      ay = ~endy,
      text = "",
      showarrow = TRUE,
      arrowcolor = ~color,
      arrowhead = ~arrowType,
      arrowsize = 1,
      arrowwidth = ~width,
      standoff = marker_size*0.9,
      startstandoff = marker_size*0.9,
      xref = "x",
      yref = "y",
      axref="x",
      ayref="y"
    )
  
  a <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
  mynetwork <- network %>% layout(xaxis = a, yaxis = a, dragmode="select",hovermode="closest")

  return(mynetwork)
  
}

fungin_make_plotdata <- function(fungin_graph) {
  require(igraph)
  L <- graph_attr(fungin_graph %>%
                    add_layout_(nicely(), component_wise()),
                  "layout")
  
  
  vis.nodes <- igraph::as_data_frame(fungin_graph,what = "vertices")
  vis.nodes$label  <- vis.nodes$name
  vis.nodes$Xn <- L[,1]
  vis.nodes$Yn <- L[,2]
  vis.nodes$color <- "lightblue"
  
  
  vis.nodes$shape <- ifelse(vis.nodes$mutated=="No data","circle",
                            ifelse(grepl("Mutated",vis.nodes$mutated), "star","square"))
  
  vis.nodes$color.background <- "grey90"
  # browser()
  maxFC<-NA
  if (!all(is.na(vis.nodes$logFC))) {
    # suppressWarnings(maxFC<-max(abs(vis.nodes$logFC), na.rm = T))
    maxFC<-max(abs(vis.nodes$logFC), na.rm = T)
  }
  if (!is.na(maxFC)) {
    require(circlize)
    fc_color_breaks <- seq(-maxFC, maxFC, length.out=6)
    # maxFC=4
    fc_color_breaks <- seq(-maxFC, maxFC, length.out=6)
    color_func <- colorRamp2(fc_color_breaks, brewer.pal(name = "PiYG",n = 6))
    vis.nodes$color.background <- color_func(vis.nodes$logFC)
  }
  vis.nodes$color.border <- "grey50"
  # vis.nodes$color.border[vis.nodes$id %in% query_genes] <- "red"
  # browser()
  
  vis.links <- igraph::as_data_frame(fungin_graph,what = "edges")
  # vis.links <- merge.data.frame(vis.links, gene_interactions, by.x=c("from","to"), by.y=c("Gene1","Gene2"), all.x=T)
  # vis.links$changetype <- ifelse(vis.links$change == 0, "Other",ifelse(vis.links$change < 0, "Inhibition","Activation"))
  # vis.links$width <- 4
  # vis.links$arrows <- "middle"
  
  interaction_colors <- c("Other"="grey70", "Inhibition"="blue2","Activation"="gold")
  # vis.links$color.background <- interaction_colors[vis.links$changetype]
  vis.links$color <- interaction_colors[vis.links$changetype]
  
  vis.links$width <- ifelse(vis.links$changetype == "Other", 1,3)
  vis.links$arrowType <- ifelse(vis.links$changetype == "Other",0,ifelse(vis.links$changetype=="Inhibits",7,2))
  # browser()
  # marker_offset <- marker_size/100
  vis.links$startx <- vis.nodes$Xn[match(vis.links$from, vis.nodes$name)]#+marker_offset
  vis.links$starty <- vis.nodes$Yn[match(vis.links$from, vis.nodes$name)]#+marker_offset
  vis.links$endx <- vis.nodes$Xn[match(vis.links$to, vis.nodes$name)]#-marker_offset
  vis.links$endy <- vis.nodes$Yn[match(vis.links$to, vis.nodes$name)]#-marker_offset
  
  
  return(list(nodes=vis.nodes, edges=vis.links))
}

make_mut_summary <- function(maf) {
  require(maftools)
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
  
  # match_idx <- match(V(my_interaction_graph)$name,gene_mut_summary$Hugo_Symbol,nomatch = 0)
  # mut_status <- gene_mut_summary[match_idx,c("Hugo_Symbol",	"AlteredSamples",	"AlteredFraction")]
  # gene_mut_summary$mutated <- ifelse(gene_mut_summary$AlteredSamples >= min_mutated_samples, 
  #                                 paste0("Mutated in ", min_mutated_samples, " or more samples"),
  #                                 NA)

  # return_df <- data.frame(name=gene_mut_summary$Hugo_Symbol, mutated=gene_mut_summary$mutated, stringsAsFactors = F)
  return(gene_mut_summary)
}


gather_FI_network_data <- function(dataDir=file.path("data","fungin"), rds_file=NULL, fromScratch=F) {
  require(igraph)
  if (is.null(rds_file)) {
    rds_file <- file.path(dataDir, "fungin.reactomefi.Rds")
  }
  
  if (!file.exists(rds_file) | fromScratch ) {
    require(openxlsx)
    if (! dir.exists(dataDir)) {
      dir.create(dataDir, recursive = T)
    }
    
    print(paste0("Getting FI data..."))
    reactome_FI_url="https://reactome.org/download/tools/ReatomeFIs/FIsInGene_020720_with_annotations.txt.zip"
    reactome_FI_file <- file.path(dataDir,basename(reactome_FI_url))
    if (!file.exists(reactome_FI_file)) {
      download.file(reactome_FI_url,reactome_FI_file)
    }
    # browser()
    print(paste0("Reading FI data..."))
    reactome_FIs <- read.table(unz(reactome_FI_file, gsub(".zip","",basename(reactome_FI_file))), sep="\t", header=T, stringsAsFactors = F)
    reactome_FIs$predicted <- grepl("predicted",reactome_FIs$Annotation)
    
    change_factor <- list(
      "|-|"= c(-1,-1),
      "<-|" = c( 1,-1),
      "|->" = c(-1, 1),
      "-|"  = c( 0,-1),
      "|-"  = c(-1, 0),
      "<->" = c( 1, 1),
      "->"  = c( 0, 1),
      "<-"  = c( 1, 0),
      "-"   = c( 0, 0)
    )
    
    include_columns=c("Annotation","predicted")
    
    print(paste0("Processing FI data (this can take 2-3 minutes)..."))
    ## Here we expand all the interactions so that interaction direction is from Column 1 to Column 2
    ## This takes a minute or two
    ## Gene interaction list
    gene_int_list <- apply(reactome_FIs, 1, function(curr_data){
      genes=curr_data[c("Gene1","Gene2")]
      change_vec=change_factor[[curr_data["Direction"]]]
      
      expanded_df <- cbind(
        do.call(rbind,lapply(change_vec, function(x){
          return_val=genes
          if(x<0){return_val=rev(genes)}
          return(return_val)
        })),
        change=change_vec,
        matrix(rep(curr_data[include_columns],each=2), 
               nrow = length(change_vec), 
               ncol=length(include_columns),
               dimnames = list(c(1:length(change_vec)), include_columns))
      )
      
      return(data.frame(expanded_df, stringsAsFactors = F))
    })
    
    ## Bind into data frame
    gene_int_mat <- do.call(rbind, gene_int_list)
    gene_interactions <- unique(gene_int_mat)
    
    
    
    # browser()
    
    print(paste0("Geting miRNA target data..."))
    ## Add miRNA target data
    mirna_data_url="https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2019/cache/download/8.0/miRTarBase_MTI.xlsx"
    mirna_data_file=file.path(dataDir, basename(mirna_data_url))
    if (!file.exists(mirna_data_file)) {
      options(timeout = max(300, getOption("timeout")))
      download.file(url=mirna_data_url,destfile = mirna_data_file)
    }
    
    print(paste0("Reading miRNA target data..."))
    mirna_data.raw <- read.xlsx(mirna_data_file)
    mirna_data <- mirna_data.raw[mirna_data.raw$`Species.(Target.Gene)`=="Homo sapiens",]
    mirna_data <- mirna_data[mirna_data$Support.Type=="Functional MTI",]
    
    print(paste0("Processing miRNA target data..."))
    mirna_interactions <- data.frame(Gene1=tolower(mirna_data$miRNA), 
                                     Gene2=mirna_data$Target.Gene,
                                     change=-1,
                                     Annotation=mirna_data$Support.Type,
                                     stringsAsFactors=F)
    mirna_interactions <- unique(mirna_interactions)
    mirna_interactions$predicted <- F
    # mirna_interactions$Annotation <- mirna_interactions$Support.Type 
    mirna_interactions$data_source <- "mirTarBase"
    # browser()
    if (! "data_source" %in% colnames(gene_interactions)) {
      gene_interactions$data_source="Reactome_FI"
    }
    gene_interactions <- rbind(gene_interactions, mirna_interactions)
    
    print(paste0("Combining interactions data..."))
    gene_interactions$change <- as.numeric(gene_interactions$change)
    gene_interactions$changetype <- ifelse(gene_interactions$change == 0, "Other",ifelse(gene_interactions$change < 0, "Inhibition","Activation"))    
    
    full_interaction_network <- igraph::graph_from_data_frame(gene_interactions)
    
    print(paste0("Saving data..."))
    save(full_interaction_network, file = rds_file)
  } else {
    # load(rds_file)    
  }
  
  return(rds_file)
}




gather_stringdb_data <- function(dataDir=file.path("data","fungin"), rds_file=NULL, 
                                 fromScratch=F, genome="hg38",
                                 stringdb_version="11") {
  
  species_taxids <- c(hg19=9606, hg38=9606, mm10=10090)
  if (! genome %in% names(species_taxids)) {
    stop(paste0("'genome' must be one of '", paste0(names(species_taxids), collapse=","),"'"))
  }
  taxid=species_taxids[genome]
  if (is.null(rds_file)) {
    rds_file <- file.path(dataDir, paste0("fungin.stringdb.",taxid,".Rds"))
  }
  
  if (!file.exists(rds_file) | fromScratch ) {
    require(STRINGdb)
    if (! dir.exists(dataDir)) {
      dir.create(dataDir, recursive = T)
    }
    
    print(paste0("Getting STRINGdb data..."))
    # browser()
    # mygenomeinfo <- detect_maf_genome(mafobj_strict)
    options(timeout = max(300, getOption("timeout")))
    string_db <- STRINGdb$new( version=stringdb_version, species=taxid, score_threshold=0, input_directory=dataDir )# ,
    # backgroundV = unique(c(mafobj_strict@maf.silent$Hugo_Symbol, mafobj_strict@data$Hugo_Symbol)))
    string_ppi_graph <- string_db$get_graph()
    
    require(org.Hs.eg.db)
    allgenes <- unique(keys(org.Hs.eg.db,keytype = "SYMBOL"))
    string_id_mapfile <- file.path(dataDir,"stringdb_id_map.txt")
    if (!file.exists(string_id_mapfile)) {
      string_id_map <- string_db$map(data.frame(gene=allgenes, stringsAsFactors = F), "gene")
      string_id_map <- string_id_map[!is.na(string_id_map$STRING_id),]
      write.table(string_id_map, file=string_id_mapfile, sep="\t", row.names=F)
    } else {
      string_id_map <- read.table(string_id_mapfile, sep="\t", header = T, stringsAsFactors = F)
    }
    
    nodenames <- names(V(string_ppi_graph))
    id2symbol <- string_id_map$gene[match(nodenames, string_id_map$STRING_id)]
    id2symbol <- ifelse(is.na(id2symbol), gsub(paste0(taxid,"\\."),"",names(V(string_ppi_graph))), id2symbol)
    full_interaction_network <- set.vertex.attribute(string_ppi_graph, "name", value=id2symbol)
    save(full_interaction_network, file = rds_file)
  } else {
    # load(rds_file)
  }
  return(rds_file)
}
