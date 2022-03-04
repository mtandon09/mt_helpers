
add_category_to_dotplot <- function(reactome_enrich_res, show_n_path = 20, id_conv_table=NULL){
#### Takes a reactomePA enrichment object and adds Reactome categories (top level of hierarchy)
#### Returns a data frame with the added info, a plot of the results, and the modified enrichment object
  
  if (class(reactome_enrich_res) %in% c("compareClusterResult")) {
    ck.reactome.mod <- reactome_enrich_res   # The enrichment results
    mydf <- attr(ck.reactome.mod, "compareClusterResult")  # This extracts the data.frame storing the info in the enrichment object
  } else if (class(reactome_enrich_res) %in% c("data.frame")) {
    required_columns <- c("ID", "geneID","Description","contrast")
    mydf <- reactome_enrich_res
    stopifnot(all(required_columns %in% colnames(mydf)))
  }
  
  reactome_data <- load_reactome_data()    # Function to load (and create) reactome annotation data
  ## Add category annotation column
  mydf$category <- as.character(reactome_data$categoryName[match(mydf$ID, reactome_data$reactomeID)])
  
  ## A few things to pretty it up for plotting later
  mydf$category[is.na(mydf$category)] <- "Unclassified"
  mydf$category <- factor(mydf$category, levels = names(sort(table(mydf$category),decreasing = TRUE)))
  
  if (!is.null(id_conv_table)) {
    print("Converting gene symbols")
    mydf$geneSymbol <- sapply(strsplit(mydf$geneID, "/"), function(x) {
      paste0(as.character(id_conv_table$SYMBOL[match(x, id_conv_table$ENTREZ)]),collapse="/")
    })
  }
  require(ggplot2)
  plotdata <- data.frame(mydf, stringsAsFactors = F)
  ## This prioritizes the pathways that occur most often, both within each category and across contrasts
  pathway_counts <- sort(rowSums(table(plotdata[,c("Description","contrast")])), decreasing = T)
  pathways_to_show <- names(pathway_counts)[1:show_n_path]
  plotdata <- mydf[mydf$Description %in% pathways_to_show,]
  plotdata$Description <- factor(plotdata$Description, levels=pathways_to_show)
  ## Set levels of factor with the sorted order because that's what ggplot uses to order stuff
  plotdata$category <- factor(plotdata$category, levels=names(sort(table(plotdata$category), decreasing = T)))
  plotdata$Description_full <- plotdata$Description 
  # browser()
  nmax=40
  plotdata$Description <- unlist(lapply(as.character(plotdata$Description),function(x){
              ifelse(nchar(x) > nmax,
                     paste0(substr(x,1,nmax),"..."),
                     x)
    })
  )
  plotdata$Description <- factor(plotdata$Description, 
                                 levels=unique(as.character(plotdata$Description[order(plotdata$category, plotdata$Description, plotdata$Count, 
                                                                                       decreasing = T)])))
  
  
  # mychar <- c("aaaa","bb","ccc","dddddddddddddd")
  
  
  ## Plot
  myplot <- ggplot(plotdata, aes(x=contrast,y=Description)) +
    facet_grid(~category, scales = "free_y", labeller = ggplot2::labeller(category = ggplot2::label_wrap_gen())) +
    geom_point(aes(shape=direction,fill=direction, size=Count),color="black") +
    scale_shape_manual(values=c(up=24, down=25, "NA"=21)) +
    ggtitle(paste0("Pathway enrichment among DE genes")) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  
  ## Add modified data frame back into enrichment object
  attr(ck.reactome.mod, "compareClusterResult") <- mydf
  
  ## Return the data frame, the plot, and the modified enrichment object
  return(list(mydf,myplot,ck.reactome.mod))
}

make_df_for_clusterprofiler <- function(limma_res, pval_thresh = 0.05, pval_col="adj.P.Val", symbol_col="gene",fc_col="logFC",unlog_fc=NULL) {
#### Takes data frame of fold changes, p values, and gene symbols (like toptable output from a limma fit)
#### and makes a data frame suitable for input into the enricher() function in clusterProfiler
  
  # limma_res <- all_res[[1]]
  # pval_thresh=0.05
  # pval_col="adj.P.Val"
  # symbol_col="gene"
  # fc_col="logFC"
  require(clusterProfiler)
  mysym <- unique(as.character(limma_res[limma_res[,pval_col] < pval_thresh, symbol_col]))
  ids_converted <- bitr(mysym,OrgDb = "org.Hs.eg.db",fromType = "SYMBOL",toType="ENTREZID")
  entrez <- as.character(ids_converted[, "ENTREZID"])
  
  fc_vals <- limma_res[match(ids_converted$SYMBOL, limma_res[,symbol_col]),fc_col]
  if (is.null(unlog_fc)) {
    if (grepl("log",fc_col)) {
      unlog_fc=T
    } else {
      unlog_fc = F
    }
  }
  
  if (unlog_fc) {
    fc_vals <- ifelse(2^fc_vals > 1, 2^fc_vals, -1/(2^fc_vals))
  }
  df_for_clusterProfiler <- data.frame(symbol = ids_converted$SYMBOL,
                                       entrez = entrez, 
                                       pval=limma_res[match(ids_converted$SYMBOL, limma_res[,symbol_col]),pval_col],
                                       fc=fc_vals)
  # browser()
  df_for_clusterProfiler$direction <- factor(ifelse(df_for_clusterProfiler$fc > 0,"up","down"),levels=c("up","down","NA"))
  
  return(list(enrichinput=df_for_clusterProfiler, symbol_converter=ids_converted))
}


find_top_level <- function(my_id, conversion_table, child_col=2, parent_col = 1) {
#### This will recursively look through the hierarchy table to find the top most level
####  Multiple top levels are possible and returned as a vector
####  Works with the hierarchy info provided by Reactome
  
  id_matches <- which(conversion_table[,child_col] %in% my_id)
  # browser()
  if (length(id_matches) > 0) {
    next_ids <- conversion_table[id_matches, parent_col]
    myparent <- find_top_level(next_ids,conversion_table)
  } else {
    myparent_id <- which(conversion_table[, parent_col] %in% my_id)
    myparent <- unique(conversion_table[myparent_id, parent_col])
  }
  return(myparent)
}

load_reactome_data <- function(flat_file = "data/reactome.MT.txt", expand_multi_parents=FALSE) {
#### Downloads pathway and hierarchy info from Reactome
#### Contructs and saves pathway info with top level pathway as the "category"
####  Will load saved file if available
  
  if (file.exists(flat_file)) {
    reactome <- read.table(flat_file, sep="\t", header = T, quote = "")
  } else {
    print("Creating Reactome category information from scratch... ")
    pathways_file="data/ReactomePathways.txt"
    hierarchy_file="data/ReactomePathwaysRelation.txt"
    my_species="Homo sapiens"
    
    if (!file.exists(pathways_file)) {
      if (!dir.exists(dirname(pathways_file))) { dir.create(dirname(pathways_file),recursive = T) }
      download.file(url="https://reactome.org/download/current/ReactomePathways.txt",destfile = pathways_file)
    }
    if (!file.exists(hierarchy_file)) {
      if (!dir.exists(dirname(hierarchy_file))) { dir.create(dirname(hierarchy_file),recursive = T) }
      download.file(url="https://reactome.org/download/current/ReactomePathwaysRelation.txt",destfile = hierarchy_file)
    }
    
    pathway_info <- read.table(pathways_file,sep="\t",quote="",header = F, stringsAsFactors = F)
    colnames(pathway_info) <- c("reactomeID","pathwayName","species")
    pathway_info <- pathway_info[pathway_info$species==my_species,]
    
    hierarchy_info <- read.table(hierarchy_file,sep="\t",quote="",header = F, stringsAsFactors = F)
    colnames(hierarchy_info) <- c("parent","child")
    hierarchy_info <- hierarchy_info[hierarchy_info$child %in% pathway_info$reactomeID,]
    
    all_parents <- lapply(pathway_info$reactomeID,find_top_level,hierarchy_info)
    
    ### Looks like 7 pathways have 2 top level pathways, so paste everything as a comma-separated list
    pathway_info$categoryID <- unlist(lapply(all_parents,paste0,collapse=","))
    
    all_parents.names <- lapply(all_parents, function(curr_id, path_info) {
      path_info[match(curr_id, path_info[,"reactomeID"]),"pathwayName"]
    }, pathway_info)
    pathway_info$categoryName <- unlist(lapply(all_parents.names,paste0,collapse=","))
    
    reactome <- pathway_info
    # save(reactome, file = "data/reactome.db.mt.Rdata")
    write.table(reactome, file = "data/reactome.MT.txt", sep="\t",row.names = F,quote=F)
  }
  
  if (expand_multi_parents) {
    reactome.expand <- do.call(rbind,apply(reactome,1,function(curr_row) {
      cats <- unlist(strsplit(as.character(curr_row["categoryID"]),","))
      
      data.frame(reactomeID=rep(curr_row["reactomeID"],length(cats)),
                 pathwayName=rep(curr_row["pathwayName"],length(cats)),
                 categoryID=cats)
    }))
    reactome.expand$categoryName <- reactome.expand$pathwayName[match(reactome.expand$categoryID, reactome.expand$reactomeID)]
    output_df <- reactome.expand
  } else {
    output_df <- reactome
  }
  
  return(output_df)
}
  
  