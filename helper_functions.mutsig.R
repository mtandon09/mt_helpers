
load_cosmic_data <- function(signatures_file=file.path("~","Documents","helper_functions","cosmic_data","COSMIC_Mutational_Signatures_v3.1.xlsx"),
                             etio_data_xlsx=file.path("~","Documents","helper_functions","cosmic_data","COSMIC_signature_etiology.xlsx"),
                             signature_type="SBS") {
  
  require(openxlsx)
  # signatures_csv <- file.path("data","cosmic","sigProfiler_exome_SBS_signatures.csv")
  # cosmic_signatures = read.table(signatures_csv, sep = ",", header = TRUE, stringsAsFactors = F)
  if (! file.exists(signatures_file)) {stop("Can't find signatures file.")}
  cosmic_signatures = read.xlsx(signatures_file, sheet=paste0(signature_type,"_GRCh37"))
  if (signature_type=="SBS") {
    cosmic_signatures$Somatic.Mutation.Type <- paste0(substr(cosmic_signatures$Subtype, 1, 1),
                                                      "[",cosmic_signatures$Type, "]",
                                                      substr(cosmic_signatures$Subtype, 3, 3))
  } else {
    cosmic_signatures$Somatic.Mutation.Type <- cosmic_signatures$Type
  }
  etiology_data <- data.frame(Etiology=NA,
                              Etiology_category = NA,
                              row.names=c(),
                              stringsAsFactors = F)
  if (file.exists(etio_data_xlsx)) {
    if (signature_type %in% getSheetNames(etio_data_xlsx)) {
      etiology_data_raw <- read.xlsx(etio_data_xlsx, sheet=1)
      etiology_data <- data.frame(Etiology=etiology_data_raw$Category_specific,
                                  Etiology_category = etiology_data_raw$Category_broad,
                                  row.names=etiology_data_raw$signature,
                                  stringsAsFactors = F)
      
      categories <- sort(unique(etiology_data$Etiology_category))
      library(ggsci)
      # category_colors <- setNames(pal_rickandmorty((palette = c("schwifty")))(length(categories)), categories)
      category_colors <- setNames(pal_aaas((palette = c("default")))(length(categories)), categories)
      # category_colors <- setNames(pal_d3((palette = c("category10")))(length(categories)), categories)
      # category_colors <- setNames(pal_futurama((palette = c("planetexpress")))(length(categories)), categories)
      # library(ggthemes)
      # category_colors <- setNames(gdocs_pal()(length(categories)), categories)
      # category_colors <- setNames(calc_pal()(length(categories)), categories)
      
      # mycat=categories[1]
      etiology_colors <- list(Etiology=unlist(lapply(categories, function(mycat) {
        subcats <- unique(etiology_data$Etiology[etiology_data$Etiology_category==mycat])
        catcolor=category_colors[mycat]
        myrgb <- colorRamp(c(catcolor, "white"))(0:length(subcats)/length(subcats))[1:(length(subcats)), , drop=F]
        rownames(myrgb) <- subcats
        hexcolor <- apply(myrgb, 1, function(x) {rgb(x[1], x[2], x[3], maxColorValue = 255)})
        return(hexcolor)
      })))
      
      
      rowOrder=order(etiology_data$Etiology_category)
      etiology_data <- etiology_data[rowOrder, ,drop=F]
      etiology_data$annotation_color <- etiology_colors$Etiology[match(etiology_data$Etiology, names(etiology_colors$Etiology))]
      etiology_data$category_color <- category_colors[match(etiology_data$Etiology_category, names(category_colors))]
    }
  }
  return(list(signatures=cosmic_signatures,etio=etiology_data))
  
}


make_tnm <- function(mymaf,use_silent_mutations=F) {
  # mymaf=maf.filtered
  require(maftools)
  genome_build_res <- detect_maf_genome(mymaf)
  genome_build <- genome_build_res[[1]]
  genome_package <- genome_build_res[[3]]
  prefix_value=ifelse(genome_build_res[[2]],"chr","")
  
  tnm = trinucleotideMatrix(maf = mymaf,
                            ref_genome = genome_package,
                            # add = add_prefix,
                            prefix = prefix_value,
                            useSyn = use_silent_mutations)
  mut_mat = t(tnm$nmf_matrix)
  
}

detect_maf_genome <- function(maf) {
  require(maftools)
  if (! "NCBI_Build" %in% colnames(maf@data)) {
    warning("No genome information in MAF obj.")
    return(NA)
  }
  
  my_genome = paste0("GRCh",gsub("GRCh","",unique(maf@data$NCBI_Build)))
  if (length(my_genome) > 1) {
    warning("Multiple genomes listed in MAF obj. Trying the first one")
    my_genome <- my_genome[1]
  }
  
  return_genome <- switch(my_genome,GRCh38="hg38",GRCh37="hg19",GRCm38="mm10", NA)
  
  my_chrs <- unique(maf@data$Chromosome)
  add_chr = sum(grepl("^chr", my_chrs)) < length(my_chrs)
  
  pkg_prefix=ifelse(return_genome=="mm10","BSgenome.Mmusculus.UCSC.","BSgenome.Hsapiens.UCSC.")
  genome_package=paste0(pkg_prefix,return_genome)
  
  return(list(genome=return_genome, add_chr=add_chr, bsgenome_pkg=genome_package))
  
}


make_signature_plot <- function(tnm, savepath=NULL, signature_type="SBS") {
  # Basically stolen from: https://github.com/UMCUGenetics/MutationalPatterns/blob/master/R/plot_96_profile.R
  require(ggplot2)
  mydf <- data.frame(substitution = gsub(".*\\[(.*)\\].*","\\1",rownames(tnm)),
                     context=gsub("\\[.*\\]","\\.",rownames(tnm)),
                     apply(tnm,2,function(x){round(x/sum(x),digits = 5)}),
                     stringsAsFactors = F)
  plotdf <- reshape2::melt(mydf, id.vars=c("substitution","context"))
  plotdf <- plotdf[!is.nan(plotdf$value),]
  ymax=max(plotdf$value)
  COLORS6 <- c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE"
  )
  # browser()
  plot_96 = ggplot(data = plotdf, aes(x = context, y = value, 
                                   fill = substitution, width = 1)) + 
    geom_bar(stat = "identity", colour = "black", size = 0.2) + 
    scale_fill_manual(values = COLORS6) + 
    facet_grid(variable ~ substitution) + ylab("Relative contribution") + 
    coord_cartesian(ylim = c(0, ymax)) + scale_y_continuous(breaks = seq(0, 
                                                                         ymax, 0.1)) + guides(fill = FALSE) + theme_bw() + 
    theme(axis.title.y = element_text(size = 12, vjust = 1), 
          axis.text.y = element_text(size = 8), axis.title.x = element_text(size = 12), 
          axis.text.x = element_text(size = 5, angle = 90, 
                                     vjust = 0.4), strip.text.x = element_text(size = 9), 
          strip.text.y = element_text(size = 9), panel.grid.major.x = element_blank(), 
          panel.spacing.x = unit(0, "lines"))
  # plot
  if (!is.null(savepath)) {
    myheight=ncol(mydf)*0.8
    ggsave(savepath,plot=plot_96,height=myheight, width=8, limitsize=F)
    return(savepath)
  } else {
    return(plot_96)
  }
}


signature_plot_colorpalette <- function(signature_type="SBS") {
  
  COLORS6 <- c(
    "#2EBAED", "#000000", "#DE1C14",
    "#D4D2D2", "#ADCC54", "#F0D0CE"
  )
  
}


make_mut_signature_heatmap <- function(mymaf,use_silent_mutations=F, clinVarNames = NULL, 
                                       data_dir="data", full_output=F,
                                       # genome_build="hg19",
                                       clin_data=NULL, clin_data_colors=NULL,
                                       signatures_file=file.path("~","Documents","helper_functions","cosmic_data","COSMIC_Mutational_Signatures_v3.1.xlsx"),
                                       etio_data_xlsx=file.path("~","Documents","helper_functions","cosmic_data","COSMIC_signature_etiology.xlsx"),
                                       savename=NULL,
                                       progress_func=NULL) {
  require(maftools)
  require(MutationalPatterns)
  require(circlize)
  require(ComplexHeatmap)
  require(RColorBrewer)
  
  # browser()
  mut_mat <-  make_tnm(mymaf, use_silent_mutations = use_silent_mutations)
  
  cosmiclist <- load_cosmic_data(signatures_file, etio_data_xlsx)
  
  # sp_url <- file.path(data_dir,"cosmic","sigProfiler_exome_SBS_signatures.csv")
  # cosmic_signatures = read.table(sp_url, sep = ",", header = TRUE, stringsAsFactors = F)
  # cosmic_signatures$Somatic.Mutation.Type <- paste0(substr(cosmic_signatures$SubType, 1, 1),
  #                                                   "[",cosmic_signatures$Type, "]",
  #                                                   substr(cosmic_signatures$SubType, 3, 3))
  # 
  cosmic_signatures <- cosmiclist[["signatures"]]
  
  # Match the order of the mutation types to MutationalPatterns standard
  new_order = match(row.names(mut_mat), cosmic_signatures$Somatic.Mutation.Type)
  # Reorder cancer signatures dataframe
  cosmic_signatures = cosmic_signatures[as.vector(new_order),]
  # Add trinucletiode changes names as row.names
  row.names(cosmic_signatures) = cosmic_signatures$Somatic.Mutation.Type
  # Keep only 96 contributions of the signatures in matrix
  cosmic_signatures = as.matrix(cosmic_signatures[,grep("SBS*", colnames(cosmic_signatures))])
  
  hclust_cosmic = cluster_signatures(cosmic_signatures, method = "average")
  # store signatures in new order
  cosmic_order = colnames(cosmic_signatures)[hclust_cosmic$order]
  # plot(hclust_cosmic)
  
  if (is.function(progress_func)) {
    progress_func(value=50, detail = "Computing similarity to COSMIC...")
  }
  cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cosmic_signatures)
  
  # fit_res <- fit_to_signatures(mut_mat, cosmic_signatures)
  
  # data_dir="data"
  # etio_data_file = file.path(data_dir, "cosmic","COSMIC_signature_etiology.xlsx")
  # etiology_data_raw <- read.xlsx(etio_data_file, sheet="final categories")
  # etiology_data <- as.character(etiology_data_raw$CATEGORY)
  # etiology_data <- data.frame(Etiology=etiology_data_raw$CATEGORY, row.names=etiology_data_raw$signature)
  etio <- cosmiclist[["etio"]]
  
  etiology_data <- data.frame(Etiology=paste0("[",etio$Etiology_category,"] ", etio$Etiology),
                              row.names=rownames(etio))
  # names(etiology_data) <- etiology_data_raw$signature
  # etiology_colors <-  list(Etiology=c("APOBEC" = "#fce116",
  #                                     "Defective DNA Mismatch Repair" = "#31A354",
  #                                     "Defective DNA Repair" = "#A1D99B",
  #                                     "Exonuclease Domain" = "#E5F5E0",
  #                                     "Exposure to Alfatoxin" = "#DE2D26",
  #                                     "Exposure to Aristolochic Acid" = "#FC9272",
  #                                     "Exposure to Haloalkanes" = "#FEE0D2",
  #                                     "Tobacco - Chewing" = "#6d3617",
  #                                     "Tobacco - Smoking" = "#a85423",
  #                                     "Tobacco - Smoking Associated" = "#d87841",
  #                                     "Prior Therapy - Alkylating Agents" = "#2171B5",
  #                                     "Prior Therapy - Platinum Drugs" = "#6BAED6",
  #                                     "Prior Therapy - Immunosuppression" = "#BDD7E7",
  #                                     "Prior Therapy Associated" = "#EFF3FF",
  #                                     "ROS Damage" = "#BCBDDC",
  #                                     "UV Light" = "#756BB1",
  #                                     "UV Light Associated" = "#a29bca",
  #                                     "Unknown" = "grey70",
  #                                     "Possible Sequencing Artifact" = "grey50")
  # )
  etiology_colors <- list(Etiology=setNames(etio$annotation_color, etiology_data$Etiology))
  plot_matrix <- t(cos_sim_samples_signatures)
  
  if (is.function(progress_func)) {
    progress_func(value=80, detail = "Making heatmap")
  }
  
  # browser()
  rowOrder=order(etiology_data)
  etiology_data <- etiology_data[rowOrder,1,drop=F]
  mycolors <- etiology_colors
  signature_anno <- rowAnnotation(df=etiology_data, 
                                  name="Signature Anno", col=etiology_colors, show_annotation_name = FALSE)
  
  plot_matrix <- plot_matrix[match(rownames(etiology_data),rownames(plot_matrix)),]
  
  myanno=NULL
  if (!is.null(clin_data)) {
    anno_data <- data.frame(clin_data[match(colnames(plot_matrix), clin_data$Tumor_Sample_Barcode),],stringsAsFactors = F)
    row.names(anno_data) <- anno_data$Tumor_Sample_Barcode
    anno_data <- anno_data[,!colnames(anno_data) %in% "Tumor_Sample_Barcode"]
    myanno <- HeatmapAnnotation(df=anno_data,col = clin_data_colors)
  }
  
  add_sample_names <- ifelse(ncol(plot_matrix)>10, F, T)
  plot_matrix[is.nan(plot_matrix)] <- 0
  myHM <- Heatmap(plot_matrix, 
                  col=colorRamp2(seq(min(plot_matrix), max(plot_matrix), length.out = 20),colorRampPalette(brewer.pal(9, "BuGn"))(20)),
                  left_annotation = signature_anno,
                  bottom_annotation = myanno,
                  cluster_rows = F, row_order = rowOrder,
                  clustering_method_rows = "median",
                  clustering_method_columns = "median",
                  heatmap_height = unit(6, "inches"),
                  heatmap_legend_param = list(
                    # at = c(-2, 0, 2),
                    # labels = c("low", "zero", "high"),
                    title = "Cosine Similarity",
                    # legend_height = unit(4, "cm"),
                    legend_direction = "horizontal"
                  ),
                  show_row_names=T, row_names_gp = gpar(fontsize = 5),
                  show_column_names = add_sample_names)
  
  if ( ! is.null(savename) ) {
    # save_name <- paste0(out_dir,"/oncoplot.",cohort_freq_thresh,".pdf")
    # browser()
    anno_height=ifelse(!is.null(clin_data), min(c(4, 0.5*ncol(clin_data))), 0)
    hm_height=max(round(0.1*nrow(plot_matrix),0),4) + anno_height
    hm_width=hm_height*0.75 + anno_height*1.2
    pdf(file = savename,height=hm_height,width=hm_width)
    draw(myHM)
    dev.off()
  }
  
  ### Return the oncoplot (if function is pointed to a variable)
  returnval <- myHM
  if (full_output) {
    output_list <- list(plot_matrix=plot_matrix,
                        signature_annotations=signature_anno,
                        etiology_data=etiology_data, 
                        heatmap_obj=myHM)
    returnval <- output_list
  }
  
  invisible(myHM)
}
