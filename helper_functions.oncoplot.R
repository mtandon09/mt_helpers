
filter_maf <- function(maf_file, flag_genes="default",save_name=NULL,no_filter=F,
                       norm_alt_max=1,t_alt_min=1,t_depth_min=20, 
                       tumor_freq_min=0.05, norm_freq_max=0.02,
                       gnomAD_AF_max=0.001, AF_max=0.001, ExAC_AF_max=0.001, 
                       n_callers=2, variant_caller=NULL) {
  # browser()
  if (length(flag_genes)==0) {
    flag_genes <- c()
  } else if (flag_genes[1]=="default") {
    flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
  }
  if (class(maf_file)=="character") {
    if (file.exists(maf_file[1])) {
      maf_df.raw <- read.table(maf_file, sep="\t", header=T, fill = T, quote="\"", stringsAsFactors = F)
    } else {
      stop("Interpreting 'maf_file' as a file name, but file not found.")
    }
  } else if (class(maf_file)=="MAF"){
    maf_df.raw <- rbind(maf_file@data, maf_file@maf.silent)
  }
  maf_df.raw <- maf_df.raw[maf_df.raw$Hugo_Symbol != "Hugo_Symbol",]
  filter_genes=!maf_df.raw$Hugo_Symbol %in% flag_genes
  maf_df.raw <- maf_df.raw[filter_genes,]
  
  
  if (!"tumor_freq" %in% colnames(maf_df.raw)) {
    maf_df.raw$tumor_freq <- as.numeric(maf_df.raw$t_alt_count)/as.numeric(maf_df.raw$t_depth)
  }  
  if (!"norm_freq" %in% colnames(maf_df.raw)) {
      maf_df.raw$norm_freq <- as.numeric(maf_df.raw$n_alt_count)/as.numeric(maf_df.raw$n_depth)
  }
  
  filter_tumor_depth=rep(TRUE,nrow(maf_df.raw))
  filter_tumor_alt=rep(TRUE,nrow(maf_df.raw))
  filter_norm_alt=rep(TRUE,nrow(maf_df.raw))
  filter_pop_freq=rep(TRUE,nrow(maf_df.raw))
  
  if (!no_filter) {
    options(warn=-1)
    filter_tumor_depth=as.numeric(maf_df.raw$t_depth) > t_depth_min
    if (!sum(is.na(maf_df.raw$norm_freq)) == nrow(maf_df.raw)){
      filter_norm_alt=maf_df.raw$norm_freq < norm_freq_max
    } 
    filter_tumor_alt=maf_df.raw$tumor_freq > tumor_freq_min
    if (! is.null(t_alt_min)){
      filter_tumor_alt <- filter_tumor_alt & maf_df.raw$t_alt_count > t_alt_min
    }
    
    filter_gnomad=rep(TRUE,nrow(maf_df.raw))
    filter_1000G=rep(TRUE,nrow(maf_df.raw))
    filter_exac=rep(TRUE,nrow(maf_df.raw))
    if (!is.null(maf_df.raw$gnomAD_AF)) {
      filter_gnomad=maf_df.raw$gnomAD_AF %in% c("-","") | is.na(maf_df.raw$gnomAD_AF) | as.numeric(maf_df.raw$gnomAD_AF) < min(gnomAD_AF_max,1)
    }
    if (!is.null(maf_df.raw$AF)) {
      filter_1000G=maf_df.raw$AF %in% c("-","") | is.na(maf_df.raw$AF)  | as.numeric(maf_df.raw$AF) < min(AF_max,1)
    }
    if (!is.null(maf_df.raw$ExAC_AF)) {
      filter_exac=maf_df.raw$ExAC_AF %in% c("-","") | is.na(maf_df.raw$ExAC_AF) | as.numeric(maf_df.raw$ExAC_AF) < min(ExAC_AF_max,1)
    }
    filter_pop_freq=filter_gnomad & filter_1000G & filter_exac
    options(warn=0)
  }
  filter_caller=rep(TRUE,nrow(maf_df.raw))
  if (! is.null(variant_caller)) {       ### Set 'variant_caller' to NULL to skip any filtering based on caller
    maf_df.raw$set[maf_df.raw$set=="" & maf_df.raw$Hugo_Symbol=="Hugo_Symbol"] <- "set"
    maf_df.raw$set[maf_df.raw$set==""] <- "N.A."
    if (variant_caller == "consensus") {   ### Set 'variant_caller' to 'consensus' to keep variants by two or more callers
      # filter_caller <- grepl("-|Intersection", maf_df.raw$set)
      filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {length(x)>=n_callers | "Intersection" %in% x}))
    } else {                             ### Set 'variant_caller' to one of the callers (mutect, mutect2, vardict, or strelka) to get only that caller
      # filter_caller <- grepl(paste0(variant_caller,"[|-]|Intersection"), maf_df.raw$set)
      filter_caller <- unlist(lapply(strsplit(maf_df.raw$set,"-"), function(x) {any(c(variant_caller,"Intersection") %in% x)}))
    }
  }
  
  maf_df.rawest <- maf_df.raw
  maf_df.raw <- maf_df.raw[filter_tumor_depth & filter_norm_alt & filter_tumor_alt & filter_pop_freq & filter_caller,]
  # browser()
  maf_df.raw <- maf_df.raw[rowSums(is.na(maf_df.raw))!=ncol(maf_df.raw),]
  
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    write.table(maf_df.raw, sep="\t", quote=F, file = save_name, row.names = F)
    print(paste0("Saving filtered maf to ",save_name))
    return(save_name)
  } else {
    return(maf_df.raw)
  }
}

filter_maf_tbl <- function(maftbl, 
                           flag_genes="default",
                           #save_name=NULL,
                           no_filter=F,
                           grep_vcf_filter_col="PASS",
                           non_silent_only=F,
                           t_alt_min=2,
                           t_alt_max=1e12,
                           t_depth_min=5,
                           t_depth_max=1e12, 
                           tumor_freq_min=0.01,
                           tumor_freq_max=1,
                           n_alt_min=0,
                           n_alt_max=1,
                           n_depth_min=0,
                           n_depth_max=1e12,
                           norm_freq_min=0,
                           norm_freq_max=0.02,
                           gnomAD_AF_min=0,
                           gnomAD_AF_max=0.001,
                           AF_min=0,
                           AF_max=0.001,
                           ExAC_AF_min=0, 
                           ExAC_AF_max=0.001, 
                           n_callers=2,
                           variant_caller=NULL) {
  
  if (length(flag_genes)==0) {
    flag_genes <- c()
  } else if (flag_genes[1]=="default") {
    flag_genes <- c("TTN","MUC16","OBSCN","AHNAK2","SYNE1","FLG","MUC5B","DNAH17","PLEC","DST","SYNE2","NEB","HSPG2","LAMA5","AHNAK","HMCN1","USH2A","DNAH11","MACF1","MUC17","DNAH5","GPR98","FAT1","PKD1","MDN1","RNF213","RYR1","DNAH2","DNAH3","DNAH8","DNAH1","DNAH9","ABCA13","SRRM2","CUBN","SPTBN5","PKHD1","LRP2","FBN3","CDH23","DNAH10","FAT4","RYR3","PKHD1L1","FAT2","CSMD1","PCNT","COL6A3","FRAS1","FCGBP","RYR2","HYDIN","XIRP2","LAMA1")
  }
  # browser()
  require(tibble)
  require(dplyr)
  df <- as_tibble(maftbl)
  maf_df.raw <- df[df$Hugo_Symbol != "Hugo_Symbol",]
  maf_df.raw <- maf_df.raw[!df$Hugo_Symbol %in% flag_genes,]
  
  if ("FILTER" %in% colnames(maf_df.raw)) {
    if (!is.null(grep_vcf_filter_col)) {
      # maf_df.raw <- maf_df.raw[grepl(grep_vcf_filter_col,pull(maf_df.raw,FILTER)),]
      maf_df.raw <- maf_df.raw[grepl(grep_vcf_filter_col,maf_df.raw$FILTER),]
    }
  } else {
    message("FILTER column not found; skipping...")
  }
  
  #Variant Classification with High/Moderate variant consequences. http://asia.ensembl.org/Help/Glossary?id=535
  vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                   "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                   "In_Frame_Ins", "Missense_Mutation")
  if ("Variant_Classification" %in% colnames(maf_df.raw)) {
    if (non_silent_only) {
      # maf_df.raw <- maf_df.raw[grepl(grep_vcf_filter_col,pull(maf_df.raw,FILTER)),]
      maf_df.raw <- maf_df.raw[maf_df.raw$Variant_Classification %in% vc.nonSilent,]
    }
  } else {
    message("Variant_Classification column not found; skipping...")
  }
  
  # if (!"tumor_freq" %in% colnames(maf_df.raw)) {
  #   # if (! all(c("t_alt_count","t_depth") %in% colnames(maf_df.raw))) {
  #     # stop("Can't find t_alt_count or t_depth columns")
  #   # }
  #   if (! "t_alt_count" %in% colnames(maf_df.raw)) {
  #     maf_df.raw$t_alt_count <- NA
  #   }
  #   if (! "t_alt_count" %in% colnames(maf_df.raw)) {
  #     maf_df.raw$t_depth <- NA
  #   }
  #   maf_df.raw$tumor_freq <- as.numeric(maf_df.raw$t_alt_count)/as.numeric(maf_df.raw$t_depth)
  # }  
  # if (!"norm_freq" %in% colnames(maf_df.raw)) {
  #   # if (! all(c("n_alt_count","n_depth") %in% colnames(maf_df.raw))) {
  #   #   maf_df.raw$norm_freq <- rep(0, nrow(maf_df.raw))
  #   # } else {
  #   #   maf_df.raw$norm_freq <- as.numeric(maf_df.raw$n_alt_count)/as.numeric(maf_df.raw$n_depth)
  #   # }
  #   if (! "n_alt_count" %in% colnames(maf_df.raw)) {
  #     maf_df.raw$t_alt_count <- NA
  #   }
  #   if (! "n_depth" %in% colnames(maf_df.raw)) {
  #     maf_df.raw$t_depth <- NA
  #   }
  #   maf_df.raw$norm_freq <- as.numeric(maf_df.raw$n_alt_count)/as.numeric(maf_df.raw$n_depth)
  # }
  
  
  maf_num_filter_columns <- list("t_alt_count"=c(min=t_alt_min, max=t_alt_max),
                                 "t_depth"=c(min=t_depth_min, max=t_depth_max),
                                 "tumor_freq"=c(min=tumor_freq_min, max=tumor_freq_max),
                                 "n_alt_count"=c(min=n_alt_min, max=n_alt_max),
                                 "n_depth"=c(min=n_depth_min, max=n_depth_max),
                                 "norm_freq"=c(min=norm_freq_min, max=norm_freq_max),
                                 "gnomAD_AF"=c(min=gnomAD_AF_min, max=gnomAD_AF_max),
                                 "AF"=c(min=AF_min, max=AF_max),
                                 "ExAC_AF"=c(min=ExAC_AF_min, max=ExAC_AF_max)
  )
  # browser()
  numfilter_columns <- names(maf_num_filter_columns)[names(maf_num_filter_columns) %in% colnames(maf_df.raw)]
  notfound <- setdiff(names(maf_num_filter_columns), numfilter_columns)
  if (length(notfound) > 0 ) {
    message(paste0("Couldn't find these columns; skipping filtering for these: ", paste0(notfound, collapse=",")))
  }
  
  return_df <- maf_df.raw
  if (length(numfilter_columns)>0) {
    all_num_filters <- lapply(numfilter_columns, function(col_name) {
      currdata <- as.numeric(pull(maf_df.raw,col_name))
      currdata[is.na(currdata)] <- 0
      filter_vec <- currdata >= maf_num_filter_columns[[col_name]]["min"] & currdata <= maf_num_filter_columns[[col_name]]["max"]
      return(filter_vec)
    })
    # browser()
    final_num_filters <- Reduce("&", all_num_filters)
    
    return_df <- maf_df.raw[final_num_filters,]
  }
  return(return_df)
}



filter_maf_chunked <- function(maf, chunk_lines=10000,...) {
  
  clindata <- NULL
  if ("MAF" %in% class(maf)) {
    filtered_df <- filter_maf_tbl(
                            rbind(maf@data, maf@maf.silent),
                            ...)
    clindata <- maf@clinical.data
    
  } else if (file.exists(maf)) {
    require(readr)
    readr_filterfunc <- function(df, pos) {
      filter_maf_tbl(df,...)
    }
    filtered_df <- read_tsv_chunked(maf,chunk_size = chunk_lines, col_types = cols(), callback = DataFrameCallback$new(readr_filterfunc), comment="#")
  } else {
    stop(paste0("Don't know what to do with input type '",class(maf),"'"))
  }
  
  maf.filtered <- read.maf(filtered_df, clinicalData = clindata)
  return(maf.filtered)

  
}


make_oncoplot <- function(maf.filtered, cohort_freq_thresh = 0.01, auto_adjust_threshold=T,
                          oncomat_only=F,
                          clin_data=NULL, clin_data_colors=NULL,
                          savename=NULL) {
  require(ComplexHeatmap)
  ### Read in MAF file
  # maf.filtered <- read.maf(maf_file)
  
  ### Structure info about the fraction of the cohort that has each gene mutated
  frac_mut <- data.frame(Hugo_Symbol=maf.filtered@gene.summary$Hugo_Symbol,
                         frac_mut=(maf.filtered@gene.summary$MutatedSamples/as.numeric(maf.filtered@summary$summary[3])),
                         stringsAsFactors = F)
  
  ngene_max=25
  target_frac = sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))]
  if (auto_adjust_threshold) {
    cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
  }
  ### Select genes based on the frequency threshold
  freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut >= cohort_freq_thresh]
  freq_genes <- freq_genes[1:ngene_max]
  if (length(freq_genes) == 0) {
    stop("No genes to plot; change the frequency threshold to include more genes.")
  }
  if (length(freq_genes) > 100) {
    target_frac = round(sort(frac_mut$frac_mut, decreasing = T)[min(ngene_max,nrow(frac_mut))],2)
    # stop(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
    warning(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
    # return(NA)
  }
  gene_list <- list(freq_genes)
  reasons <- paste0("Cohort Freq > ",cohort_freq_thresh)
  
  ### Collect genes to plot
  genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), reason=c())
  for (i in 1:length(gene_list)) {
    if (is.na(gene_list[[i]][1])) {
      next
    }
    genes_for_oncoplot <- rbind(genes_for_oncoplot,
                                data.frame(Hugo_Symbol=gene_list[[i]],
                                           reason=reasons[i]))
  }
  genes_for_oncoplot <- cbind(genes_for_oncoplot,
                              frac=frac_mut$frac_mut[match(genes_for_oncoplot$Hugo_Symbol, frac_mut$Hugo_Symbol)])
  
  genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$reason, -genes_for_oncoplot$frac),]
  
  ### Split the oncoplot based on the reason for picking the gene
  ###   Here, we're only picked based on the frequency
  ###   But this framework is useful for plotting genes picked using various criteria
  split_idx=genes_for_oncoplot$reason
  split_colors <- rainbow(length(levels(split_idx)))
  names(split_colors) <- as.character(genes_for_oncoplot$reason[!duplicated(genes_for_oncoplot$reason)])
  split_colors <- list(Reason=split_colors)
  
  # source("scripts/helper_functions.oncoplot.R")
  ### Make matrix to plot, and order it correctly
  oncomat <- createOncoMatrix(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = F)$oncoMatrix
  oncomat <- oncomat[match(genes_for_oncoplot$Hugo_Symbol,rownames(oncomat)), ]
  onco_genes <- rownames(oncomat)
  
  if (oncomat_only) {
    return(oncomat)
  }
  oncomat.plot <- oncomat
  
  ### Set the height of the plot based on number of genes
  onco_height=NULL
  if (is.null(onco_height)) {
    onco_height=max(round(0.2*nrow(oncomat.plot),0),5)
  }
  
  ### Make the mutation type names prettier by removing the underscore
  # my_mut_col <- mutation_colors
  # names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
  oncomat.plot <- gsub("_"," ",oncomat.plot)
  
  ### Column labels get cluttered if too many samples
  show_sample_names=T
  if (ncol(oncomat.plot) > 20) {
    show_sample_names=F
  }
  
  # browser()
  myanno=NULL
  if (!is.null(clin_data)) {
    
    myanno <- make_column_annotation(clin_data,colnames(oncomat.plot), clin_data_colors)
    # print(myanno)
  }
  
  ## Show total burden for top annotation
  variant_type_data <- data.frame(maf.filtered@variant.classification.summary)
  rownames(variant_type_data) <- variant_type_data$Tumor_Sample_Barcode
  colnames(variant_type_data) <- gsub("_"," ",colnames(variant_type_data))
  variant_type_data <- variant_type_data[,c(-1,-ncol(variant_type_data))]
  variant_type_data <- variant_type_data[match(colnames(oncomat.plot), rownames(variant_type_data)),
                                         rev(order(colSums(variant_type_data)))]
  # browser()
  var_anno_colors <- mutation_colors[match(colnames(variant_type_data), names(mutation_colors))]
  top_ha = HeatmapAnnotation("Total\nMutations" = anno_barplot(variant_type_data, gp = gpar(fill = var_anno_colors), border = F),
                             annotation_name_side = "left",annotation_name_rot=90,annotation_name_gp = gpar(cex=0.7))
  
  # browser()
  
  pct_anno <- paste0(prettyNum(frac_mut$frac_mut[match(onco_genes, frac_mut$Hugo_Symbol)]*100,digits=1),"%")
  left_ha = rowAnnotation("Cohort Pct"=anno_text(pct_anno,gp = gpar(cex=0.7)), show_annotation_name=F)
  # print(oncomat.plot)
  ### Make the oncoplot
  onco_base_default <- oncoPrint(oncomat.plot, alter_fun = alter_fun, col=mutation_colors, row_order=1:nrow(oncomat.plot),
                                 name="oncoplot",
                                 show_pct = F,
                                 row_split=split_idx,
                                 row_title = NULL,
                                 bottom_annotation = myanno,
                                 top_annotation = top_ha,
                                 left_annotation = left_ha,
                                 show_column_names = show_sample_names)#,
  
  if ( ! is.null(savename) ) {
    # save_name <- paste0(out_dir,"/oncoplot.",cohort_freq_thresh,".pdf")
    onco_height=max(round(0.15*nrow(oncomat.plot),0),6)
    onco_width=onco_height*0.75
    pdf(file = savename,height=onco_height,width=onco_width)
    draw(onco_base_default)
    dev.off()
  }
  
  ### Return the oncoplot (if function is pointed to a variable)
  invisible(onco_base_default)
}


make_oncoplot2 <- function(maf.filtered, cohort_freq_thresh = 0.1,
                          use_clinvar_anno=F, title_text="",
                          genes_to_plot=NULL, include_all=F,
                          custom_column_order=NULL, full_output=F,
                          savename=NULL) {
  
  ### Read in MAF file
  # maf.filtered <- read.maf(maf_file)
  
  if (! is.null(genes_to_plot)) {
    if (class(genes_to_plot)=="character") {
      if (length(genes_to_plot)==1) {
        ## Then it's either a file name or a single gene
        if (file.exists(genes_to_plot)) {
          ### Need to parse file type and read accordingly; assuming tsv for now
          gene_data <- read.table(genes_to_plot,sep="\t", header=T)
          if (sum(c("Hugo_Symbol","Reason") %in% colnames(gene_data)) != 2) {
            stop("Can't find Hugo Symbol or Reason in custom gene input.")
          }
          genes_for_oncoplot <- gene_data
        } else {
          stop(paste0("Can't find file: ",genes_to_plot))
        }
      } else {
        genes_for_oncoplot <- data.frame(Hugo_Symbol=genes_to_plot,Reason="Selected Genes")
      }
    } else if (class(genes_to_plot)=="data.frame") {
      genes_for_oncoplot <- genes_to_plot
      if (! "Reason" %in% colnames(genes_for_oncoplot)) {
        genes_for_oncoplot$Reason <- "Selected Genes"
      }
    } else {
      stop(paste0("Don't know what to do with 'genes_to_plot' of class: ",class(genes_to_plot)))
    }
    
    genes_for_oncoplot <- genes_for_oncoplot[,c("Hugo_Symbol","Reason")]
    genes_for_oncoplot <- data.frame(apply(genes_for_oncoplot,2,as.character), stringsAsFactors = F)
    genes_for_oncoplot <- genes_for_oncoplot[genes_for_oncoplot$Hugo_Symbol %in% maf.filtered@gene.summary$Hugo_Symbol, ]
  } else {
    genes_for_oncoplot <- data.frame(Hugo_Symbol=c(), Reason=c())
  }
  
  frac_mut <- data.frame(Hugo_Symbol=maf.filtered@gene.summary$Hugo_Symbol,
                         frac_mut=(maf.filtered@gene.summary$MutatedSamples/as.numeric(maf.filtered@summary$summary[3])),
                         stringsAsFactors = F)
  if (! is.null(cohort_freq_thresh)) {
    ### Structure info about the fraction of the cohort that has each gene mutated
    # browser()
    # target_frac = sort(frac_mut$frac_mut, decreasing = T)[min(50,nrow(frac_mut))]
    # cohort_freq_thresh <- max(c(cohort_freq_thresh,target_frac))
    ### Select genes based on the frequency threshold
    freq_genes <- frac_mut$Hugo_Symbol[frac_mut$frac_mut > cohort_freq_thresh]
    if (length(freq_genes) == 0) {
      stop("No genes to plot; change the frequency threshold to include more genes.")
    }
    if (length(freq_genes) > 200) {
      target_frac = round(sort(frac_mut$frac_mut, decreasing = T)[min(50,nrow(frac_mut))],2)
      stop(paste0("Too many genes for oncoplot. Trying setting the cohort mutated fraction to > ", target_frac))
      # return(NA)
    }
    gene_list <- list(freq_genes)
    reasons <- paste0("Cohort Freq > ",cohort_freq_thresh)
    
    ### Collect genes to plot
    
    for (i in 1:length(gene_list)) {
      if (is.na(gene_list[[i]][1])) {
        next
      }
      genes_for_oncoplot <- rbind(genes_for_oncoplot,
                                  data.frame(Hugo_Symbol=gene_list[[i]],
                                             Reason=reasons[i]))
    }
    genes_for_oncoplot <- cbind(genes_for_oncoplot,
                                frac=frac_mut$frac_mut[match(genes_for_oncoplot$Hugo_Symbol, frac_mut$Hugo_Symbol)])
    
    genes_for_oncoplot <- genes_for_oncoplot[order(genes_for_oncoplot$Reason, -genes_for_oncoplot$frac),]
  }
  # browser()
  ### Split the oncoplot based on the reason for picking the gene
  ###   Here, we're only picked based on the frequency
  ###   But this framework is useful for plotting genes picked using various criteria
  split_idx=factor(genes_for_oncoplot$Reason)
  split_colors <- rainbow(length(levels(split_idx)))
  names(split_colors) <- as.character(genes_for_oncoplot$Reason[!duplicated(genes_for_oncoplot$Reason)])
  split_colors <- list(Reason=split_colors)
  
  # source("scripts/helper_functions.oncoplot.R")
  ### Make matrix to plot, and order it correctly
  if (use_clinvar_anno) {
    oncomat <- createOncoMatrix_CLINSIG(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = include_all)$oncoMatrix
  } else {
    oncomat <- createOncoMatrix(maf.filtered, g=genes_for_oncoplot$Hugo_Symbol, add_missing = include_all)$oncoMatrix
  }
  oncomat <- oncomat[match(genes_for_oncoplot$Hugo_Symbol,rownames(oncomat)), ]
  onco_genes <- rownames(oncomat)
  
  # browser()
  if (include_all) {
    ### createOncoMatrix drops empty samples, so this adds them back in
    # all_wes_samples <- as.character(sample_info.exome$Tumor_Sample_Barcode[!is.na(sample_info.exome$Tumor_Sample_Barcode)])
    all_wes_samples <- levels(maf.filtered@variants.per.sample$Tumor_Sample_Barcode)
    extra_samples <- setdiff(all_wes_samples, colnames(oncomat) )
    print(paste0("Adding back ", length(extra_samples), " samples with no reported mutations..."))
    empty_data <- matrix(data = "", nrow=nrow(oncomat), ncol=length(extra_samples), dimnames=list(rownames(oncomat), extra_samples))
    oncomat <- cbind(oncomat, empty_data)
  }
  
  if (!is.null(custom_column_order)) {
    custom_order <- match(custom_column_order, colnames(oncomat), nomatch=0)
    oncomat <- oncomat[,custom_order]
  }
  oncomat.plot <- oncomat
  
  ### Set the height of the plot based on number of genes
  onco_height=NULL
  if (is.null(onco_height)) {
    onco_height=max(round(0.2*nrow(oncomat.plot),0),6)
  }
  
  ### Make the mutation type names prettier by removing the underscore
  # my_mut_col <- mutation_colors
  # names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
  oncomat.plot <- gsub("_"," ",oncomat.plot)
  
  ### Column labels get cluttered if too many samples
  show_sample_names=T
  if (ncol(oncomat.plot) > 20) {
    show_sample_names=F
  }
  
  # if (show_burden) {
  variant_type_data <- data.frame(maf.filtered@variant.classification.summary)
  unmut_samples <- setdiff(levels(variant_type_data$Tumor_Sample_Barcode),variant_type_data$Tumor_Sample_Barcode)
  
  if (length(unmut_samples) > 0) {
    unmut_data <- data.frame(unmut_samples, matrix(0, nrow=length(unmut_samples), ncol=ncol(variant_type_data)-1))
    colnames(unmut_data) <- colnames(variant_type_data)
    variant_type_data <- rbind(variant_type_data,unmut_data)
  }
  
  rownames(variant_type_data) <- variant_type_data$Tumor_Sample_Barcode
  colnames(variant_type_data) <- gsub("_"," ",colnames(variant_type_data))
  variant_type_data <- variant_type_data[,c(-1,-ncol(variant_type_data))]
  variant_type_data <- variant_type_data[match(colnames(oncomat.plot), rownames(variant_type_data)),
                                         rev(order(colSums(variant_type_data)))]
  # browser()
  var_anno_colors <- mutation_colors[match(colnames(variant_type_data), names(mutation_colors))]
  ha = HeatmapAnnotation(Burden = anno_barplot(variant_type_data, gp = gpar(fill = var_anno_colors), border = F))
  # }
  
  pct_anno <- paste0(prettyNum(frac_mut$frac_mut[match(onco_genes, frac_mut$Hugo_Symbol)]*100,digits=1),"%")
  left_ha = rowAnnotation("Cohort Pct"=anno_text(pct_anno,gp = gpar(cex=0.7)), show_annotation_name=F)
  # print(oncomat.plot)
  ### Make the oncoplot
  # browser()
  
  onco_base_default <- oncoPrint(oncomat.plot, alter_fun = alter_fun, col=mutation_colors, row_order=1:nrow(oncomat.plot),
                                 name="oncoplot",
                                 row_title=title_text,
                                 show_pct = F,
                                 row_split=split_idx,
                                 # row_title = NULL,
                                 left_annotation = rowAnnotation(Reason = split_idx, col=split_colors, annotation_width = unit(0.3, "mm"), show_annotation_name=F),
                                 right_annotation = left_ha,
                                 top_annotation = ha,
                                 show_column_names = show_sample_names)#,
  # column_names_rot = 30,
  # column_gap = unit(0.0001,"npc"),
  # width = unit(0.75, "npc"))
  # print(dim(oncomat))
  ### Save the oncoplot
  if (full_output) {
    return_val <- list(onco_base_default, oncomat.plot)
  } else {
    return_val <- onco_base_default
  }
  
  if ( ! is.null(savename) ) {
    # save_name <- paste0(out_dir,"/oncoplot.",cohort_freq_thresh,".pdf")
    onco_height=max(round(0.15*nrow(oncomat.plot),0),6)
    onco_width=onco_height*0.75
    pdf(file = savename,height=onco_height,width=onco_width)
    draw(onco_base_default)
    dev.off()
  }
  
  invisible(return_val)
}


make_burden_plot <- function(maf.filtered, plotType=NULL, mb_covered=NULL, save_data_to_file=NULL) {
  
  require(dplyr)
  num_var_data <- maf.filtered@variants.per.sample
  colnames(num_var_data) <- c("Tumor_Sample_Barcode","Variants_filtered")
  num_var_data$mut_burden_count <- num_var_data$Variants_filtered
  num_var_data$mut_burden <- num_var_data$mut_burden_count
  y_label_text="Mutation Count"
  if (is.numeric(mb_covered)) {
    print("Normalizing mutation count by covered bases...")
    num_var_data$mut_burden <- num_var_data$mut_burden/mb_covered
    y_label_text="Mutation Burden (mutations/Mb)"
  }
  
  nsamples=nrow(num_var_data)
  if (is.null(plotType)) {
    plotType <- ifelse(nrow(num_var_data) > 15, "Dotplot", "Barplot")
    print(paste0("Using plot type: ", plotType))
  }
  # browser()
  ## Re-jigger the factor levels so they're ordered by decreasing mutation burden (for the plotting)
  num_var_data$Tumor_Sample_Barcode <- factor(num_var_data$Tumor_Sample_Barcode,
                                              levels=num_var_data$Tumor_Sample_Barcode[order(num_var_data$mut_burden, decreasing = T)])
  
  ########################################################
  #### 5. Generate plots for mutation burden
  
  ## Pick colors
  median_mut_burdens <- num_var_data %>% summarise(median=median(mut_burden))
  
  num_var_data$xlabel <- factor(num_var_data$xlabel,
                                levels=num_var_data$xlabel[order(num_var_data$mut_burden, decreasing = T)])
  num_var_data$hoverlabel <- paste0("Sample: ",num_var_data$Tumor_Sample_Barcode,"\nMutations: ", num_var_data$mut_burden)
  
  ### Mutation burden stacked with variant classification counts
  ### Works better for smaller cohorts, say < 20
  variant_type_per_sample <- as.data.frame(maf.filtered@variant.classification.summary)
  var_type.melt <- reshape2::melt(variant_type_per_sample, id.vars="Tumor_Sample_Barcode",variable.name="classification",value.name="mutation_count")
  var_type.melt$mut_burden <- var_type.melt$mutation_count
  if (is.numeric(mb_covered)) {
    var_type.melt$mut_burden <- var_type.melt$mut_burden/mb_covered
  }
  median_mut_burdens <- data.frame(median=median(var_type.melt[var_type.melt$classification== "total","mut_burden"]))
  
  plotdata <- var_type.melt[var_type.melt$classification != "total",]
  plotdata$Tumor_Sample_Barcode <- factor(as.character(plotdata$Tumor_Sample_Barcode),
                                          levels=variant_type_per_sample$Tumor_Sample_Barcode[order(variant_type_per_sample$total, decreasing = T)])
  plotdata$classification <- gsub("_"," ",plotdata$classification)
  
  class_means <- plotdata %>% group_by(classification) %>% summarise(mean=mean(mut_burden))
  plotdata$classification <- factor(as.character(plotdata$classification),
                                    levels=class_means$classification[order(class_means$mean, decreasing = F)])
  
  my_class_colors <- mutation_colors
  
  plotdata$hoverlabel <- paste0("Sample: ",plotdata$Tumor_Sample_Barcode,"\nMutations: ", plotdata$mut_burden)
  
  if (plotType=="Barplot") {
    if (length(unique(plotdata$Tumor_Sample_Barcode)) <= 20) {
      xaxis_text <- element_text(angle=30, hjust=1)
    } else {
      xaxis_text <- element_blank()
    }
    burden_plot <- ggplot(plotdata, aes(x=Tumor_Sample_Barcode, y=mut_burden, text=hoverlabel)) +
      geom_bar(aes(fill=classification), stat="identity",width=1,size=0.3, color="black") +
      scale_fill_manual(values=my_class_colors) +
      theme_linedraw(base_size = 12) +
      xlab("") + ylab(y_label_text) +
      # geom_hline(data = median_mut_burdens, aes(yintercept=median),linetype="dashed", color="grey60") +  ### Is screwed up with ggplotly
      theme(
        axis.text.x = xaxis_text,
        axis.ticks.x = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=8),
        legend.title = element_blank(),
        legend.key.height = unit(0.01,"npc"),
        legend.key.width =  unit(0.02,"npc"),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank())
    if (nsamples > 1) {
      burden_plot <- burden_plot + geom_hline(data = median_mut_burdens, aes(yintercept=median),linetype="dashed", color="grey60")
    }
  } else {
  
    require(ggbeeswarm)
    ### Mutation Burden - Scatter/Dot plot
    ### Works better for larger cohorts
    alpha_val=1
    point_cex=2
    if (nrow(num_var_data) > 200) {
      alpha_val=0.5
    } else if (nrow(num_var_data) > 20) {
      alpha_val=0.8
    }
    burden_plot <- ggplot(num_var_data, aes(x=1, y=mut_burden, text=hoverlabel)) +
      # geom_beeswarm(color="blue3",cex=2,size=5,dodge.width=0.2,priority='density', alpha=alpha_val) +
      geom_quasirandom(color="blue3",width=0.3,size=5,alpha=alpha_val, method="quasirandom", bandwidth=0.1) +
      # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
      #              geom = "crossbar", width = 0.7, color="gray70", size = 0.2) +
      scale_y_log10()+
      theme_linedraw(base_size = 12) +
      ggtitle("Mutation Burden") +
      ylab(y_label_text) + xlab("") +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  
  
  ## Write data to a file for external plotting if desired
  if (!is.null(save_data_to_file)) {
    if (dir.exists(dirname(save_data_to_file))) {
      outdata <- as.data.frame(maf.filtered@variant.classification.summary)
      outdata$total_per_mb <- outdata$total/mb_covered
      outdata$mb_covered <- mb_covered
      print(paste0("Saving plot data to ", save_data_to_file))
      write.table(outdata, file = save_data_to_file, sep="\t", quote=F,row.names = F)
    } else {
      warning("Path for output data file not found! Skipping...")
    }
  }
  
  
  return(burden_plot)
  
}


plot_silent_nonsilent <- function(mymaf, savename=NULL, returndata=F) {
  nonsilent_summary <- mymaf@variant.classification.summary[,c("Tumor_Sample_Barcode","total")]
  nonsilent_summary$type <- "Non-Silent"
  silent_classif_data <- mymaf@maf.silent %>% group_by(Tumor_Sample_Barcode, Variant_Classification) %>% summarise(count=n())
  silent_classif_data <- reshape2::dcast(silent_classif_data,Tumor_Sample_Barcode ~ Variant_Classification, value.var = "count")
  silent_summary <- data.frame(Tumor_Sample_Barcode = silent_classif_data$Tumor_Sample_Barcode,
                               total = rowSums(silent_classif_data[,-1], na.rm=T),
                               type = "Silent"
  )
  
  # browser()
  plotdata <- rbind(nonsilent_summary, silent_summary)
  tots <- plotdata %>% group_by(Tumor_Sample_Barcode) %>% summarise(tot=sum(total))
  plotdata$Tumor_Sample_Barcode <- factor(plotdata$Tumor_Sample_Barcode,
                                          levels=as.character(tots$Tumor_Sample_Barcode)[order(tots$tot,decreasing = T)])
  myplot <- ggplot(plotdata, aes(x=Tumor_Sample_Barcode, y = total, fill=type)) + 
    geom_col() + scale_fill_brewer(palette="Set1") +
    theme_linedraw(base_size = 12) +
    xlab("")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle=30, hjust = 1, size=10),
          axis.ticks.x = element_blank())
  
  if (!is.null(savename)) {
    if (! dir.exists(dirname(savename))) {dir.create(dirname(savename), recursive = T)}
    ggsave(savename,width=6, height=6)
  }
  return_val <- myplot
  if (returndata) {
    return_val <- list(plot=myplot, data=plotdata)
  }
  return(return_val)
}


make_overlap_plot <- function(mymaf, use_silent_mutations=F,
                              summarize_by="gene",
                              plotType=c("ribbon","heatmap"), 
                              savename="overlap_plot.pdf",
                              savewidth=8, saveheight=8) {
  
  require(dplyr)
  mafdata <- mymaf@data
  if (use_silent_mutations) {
    mafdata <- rbind(mafdata, mymaf@maf.silent)
  }
  
  mafdata_by_sample <- split(mafdata, mafdata$Tumor_Sample_Barcode)
  
  if (summarize_by=="gene") {
    id_cols <- c("Chromosome","Start_Position","End_Position","Hugo_Symbol")
  } else if (summarize_by=="mutation") {
    id_cols <- c("Chromosome","Start_Position","End_Position","Hugo_Symbol","HGVSp_Short")
  } else {
    stop("'summarize_by' must be either 'gene' or 'mutation'")
  }
  
  mutations_list <- lapply(mafdata_by_sample, function(currmafdata, mycols) {
    id_string <- apply(currmafdata[,..mycols],1,paste0, collapse="_")
    return(id_string)
  }, id_cols)
  
  pw_combinations <- matrix(0,nrow = length(mutations_list),ncol = length(mutations_list))
  colnames(pw_combinations) <- names(mutations_list)
  rownames(pw_combinations) <- names(mutations_list)
  for ( row_idx in 1:nrow(pw_combinations) ) {
    for (col_idx in 1:ncol(pw_combinations) ) {
      if (!row_idx==col_idx) {
        # overlap_val=sum(mutations_list[[ rownames(pw_combinations)[row_idx] ]]$id_str %in% mutations_list[[colnames(pw_combinations)[col_idx] ]]$id_str)
        overlap_val=length(intersect(mutations_list[[ rownames(pw_combinations)[row_idx] ]],mutations_list[[colnames(pw_combinations)[col_idx] ]]))
        pw_combinations[row_idx,col_idx]=overlap_val
      }
    }
  }
  
  
  hm_data <- pw_combinations
  
  pdf(savename, width=savewidth, height=saveheight)
  if ("heatmap" %in% plotType) {
    require(ComplexHeatmap)
    # pheatmap(hm_data,cluster_rows = T, cluster_cols = T,
    #                clustering_method = "complete",
    #                main = "Hierarchically Clustered")
    myhm <- Heatmap(hm_data,
            cluster_rows = T, cluster_columns = T,
            column_title = paste0("Number of ",ifelse(summarize_by=="gene","mutated genes","variants")," in common between samples (pair-wise)")
            )
    draw(myhm)
  }
  
  if ("ribbon" %in% plotType) {
    require(circlize)
    count_clustering <- hclust(dist(pw_combinations), method = "complete")
    cluster_order <- count_clustering$labels[count_clustering$order]
    
    chordData <- pw_combinations
    grid.col <- rainbow(n=nrow(chordData))
    names(grid.col) <- cluster_order
    
    chordDiagram(chordData,order=cluster_order,grid.col=grid.col,
                 annotationTrack = c("grid","axis"), 
                 preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(chordData))))))
    circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                  facing = "clockwise", niceFacing = TRUE, adj = c(-0.2, 0))
    }, bg.border = NA) # here set bg.border to NA is important
  }
  dev.off()
  invisible(NULL)
}





make_single_ribbon_plot <- function(maf, onco_genes=NULL, save_name=NULL, ribbon_color=NULL, 
                                    pval_high=0.1,  ## All interactions with less than this p-value will be shown
                                    pval_low=0.05,  ## Links with p-value less than this will be highlighted with a dashed border
                                    plot_frac_mut_axis=TRUE,  ## Whether or not to draw a numerical axis on the perimeter
                                    rotate_plot_degrees=0,   ## For custom rotation
                                    shrink_factor=1.3, # Higher = more shrinkage; use to control whitespace (or lack thereof) around figure 
                                    scale_ribbon_to_fracmut=TRUE,  ## Whether or not to scale ribbon widths to their frequency
                                    sig_colors=NULL,   ## Vector of 4 colors for coloring significance
                                    gene_colors=NULL   ## color(s) for gene blocks
) {
  # pval_low <- 0.05
  # browser()
  require(circlize)
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    plot_file <- gsub(".pdf",".interactions.pdf",save_name)
    pdf(file = plot_file,height=5,width=5)
  } else {
    pdf(file = NULL)
  }
  # browser()
  # if (is.null(onco_genes)) {
  #   onco_genes = maf@gene.summary$Hugo_Symbol
  # }

  som_int <-  somaticInteractions(maf = maf, genes=onco_genes, pvalue = c(pval_low, pval_high))
  dev.off()
  # browser()
  cooccur_data <- som_int
  cooccur_data$pair_string <- apply(cooccur_data[,1:2], 1, function(x) {paste0(sort(x), collapse="_")})
  cooccur_data$popfrac <- NA
  cooccur_idx <- cooccur_data$Event=="Co_Occurence"
  mut_excl_idx = cooccur_data$Event=="Mutually_Exclusive"
  
  if (scale_ribbon_to_fracmut) {
    cooccur_data$popfrac1[cooccur_idx] <- unlist(cooccur_data[cooccur_idx,"11"]/as.numeric(maf@summary$summary[3]))
    cooccur_data$popfrac2[cooccur_idx] <- cooccur_data$popfrac1[cooccur_idx]
    cooccur_data$popfrac1[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"10"]/as.numeric(maf@summary$summary[3]))
    cooccur_data$popfrac2[mut_excl_idx] <- unlist(cooccur_data[mut_excl_idx,"01"]/as.numeric(maf@summary$summary[3]))
  } else {
    cooccur_data$popfrac1 <- 1
    cooccur_data$popfrac2 <- 1
  }
  chord_data <- cooccur_data[,c("gene1","gene2","popfrac1","popfrac2","pValue","Event")]
  chord_data[which(is.na(chord_data[,3])),3] <- 0
  
  if (is.null(sig_colors)) {
    sig_colors = RColorBrewer::brewer.pal(5, "BrBG")
    sig_colors <- sig_colors[-3]
  }
  names(sig_colors) <- paste0(rep(c("Co-occurence", "Mutually Exclusive"), each=2), " p-val < ",c(pval_low, pval_high, pval_high, pval_low ))
  color_legend <- Legend(labels=names(sig_colors),
                         legend_gp = gpar(fill = sig_colors, col=sig_colors),background = sig_colors,size=unit(0.08,"npc"),
                         type="points",direction="vertical")
  
  
  chord_data$color_category <- paste0(ifelse(cooccur_data$Event=="Co_Occurence", "Co-occurence", "Mutually Exclusive"),
                                      paste0( " p-val < ", ifelse(cooccur_data$pValue < pval_low, pval_low, pval_high)))
  chord_data$color_val <- sig_colors[chord_data$color_category]
  
  require(RColorBrewer)
  # browser()
  interacting_genes <- unique(unlist(chord_data[,1:2]))
  if (is.null(gene_colors)) {
    gene_colors <- colorRampPalette(brewer.pal(8,"Accent"))(length(interacting_genes))
    # gene_colors <- colorRampPalette(brewer.pal(8,"Dark2"))(length(interacting_genes))
    # gene_colors <- colorRampPalette(brewer.pal(8,"Set1"))(length(interacting_genes))
    # gene_colors <- rainbow((length(interacting_genes)))
  }
  if (length(gene_colors) != length(interacting_genes)) {
    # tmpcolors <- rep("grey90", length(interacting_genes))
    tmpcolors <- rep(gene_colors, length(interacting_genes))
    # tmpcolors[1:length(interacting_genes)] <- gene_colors
    gene_colors <- tmpcolors
  }
  names(gene_colors) <- interacting_genes
  
  if (!is.null(save_name)) {
    if (! dir.exists(dirname(save_name))) {
      dir.create(dirname(save_name))
    }
    pdf(file = save_name,height=7,width=7)
  }

  circos.clear()
  circos.par(canvas.xlim=c(-shrink_factor,shrink_factor),
             canvas.ylim=c(-shrink_factor,shrink_factor),
             start.degree = rotate_plot_degrees,
             message=F)
  # circos.par$message = FALSE
  chordDiagram(chord_data[,1:4],grid.col = gene_colors,
               annotationTrack = c("grid",ifelse(plot_frac_mut_axis, "axis", "")),
               col=chord_data$color_val,
               # transparency = link_alpha, 
               link.lty = 0,
               link.border = "black",
               link.sort = TRUE)
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
                facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0.5))
  }, bg.border = NA) # here set bg.border to NA is important
  draw(color_legend, x = unit(0.05, "npc"), y = unit(0.5, "npc"), just = c("left"))
  # draw(line_legend, x = unit(0.5, "npc"), y = unit(0.97, "npc"), just = c("center"))
  
  if (!is.null(save_name)) {
    dev.off()
  }
  
  
}



detect_maf_genome <- function(maf) {
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


### Cretaes matrix for oncoplot from maf file
### Adapted from maftools: https://github.com/PoisonAlien/maftools/blob/master/R/oncomatrix.R
createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){
  
  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }
  
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)
  
  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      vc = c("")
      names(vc) = 0
      
      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }
  
  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = unique(g))]
  }
  
  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x){
                                x = unique(as.character(x))
                                xad = x[x %in% c('Amp', 'Del')]
                                xvc = x[!x %in% c('Amp', 'Del')]
                                
                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                  # xvc = paste0(xvc, collapse="|")
                                }
                                
                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)
  
  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])
  
  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)
  
  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  
  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  
  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
  
  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  
  
  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId
    
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]
    
    mdf = mdf[, -ncol(mdf)]
    
    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
    
    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
    
    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy
    
    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}

### Cretaes matrix for oncoplot from maf file annotated with ClinVar pathogenicity
### Adapted from maftools: https://github.com/PoisonAlien/maftools/blob/master/R/oncomatrix.R
### Modified to capture ClinVar annotations
createOncoMatrix_CLINSIG = function(m, g = NULL, chatty = TRUE, add_missing = FALSE){
  
  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }
  
  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)
  
  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])
      
      vc = c("")
      names(vc) = 0
      
      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }
  
  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }
  # browser()
  # mafdata <- data.table::data.table(subMaf[,.(Hugo_Symbol, Tumor_Sample_Barcode, CLIN_SIG, Variant_Classification)] %>% 
  #                   mutate(oncomat_label=paste0(ifelse(length(unique(Variant_Classification))>1, "Multi_Hit",unique(Variant_Classification)), 
  #                                               ifelse(grepl("pathogenic",CLIN_SIG),";pathogenic",""))))
  mafdata <- data.table::data.table(
    
    subMaf[,.(Hugo_Symbol, Tumor_Sample_Barcode, CLIN_SIG, Variant_Classification)] %>%
      group_by(Hugo_Symbol, Tumor_Sample_Barcode) %>%
      mutate(oncomat_label=paste0(ifelse(length(unique(as.character(Variant_Classification)))>1, "Multi_Hit",unique(as.character(Variant_Classification))),
                                  # mutate(oncomat_label=paste0(length(unique(as.character(Variant_Classification))),
                                  ifelse(grepl("pathogenic",CLIN_SIG),";Pathogenic",
                                         ifelse(grepl("uncertain",CLIN_SIG), ";VUS",""))))
  )
  oncomat = data.table::dcast(data = mafdata, formula = Hugo_Symbol ~ Tumor_Sample_Barcode, value.var="oncomat_label", fill='', drop=F,fun.aggregate=unique)
  
  
  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])
  
  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)
  
  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)
  
  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }
  
  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
  
  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)
  
  
  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId
    
    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]
    
    mdf = mdf[, -ncol(mdf)]
    
    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix
    
    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort
    
    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy
    
    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]
    
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }
}

make_variant_table <- function(maf.filter, use_syn=F, extra_cols=c()) {
  
  output_data <- maf.filter@data
  if (use_syn) {
    output_data <- rbind(output_data, maf.filter@maf.silent)
  }
  
  output_data$tumor_genotype <- apply(output_data[,c("Tumor_Seq_Allele1","Tumor_Seq_Allele2")], 1, paste, collapse="/")
  output_data$normal_genotype <- apply(output_data[,c("Match_Norm_Seq_Allele1","Match_Norm_Seq_Allele2")], 1, paste, collapse="/")
  
  # pheno_info <- sample_info.exome[match(output_data$Tumor_Sample_Barcode, sample_info.exome$Tumor_Sample_Barcode),]
  # pheno_info <- cbind(pheno_info[,"Tumor_Sample_Barcode"],pheno_info[,-c("Tumor_Sample_Barcode")])
  # pheno_columns <- colnames(pheno_info)
  # names(pheno_columns) <- make.names(pheno_columns, unique = T)
  if (! "tumor_freq" %in% colnames(output_data)) {
    # browser()
    if (all(c("t_depth","t_alt_count")%in% colnames(output_data))) {
      output_data$tumor_freq <- as.numeric(as.character(output_data$t_alt_count))/as.numeric(as.character(output_data$t_depth))
    }
  }
  # output_data <- cbind(output_data,pheno_info)
  cols_for_table <- c("Hugo Symbol" = "Hugo_Symbol",
                      "Sample ID" = "Tumor_Sample_Barcode",
                      "Variant Classification"="Variant_Classification",
                      "Variant Type"="Variant_Type",
                      "Consequence"="Consequence",
                      # pheno_columns,
                      "Chromosome"="Chromosome","Start Position" ="Start_Position","End Position"="End_Position","Strand"="Strand",
                      "Reference Allele"="Reference_Allele",
                      "Tumor Genotype"="tumor_genotype",
                      "Normal Genotype"="normal_genotype",
                      "Known Effects ClinVar"="CLIN_SIG",
                      "Transcript Change"="HGVSc",
                      "Protein Change"="HGVSp_Short",
                      "Normal Depth"="n_depth",
                      "Normal Ref Depth"="n_ref_count",
                      "Normal Alt Depth"="n_alt_count",
                      "Tumor Depth"="t_depth",
                      "Tumor Ref Depth"="t_ref_count",
                      "Tumor Alt Depth"="t_alt_count",
                      "Tumor Alt Frequency"="tumor_freq",
                      "Existing Annotation"="Existing_variation",
                      "gnomAD Frequency"="gnomAD_AF",
                      "ExAC Frequency"="ExAC_AF",
                      "1000Genomes Frequency"="AF",
                      "Effect Prediction - SIFT"="SIFT",
                      "Effect Prediction - PolyPhen"="PolyPhen"
  )
  # browser()
  cols_for_table <- c(cols_for_table, extra_cols)
  # variant_info <- as.data.frame(output_data)[,cols_for_table]
  norm_info_cols <- grep("^n_",cols_for_table, value=T)
  # mydat <- apply(output_data[,..norm_info_cols],2,function(x){as.numeric(x)})
  if (sum(rowSums(apply(output_data[,..norm_info_cols],2,function(x){as.numeric(x)}), na.rm=T), na.rm = T)==0) {
    cols_for_table <- cols_for_table[!cols_for_table %in% norm_info_cols]
  }
  output_cols <- colnames(output_data)[match(cols_for_table, colnames(output_data), nomatch=0)]
  not_output <- cols_for_table[!cols_for_table %in% output_cols]
  if (length(not_output) > 0) {
    print("Not outputting these columsn: ")
    print(not_output)
  }
  variant_info <- as.data.frame(output_data)[,output_cols]
  colnames(variant_info) <- names(cols_for_table)[match(colnames(variant_info),cols_for_table)]
  return(variant_info)
}


compute_exome_coverage <- function(targets_bed_file, out_file=NULL) {
  ##### This function will read the target regions BED file and
  #####  compute the sum of the lengths of the regions
  require(GenomicRanges)
  
  ## This bit will only read in the first three columns
  num_fields <- max(count.fields(targets_bed_file, sep = "\t"))
  my_classes <- c("character","integer","integer", rep("NULL", num_fields-3))
  
  ## Read the BED file as a table
  bed_data <- read.table(targets_bed_file, sep="\t",colClasses = my_classes, 
                         stringsAsFactors = F)
  colnames(bed_data) <- c("chr","start","end")
  
  ## Convert to a GenomicRanges object
  bed.gr <- makeGRangesFromDataFrame(bed_data)
  
  ## Collapse any overlapping features
  bed.gr.collapse <- reduce(bed.gr)
  
  ## Sum up the widths of each range
  total_exome_coverage = sum(width(bed.gr.collapse))
  
  if (! is.null(out_file)) {
    ## Write to a file
    write.table(total_exome_coverage,file = out_file, col.names =F, row.names = F)
  } else {
    return(total_exome_coverage)
  }
}




### Define colors for mutation types
mutation_colors <- c(Nonsense_Mutation="#ad7aff",Missense_Mutation="#377EB8",Frame_Shift_Del="#4DAF4A",
         In_Frame_Ins="#ff008c",Splice_Site="#FF7F00",Multi_Hit="#FFFF33",Frame_Shift_Ins="#A65628",
         In_Frame_Del="#f781bf",Translation_Start_Site="#400085",Nonstop_Mutation="#b68dfc",
         Amp="green2",Del="darkred",
         no_variants="#d6d6d6", Pathogenic="black",VUS="grey50")
names(mutation_colors) <- gsub("_"," ",names(mutation_colors))
### List defining functions for color and shape of cells in oncoplot
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # "0" = function(x, y, w, h) {
  #   grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
  #             gp = gpar(fill = "#CCCCCC", col = NA))
  # },
  "Nonsense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Nonsense Mutation"], col = NA))
  },
  "Missense Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Missense Mutation"], col = NA))
  },
  "Frame Shift Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["Frame Shift Del"], col = NA))
  },
  "In Frame Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["In Frame Ins"], col = NA))
  },
  "Splice Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Splice Site"], col = NA))
  },
  "Multi Hit" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Multi Hit"], col = NA))
  },
  "Frame Shift Ins" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Frame Shift Ins"], col = NA))
  },
  "In Frame Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = mutation_colors["In Frame Del"], col = NA))
  },
  "Nonstop Mutation" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Nonstop Mutation"], col = NA))
  },
  "Translation Start Site" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = mutation_colors["Translation Start Site"], col = NA))
  },
  "Amp" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.25, "mm"), h-unit(0.25, "mm"),
              gp = gpar(fill = mutation_colors["Amp"], col = NA))
  },
  "Del" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.25, "mm"), h-unit(0.25, "mm"),
              gp = gpar(fill = mutation_colors["Del"], col = NA))
  },
  "no variants" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              # gp = gpar(fill = "#e0e0e0", col = NA))
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  "Pathogenic" = function(x, y, w, h) {
    # grid.points(x, y, pch = 18, size=w, gp=gpar(col=col["pathogenic"]))
    # grid.rect(x, y, w*0.7, h*0.2,
    #           gp = gpar(fill = col["pathogenic"], col = NA))
    # grid.rect(x, y, w*0.1, h*0.7,
    #           gp = gpar(fill = col["pathogenic"], col = NA))
    grid.rect(x, y, w*0.8, h*0.8,
              gp = gpar(col = mutation_colors["Pathogenic"], fill = NA, lwd=5))
  },
  "VUS" = function(x, y, w, h) {
    # grid.points(x, y, pch = 3, size=w,gp=gpar(col=col["VUS"], lwd=3))
    # grid.rect(x, y, w*0.2, h-unit(0.5, "mm"),
    #           gp = gpar(fill = col["VUS"], col = NA))
    grid.rect(x, y, w*0.8, h*0.8,
              gp = gpar(col = mutation_colors["VUS"], fill = NA, lwd=5))
  }
)


