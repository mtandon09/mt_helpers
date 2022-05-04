get_sig_genes <- function(limma_df, pval = 0.05, topn = NA, name_with_symbol = FALSE, 
                          direction = "all", gene_col = NULL, fc_col = "logFC", sig_col = "adj.P.Val") {
  
  if (is.null(gene_col)) {
    limma_df <- cbind(rowid=rownames(limma_df), limma_df)
    gene_col="rowid"
  }
  
  if (! all(c(gene_col, fc_col, sig_col) %in% colnames(limma_df))){
    stop("Couldn't find columns named: ", paste0(setdiff(c(gene_col, fc_col, sig_col), colnames(limma_df)), collapse=","))
  }
  
  limma_df <- limma_df[order(limma_df[,sig_col]), ]
  # browser()
  if (name_with_symbol) {
    limma_df[,gene_col] <- make.names(sapply(strsplit(limma_df[,gene_col],"\\|"),"[[",2), unique = TRUE)
  }
  if (direction == "up") {
    limma_df <- limma_df[limma_df[, fc_col] > 0, ]
  } else if (direction == "down") {
    limma_df <- limma_df[limma_df[, fc_col] < 0, ]
  } else if (direction %in% c("all","any")) {
    # limma_df <- limma_df
  } else {
    stop("Unknown 'direction' supplied.")
  }
  # browser()
  limma_df <- limma_df[order(limma_df[,sig_col], decreasing = F), ]
  
  return_vector <- as.vector(limma_df[limma_df[,sig_col] < pval, fc_col])
  names(return_vector) <- limma_df[limma_df[,sig_col] < pval, gene_col]
  
  if (!is.na(topn)) {
    if (length(return_vector) > topn) {
      return_vector <- return_vector[1:topn]
    }
  }
  
  return(return_vector)
  
}

makeVolcano <- function(myresults,labelTop=TRUE,nLabel=20,pval_cutoff = 0.05, pvalType = "padj", geneLabel = NULL,
                        saveplot=FALSE,savename="volcano.png",title_suffix="",colorpalette = "Set1") {
  
  require(reshape2)
  require(ggplot2)
  require(ggrepel)
  require(RColorBrewer)
  
  myresults <- myresults[order(myresults[,pvalType], decreasing = FALSE),]
  
  if (! pvalType %in% colnames(myresults)) {
    stop(paste0("pvalType must be one of columns provided"))
  }
  
  if (grepl("adj",pvalType)) {
    ylabel <- "-log(corrected p-value)"
  } else {
    ylabel <- "-log(p-value)"
  }
  
  
  
  # browser()
  fc_col <- grep("^log",tolower(colnames(myresults)))[1]
  if (is.na(fc_col)) {
    fc_col <- grep("fc",tolower(colnames(myresults)))[1]
    if (is.na(fc_col)) {
      stop("Can't find fold-change column")
    } else {
      print("Found fold-change column, converting to log scale...")
      logfc <- log2(ifelse(myresults[,fc_col] < 0, myresults[,fc_col], 1/abs(myresuls[,fc_col])))
    }
  } else {
    print("Found fold-change column, already in log scale...")
    logfc <- myresults[,fc_col]
  }
  
  if (is.null(geneLabel)) {
    gene_names <- rownames(myresults)
  } else {
    gene_names <- myresults[,geneLabel]
  }
  
  myresults.anno <- data.frame(gene=gene_names,
                               logFC=logfc,
                               pval=-log10(myresults[,pvalType]),
                               sig <- myresults[,pvalType] < pval_cutoff)
  # myresults.anno <- myresults.anno[order(myresults.anno$pval),]
  # fillColors <- colorRampPalette(brewer.pal(8,colorpalette))(2)
  fillColors <- brewer.pal(3,colorpalette)[1:2]
  names(fillColors) <- unique(myresults.anno$sig)
  
  
  myplot <- ggplot(myresults.anno,aes(x=logFC,y=pval)) +
    geom_point(aes(color=sig),alpha=0.4) +
    xlab("log(Fold-Change)") + ylab(ylabel) + ggtitle(paste0("Volcano Plot: ",title_suffix)) +
    scale_color_manual(values=fillColors) +
    theme_bw()
  if (labelTop) {
    myresults.label <- myresults.anno[1:nLabel,c("logFC","pval","gene")]
    head(myresults.label)
    myplot <- myplot + geom_label_repel(data=myresults.label,aes(x=logFC,y=pval,label=gene))
  }
  if (saveplot) {
    ggsave(filename=savename,plot=myplot, width=7, height=6)
  } else {
    print(myplot)
  }
  
}

# makePCAplot <- function(expr_vals,pheno.data,topNgenes=NULL, sampleLabel = FALSE,
#                         whichPCs=c(1,2),fillVar=NULL, shapeVar=NULL, fillColors=NULL,
#                         return_plotly=F,
#                         drawEllipses=FALSE,saveplot=FALSE,savefilename,savedata=FALSE,
#                         point_size=5, point_alpha=0.9, point_outline_color=NA, label_size=0.8) {
#   require(ggplot2)
#   require(RColorBrewer)
#   require(ggrepel)
#   if (missing(topNgenes)) {
#     topNgenes <- nrow(expr_vals)
#   }
#   
#   if (missing(fillVar)) {
#     if ("condition" %in% colnames(pheno.data)) {
#       fillVar = "condition"
#     } else {
#       fillVar = colnames(pheno.data)[1]
#     }
#     warning(paste0("No coloring variable specified. Using ",fillVar," to color samples..."))
#   }
#   
#   # pheno.data$fillVar <- pheno.data[,fillVar] 
#   if (is.null(fillColors)) {
#     # browser()
#     mygroups <- pheno.data[,fillVar]
#     numgroups <- length(unique(mygroups))
#     colorpalette <- "Paired"
#     fillColors <- colorRampPalette(brewer.pal(numgroups,colorpalette))(numgroups)
#     names(fillColors) <- unique(mygroups)
#   }
#   
#   if (is.null(shapeVar)) {
#     # pheno.data$shapeVar = rep(15,nrow(pheno.data))
#     # shapeVar = fillVar
#   }
#   
#   # pheno.data$shapeVar = pheno.data[,shapeVar]
#   n_shapes <- length(unique(pheno.data[,shapeVar]))
#   if (n_shapes > 5) {
#     warning("Can't plot more than 5 different values for 'shapeVar'")
#     pheno.data$shapeVar = "datapoint"
#     shapeVar=NULL
#   }
#   
#   pcadata <- expr_vals
#   gr <- pheno.data[match(colnames(pcadata),rownames(pheno.data)),]
#   # rownames(gr) <- colnames(pcadata)
#   # colnames(gr)[colnames(gr) == "sample_id"] <- "id"
#   
#   ## Compute PCA and variance contributions
#   no_var_rows <- which(apply(pcadata, 1, var)==0)
#   if (length(no_var_rows) > 0) {
#     print("Some rows have zero variance; turning off unit scaling of PCs...")
#     scale_prcmp = FALSE
#   } else {
#     scale_prcmp = TRUE
#   }
#   pca <- prcomp(t(pcadata), scale=scale_prcmp)
#   eigs <- pca$sdev^2
#   var_explained <- paste(sprintf(eigs/sum(eigs)*100,fmt = '%.2f'),"%")
#   
#   # browser()
#   ## Plot PCA
#   library(ggplot2)
#   library(reshape2)
#   library(RColorBrewer)
#   plotdata <- data.frame(cbind(pca$x[,whichPCs],gr))
#   colnames(plotdata)[1:2] <- paste0("PC",whichPCs)
#   plotdata.melt <- melt(plotdata, id.vars = colnames(gr))
# 
#   pcaplot <- ggplot(plotdata, aes_string(x = paste0("PC",whichPCs[1]), y = paste0("PC",whichPCs[2]))) +
#     # geom_point(aes_string(shape=shapeVar), size = point_size, alpha=point_alpha, color=NA) + 
#     # geom_point(aes_string(shape=shapeVar, fill=fillVar),size = point_size, alpha=point_alpha, color="black") + 
#     geom_point(aes_string(shape=shapeVar, fill=fillVar),size = point_size, alpha=point_alpha, color="black") + 
#     # scale_color_brewer(palette = "Paired") +
#     scale_fill_manual(values = fillColors) +
#     # scale_color_manual(values = fillColors) +
#     # scale_shape_manual(values = all_shapes) +
#     xlab(paste0("PC1 (Variance explained: ", var_explained[1],")")) + 
#     ylab(paste0("PC2 (Variance explained: ", var_explained[2],")")) +
#     theme_linedraw(base_size = 14)
#   if (! is.null(shapeVar)) {
#     all_shapes <- c(21:25)
#     names(all_shapes) <- unique(as.character(pheno.data[,shapeVar]))
#     all_shapes <- all_shapes[!is.na(names(all_shapes))]
#     
#     pcaplot <- pcaplot + scale_shape_manual(values = all_shapes)
#   }
#   # pcaplot
#   if (sampleLabel) {
#     pcaplot <- pcaplot + geom_label_repel(aes(label = sample_id), size=label_size)
#   }
#   if (drawEllipses) {
#     pcaplot <- pcaplot + stat_ellipse(type="norm")
#   }
#   
#   # print(pcaplot)
#   if (saveplot || !missing(savefilename)) {
#     if (missing(savefilename)) {
#       savefilename = "pcaplot.pdf"
#     }
#     ggsave(savefilename,plot=pcaplot,units="in",width=10,height=7)
#   } else {
#     # print(pcaplot)
#   }
#   # browser()
#   if (savedata) {
#     write.table(pcadata,file=paste0(savefilename,".data.txt"),row.names=FALSE,col.names = TRUE,sep="\t")
#   }
#   return(pcaplot)
# } 


makePCAplot_MT <- function(expr_vals,pheno.data=NULL,topNgenes=NULL, sampleLabel = FALSE, labelColumn=NULL,
                           biplot=FALSE, n_loadings=10,
                           make_plotly = F,
                           whichPCs=c(1,2),
                           pointsize=6, auto_adjust_pointsize=T, bordersize=2,
                           fillVar=NULL, fillColors=NULL, 
                           shapeVar=NULL, shapeVals=NULL, 
                           drawEllipses=FALSE,saveplot=FALSE,savefilename,savedata=FALSE) {
  require(ggplot2)
  require(RColorBrewer)
  require(ggrepel)
  
  if (missing(pheno.data) || is.null(pheno.data)) {
    pheno.data <- data.frame(sampleID=colnames(expr_vals),
                             fillVar=NA,
                             shapeVar=NA,
                             row.names = colnames(expr_vals))
  }
  if (missing(topNgenes) || is.null(topNgenes)) {
    topNgenes <- nrow(expr_vals)
  }
  myind <- order(apply(expr_vals, 1, var), decreasing = T)
  myind <- myind[1:min(c(topNgenes, nrow(expr_vals)))]
  expr_vals <- expr_vals[myind,]
  
  # browser()
  # if (missing(fillVar) || is.null(fillVar)) {
  #   if ("condition" %in% colnames(pheno.data)) {
  #     fillVar = "condition"
  #   }
  #   # warning(paste0("No coloring variable specified. Using ",fillVar," to color samples..."))
  # }
  # if (! fillVar %in% colnames(pheno.data)) {
  #   print(paste0(fillVar," not found in phenotype data.  Using this instead: ",colnames(pheno.data)[1]))
  #   fillVar = colnames(pheno.data)[1]
  # }
  if (is.null(fillVar)) {
    fillVar="fillVar"
  }
  continuous_color=F
  if (missing(fillColors) || is.null(fillColors)) {
    colordata <- pheno.data[,fillVar]
    if (is.numeric(colordata)) {
      continuous_color=T
      colorpalette <- "PiYG"
      if (min(colordata, na.rm=T) >= 0) {
        colorpalette <- "Blues"
      }
      fillColors <- setNames(brewer.pal(9,colorpalette)[c(1,9)],c("low","high"))
      
    } else {
      numgroups <- length(unique(colordata))
      if (numgroups < 13) {
        colorpalette <- "Paired"
        fillColors <- brewer.pal(12,colorpalette)[1:numgroups]
      } else {
        fillColors <- rainbow(numgroups)
      }
      names(fillColors) <- unique(colordata)
      
    }
  }
  
  if (missing(shapeVar) || is.null(shapeVar)) {
    # pheno.data$shapeVar = rep(15,nrow(pheno.data))
    pheno.data$shapeVar = pheno.data[,fillVar]
  } else {
    pheno.data$shapeVar = pheno.data[,shapeVar]
  }
  if (missing(shapeVals) || is.null(shapeVals)) {
    # shapeVals=c(0:14)[1:length(unique(pheno.data$shapeVar))]
    shapeVals=c(21:26)[1:length(unique(pheno.data$shapeVar))]
  }
  if (is.null(names(shapeVals))) {
    names(shapeVals) <- sort(unique(pheno.data$shapeVar))
  }
  
  if (missing(labelColumn) || is.null(labelColumn)) {
    pheno.data$sample_id <- rownames(pheno.data)
  } else {
    pheno.data$sample_id <- pheno.data[,labelColumn]
  }
  rownames(pheno.data) <- pheno.data$sample_id
  
  pcadata <- expr_vals
  gr <- pheno.data[match(colnames(pcadata),rownames(pheno.data)),]
  
  ## Compute PCA and variance contributions
  # no_var_rows <- which(apply(pcadata, 1, var)==0)
  # if (length(no_var_rows) > 0) {
  no_var_rows <- unlist(apply(pcadata, 1, var)==0)
  if (sum(no_var_rows) > 0) {
    # browser()
    print("Some rows have zero variance; turning off unit scaling of PCs...")
    scale_prcmp = FALSE
  } else {
    scale_prcmp = TRUE
  }
  pca <- prcomp(t(pcadata), scale=scale_prcmp)
  eigs <- pca$sdev^2
  var_explained <- paste(sprintf(eigs/sum(eigs)*100,fmt = '%.2f'),"%")
  
  if (biplot) {
    # pcaplot <- ggbiplot_MT(pca)
    # browser()
    if (missing(n_loadings) || as.integer(n_loadings) > nrow(pca$rotation)) {
      # loadingsInd <- 1:as.integer(n_loadings)
      n_loadings=as.integer(n_loadings)
    }
    loadingsInd <- 1:as.integer(n_loadings)
    nobs.factor <- sqrt(nrow(pca$x) - 1)
    d <- pca$sdev
    u <- sweep(pca$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pca$rotation[loadingsInd,]
    var.scale=as.integer(scale_prcmp)
    
    whichPCs <- pmin(whichPCs, ncol(pcadata))
    
    v <- sweep(v, 2, d^var.scale, FUN = "*")
    # df.u <- as.data.frame(sweep(u[, whichPCs], 2, d[whichPCs]^as.integer(scale_prcmp), 
    #                             FUN = "*"))
    # df.u <- df.u * nobs.factor
    # 
    df.v <- as.data.frame(v[, whichPCs])
    colnames(df.v) <- c("xvar", "yvar")
    
    circle.prob = 0.69
    varname.adjust = 1.5
    
    # r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans((pcadata*nobs.factor)^2))^(1/4)
    r <- sqrt(qchisq(circle.prob, df = 2)) * 
      ifelse(is.infinite(max(prod(colMeans((pcadata*nobs.factor)^2)))),
             .Machine$integer.max,
             max(prod(colMeans((pcadata*nobs.factor)^2))))^(1/4)
    v.scale <- rowSums(v^2)
    df.v <- r * df.v/sqrt(max(v.scale))
    df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
    df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
    
    df.v$varname <- rownames(v)
    
    
  }
  
  # browser()
  ## Plot PCA
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  plotdata <- data.frame(cbind(pca$x[,whichPCs],gr))
  colnames(plotdata)[1:2] <- paste0("PC",whichPCs)
  
  if (auto_adjust_pointsize) {
    npoints <- nrow(plotdata)
    pointsize <- ifelse(npoints > 1e4, 0.5,
                        ifelse(npoints > 1000, 1, 
                               ifelse(npoints > 100, 2, 
                                      ifelse(npoints > 10, 4, pointsize))))
  }
  pcaplot <- ggplot(plotdata, aes_string(x = paste0("PC",whichPCs[1]), y = paste0("PC",whichPCs[2]))) +
    geom_point(aes_string(shape=shapeVar, fill = fillVar),color="grey10",size = pointsize, alpha=0.6, stroke = bordersize) + 
    # guides(color = guide_legend(override.aes = list(size=3, alpha=1, shape=15))) +
    # scale_color_brewer(palette = "Paired") +
    # scale_color_manual(values = fillColors) +
    scale_shape_manual(values = shapeVals) +
    xlab(paste0("PC1 (Variance explained: ", var_explained[1],")")) + 
    ylab(paste0("PC2 (Variance explained: ", var_explained[2],")")) +
    ggtitle(paste0("PCA plot (n=",ncol(pcadata),", ",topNgenes," genes)")) +
    theme_linedraw(base_size = 14)
  
  if (continuous_color) {
    # pcaplot <- pcaplot + scale_color_gradient(low=fillColors[1], high=fillColors[2])
    pcaplot <- pcaplot + scale_fill_gradient(low=fillColors[1], high=fillColors[2])
  } else {
    # pcaplot <- pcaplot + scale_color_manual(values = fillColors) +
      # guides(color = guide_legend(override.aes = list(size=3, alpha=1, shape=15)))
    pcaplot <- pcaplot + scale_fill_manual(values = fillColors) +
      guides(fill = guide_legend(override.aes = list(size=8, alpha=1, shape=22)))
  }
  
  if (sampleLabel && !make_plotly) {
    pcaplot <- pcaplot + geom_label_repel(aes(label = sample_id))
  }
  if (drawEllipses) {
    pcaplot <- pcaplot + stat_ellipse(type="norm")
  }
  
  if (biplot) {
    pcaplot <- pcaplot + geom_segment(data = df.v, 
                                      aes(x = 0, y = 0, xend = xvar, yend = yvar), 
                                      size=1.3,
                                      arrow = arrow(length = unit(0.05,"npc"), type="open"),
                                      color = "red3")
    if (!make_plotly) {
      varname.size = 3
      pcaplot <- pcaplot + geom_text_repel(data = df.v, 
                                           aes(label = varname, x = xvar, y = yvar, angle = 0), 
                                           color = "red3", 
                                           size = varname.size, 
                                           segment.color="red3",
                                           segment.size=0.2,force=50,segment.alpha=0.3)
    }
  }
  
  if (make_plotly) {
    require(plotly)
    require(htmlwidgets)
    pcaplot <- ggplotly(pcaplot)
  }
  
  if (saveplot || !missing(savefilename)) {
    if (missing(savefilename)) {
      savefilename = "pcaplot.pdf"
    }
    if (make_plotly) {
      # savefilename <- gsub("/","//",savefilename)
      # saveWidget(pcaplot, paste0(gsub("\\.html$","",savefilename), ".html"))
      saveWidgetFix(pcaplot, paste0(gsub("\\.html$","",savefilename), ".html"))
    } else {
      ggsave(savefilename,plot=pcaplot,units="in",width=10,height=7)
    }
  } else {
    print(pcaplot)
  }
  # browser()
  # if (savedata) {
  #   #write.table(pcadata,file=paste0(savefilename,".data.txt"),row.names=FALSE,col.names = TRUE,sep="\t")
  #   write.table(plotdata,file=paste0(savefilename,".plot.txt"),row.names=FALSE,col.names = TRUE,sep="\t")
  # }
  
  invisible(pcaplot)
} 


saveWidgetFix <- function (widget,file,...) {
  ########### FROM: https://github.com/ramnathv/htmlwidgets/issues/299 ###########
  ## A wrapper to saveWidget which compensates for arguable BUG in
  ## saveWidget which requires `file` to be in current working
  ## directory.
  wd<-getwd()
  on.exit(setwd(wd))
  outDir<-dirname(file)
  file<-basename(file)
  setwd(outDir);
  saveWidget(widget,file=file,...)
}


# ggbiplot_MT <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, n_loadings=10,
#                          obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
#                          ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
#                          alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
#                          varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
#                          pointsize = 5,
#                          ...) 
# {
#   library(ggplot2)
#   library(plyr)
#   library(scales)
#   library(grid)
#   library(ggrepel)
#   stopifnot(length(choices) == 2)
#   # browser()
#   
#   
#   # if ((max(loadingsInd) > nrow(pcobj$rotation)) || (missing(loadingsInd))) {
#   if (missing(n_loadings) || as.integer(n_loadings) > nrow(pcobj$rotation)) {
#     # loadingsInd <- 1:nrow(pcobj$rotation)
#     loadingsInd <- 1:n_loadings
#   }
#   
#   if (inherits(pcobj, "prcomp")) {
#     nobs.factor <- sqrt(nrow(pcobj$x) - 1)
#     d <- pcobj$sdev
#     u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
#     # browser()
#     v <- pcobj$rotation[loadingsInd,]
#   } else if (inherits(pcobj, "princomp")) {
#     nobs.factor <- sqrt(pcobj$n.obs)
#     d <- pcobj$sdev
#     u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
#     v <- pcobj$loadings
#   } else if (inherits(pcobj, "PCA")) {
#     nobs.factor <- sqrt(nrow(pcobj$call$X))
#     d <- unlist(sqrt(pcobj$eig)[1])
#     u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
#     v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
#                                                   1]), FUN = "/")
#   } else if (inherits(pcobj, "lda")) {
#     nobs.factor <- sqrt(pcobj$N)
#     d <- pcobj$svd
#     u <- predict(pcobj)$x/nobs.factor
#     v <- pcobj$scaling
#     d.total <- sum(d^2)
#   } else {
#     stop("Expected a object of class prcomp, princomp, PCA, or lda")
#   }
#   
#   choices <- pmin(choices, ncol(u))
#   df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
#                               FUN = "*"))
#   v <- sweep(v, 2, d^var.scale, FUN = "*")
#   df.v <- as.data.frame(v[, choices])
#   names(df.u) <- c("xvar", "yvar")
#   names(df.v) <- names(df.u)
#   if (pc.biplot) {
#     df.u <- df.u * nobs.factor
#   }
#   r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
#   v.scale <- rowSums(v^2)
#   df.v <- r * df.v/sqrt(max(v.scale))
#   if (obs.scale == 0) {
#     u.axis.labs <- paste("standardized PC", choices, sep = "")
#   }
#   else {
#     u.axis.labs <- paste("PC", choices, sep = "")
#   }
#   
#   u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
#                                             100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
#   if (!is.null(labels)) {
#     df.u$labels <- labels
#   }
#   if (!is.null(groups)) {
#     df.u$groups <- groups
#   }
#   if (varname.abbrev) {
#     df.v$varname <- abbreviate(rownames(v))
#   } else {
#     df.v$varname <- rownames(v)
#   }
#   df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
#   df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
#   g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
#     ylab(u.axis.labs[2]) + coord_equal()
#   if (var.axes) {
#     if (circle) {
#       theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
#                                                 length = 50))
#       circle <- data.frame(xvar = r * cos(theta), yvar = r * 
#                              sin(theta))
#       g <- g + geom_path(data = circle, color = muted("white"), 
#                          size = 1/2, alpha = 1/3)
#     }
#     g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
#                                            xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
#                                                                                                   "picas")), color = muted("red"))
#   }
#   if (!is.null(df.u$labels)) {
#     if (!is.null(df.u$groups)) {
#       # g <- g + geom_text(aes(label = labels, color = groups), 
#       #                    size = labels.size)
#       g <- g + geom_text_repel(aes(label = labels, color = groups), 
#                                size = labels.size)
#     } else {
#       # g <- g + geom_text(aes(label = labels), size = labels.size)
#       g <- g + geom_text_repel(aes(label = labels), size = labels.size)
#     }
#   } else {
#     if (!is.null(df.u$groups)) {
#       g <- g + geom_point(aes(fill = groups), alpha = alpha, size=pointsize,shape=21)
#     } else {
#       g <- g + geom_point(alpha = alpha)
#     }
#   }
#   if (!is.null(df.u$groups) && ellipse) {
#     theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
#     circle <- cbind(cos(theta), sin(theta))
#     ell <- ddply(df.u, "groups", function(x) {
#       if (nrow(x) <= 2) {
#         return(NULL)
#       }
#       sigma <- var(cbind(x$xvar, x$yvar))
#       mu <- c(mean(x$xvar), mean(x$yvar))
#       ed <- sqrt(qchisq(ellipse.prob, df = 2))
#       data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
#                        mu, FUN = "+"), groups = x$groups[1])
#     })
#     names(ell)[1:2] <- c("xvar", "yvar")
#     g <- g + geom_path(data = ell, aes(color = groups, group = groups))
#   }
#   if (var.axes) {
#     # g <- g + geom_text(data = df.v, aes(label = varname, 
#     #                                     x = xvar, y = yvar, angle = angle, hjust = hjust), 
#     #                    color = "darkred", size = varname.size)
#     g <- g + geom_text_repel(data = df.v, aes(label = varname, 
#                                               x = xvar, y = yvar, angle = 0), 
#                              color = "black", size = varname.size,segment.color="grey10",segment.size=0.2,force=20,segment.alpha=0.3)
#   }
#   return(g)
# }





makeVenn <- function(resultsdf,dds.norm,pval_cutoff=0.01,fc_cutoff=0.05,
                     save_name="venn.tiff",write_overlaps=FALSE) {
  require(VennDiagram)
  source("scripts/helper_functions.R")

  allresults.sig <- allresults[complete.cases(resultsdf[,grep("padj",colnames(resultsdf))]),]

  if (sign(fc_cutoff) > 0) {
    venngenes <- list(Cluster1=as.character(allresults.sig$gene[allresults.sig$log2FC_Cluster_1_vs_0 > fc_cutoff & allresults.sig$padj_Cluster_1_vs_0 < pval_cutoff]),
                      Cluster2=as.character(allresults.sig$gene[allresults.sig$log2FC_Cluster_2_vs_0 > fc_cutoff & allresults.sig$padj_Cluster_2_vs_0 < pval_cutoff]),
                      Cluster3=as.character(allresults.sig$gene[allresults.sig$log2FC_Cluster_3_vs_0 > fc_cutoff & allresults.sig$padj_Cluster_3_vs_0 < pval_cutoff]))
    fc_char = ">"
  } else if (sign(fc_cutoff) < 0) {
    venngenes <- list(Cluster1=as.character(allresults.sig$gene[allresults.sig$log2FC_Cluster_1_vs_0 < fc_cutoff & allresults.sig$padj_Cluster_1_vs_0 < pval_cutoff]),
                      Cluster2=as.character(allresults.sig$gene[allresults.sig$log2FC_Cluster_2_vs_0 < fc_cutoff & allresults.sig$padj_Cluster_2_vs_0 < pval_cutoff]),
                      Cluster3=as.character(allresults.sig$gene[allresults.sig$log2FC_Cluster_3_vs_0 < fc_cutoff & allresults.sig$padj_Cluster_3_vs_0 < pval_cutoff]))
    fc_char = "<"
  } else {
    stop("FoldChange cutoff must be non-zero")
  }
  # mycolor = colorRampPalette(rev(brewer.pal(n = length(venngenes), name ="Set1")))(length(venngenes))
  mycolor = anno_col[["Cluster"]][-1]
  if (missing(save_name)) {
      save_name = paste0(resultsdir,"/venn-cluster-gene-lists_p_",pval_cutoff,"_fc_",fc_cutoff,".tiff")
  }
  venn.diagram(venngenes,filename=save_name,
               main=paste0(length(unique(unlist(venngenes)))," Genes with log2FC ",fc_char," ",fc_cutoff,"and p < ",pval_cutoff,"\nin each cluster"),
               imagetype = "tiff",height=5,width=5,units="in",
               cat.col=mycolor,fill=mycolor,cat.cex=1.1,cex=1.3,
               cat.dist=c(0.05,0.05,0.05),
               # cat.pos=c(337,20,345,15),
               margin=0.05)
  
  if (write_overlaps) {
    outfile = paste0(resultsdir,"/venn-cluster-overlaps_p_",pval_cutoff,"_fc_",fc_cutoff,".txt")
    write_set_overlaps(venngenes,outfile_name=outfile)
  }
  junk <- dir(path=resultsdir, pattern=".log$",full.names=TRUE)
  file.remove(junk)
  return(venngenes)
}

st.err <- function(x) {
  sd(x)/sqrt(length(x))
}

writeNormalizedExprs <- function(mydds,normType="rlog",phenoFile = TRUE, resultsdir = ".", prefix="rnaseq") {
  require(DESeq2)
  if (normType == "rlog") {
    dds.norm <- rlog(mydds)
  } else if (normType == "vst") {
    dds.norm <- varianceStabilizingTransformation(mydds)
  } else if (normType == "raw") {
    dds.norm <- mydds
  } else {
    stop("Invalid normalization type. Must be one of c(rlog, vst, raw).")
  }
  
  exprsFile <- paste0(resultsdir,"/",prefix,".",normType,".txt")
  
  if (normType == "raw") {
    myexprs <- data.frame(cbind(gene=rownames(dds.norm),counts(dds.norm)))
  } else {
    myexprs <- data.frame(cbind(gene=rownames(dds.norm),assay(dds.norm)))
  }
  
  write.table(myexprs,file=exprsFile,sep="\t",row.names=FALSE,quote=FALSE)
  
  if (phenoFile) {
    myexprs.pheno <- data.frame(sampleID=rownames(colData(dds.norm)),colData(dds.norm))
    exprsFile.pheno <- paste0(resultsdir,"/",prefix,".pheno.txt")
    write.table(myexprs.pheno,file=exprsFile.pheno,sep="\t",row.names=FALSE,quote=FALSE)
  }
}


ggbiplot_MT <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                         obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                         ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                         alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                         varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                         loadingsInd, pointsize = 5,
                         ...) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  library(ggrepel)
  stopifnot(length(choices) == 2)
  # browser()
  if ((max(loadingsInd) > nrow(pcobj$rotation)) || (missing(loadingsInd))) {
    loadingsInd <- 1:nrow(pcobj$rotation)
  }
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    # browser()
    v <- pcobj$rotation[loadingsInd,]
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                  "picas")), color = muted("red"))
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      # g <- g + geom_text(aes(label = labels, color = groups), 
      #                    size = labels.size)
      g <- g + geom_text_repel(aes(label = labels, color = groups), 
                               size = labels.size)
    }
    else {
      # g <- g + geom_text(aes(label = labels), size = labels.size)
      g <- g + geom_text_repel(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(fill = groups), alpha = alpha, size=pointsize,shape=21)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    # g <- g + geom_text(data = df.v, aes(label = varname, 
    #                                     x = xvar, y = yvar, angle = angle, hjust = hjust), 
    #                    color = "darkred", size = varname.size)
    g <- g + geom_text_repel(data = df.v, aes(label = varname, 
                                              x = xvar, y = yvar, angle = 0), 
                             color = "black", size = varname.size,segment.color="grey10",segment.size=0.2,force=20,segment.alpha=0.3)
  }
  return(g)
}
