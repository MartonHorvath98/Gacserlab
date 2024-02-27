make_deseq <- function(matrix, coldata, design){
  .design <- reformulate(design)
  deseq <- DESeqDataSetFromMatrix(matrix, colData = coldata, design = .design)
  dds <- DESeq(deseq)
  subset <- rowSums(cpm(counts(dds), normalized = T) >1 ) >= 3
  dds <- dds[subset,]
  rld <- rlog(dds, blind = F)
  
  return(list(deseq = deseq, dds = dds, rld = rld))
  rm(list(deseq, dds, subset, rld))
}

make_pca <- function(rld, group, labs, cols){
  tmp <- prcomp(t(assay(rld)), center =T, scale. = TRUE)
  plot <- ggbiplot(pcobj = tmp, choices = 1:2, scale = 1,
                   groups = group, var.axes = F, circle = T) +
    geom_point(size = 5, aes(color = group)) +
    geom_label_repel(aes(label = rld@colData@rownames, color = group), 
                     size = 5, point.padding = 5, label.padding = .5, 
                     show.legend = F) + 
    scale_color_manual(name = "Conditions",
                       labels = labs,
                       values = cols) +
    theme(axis.text = element_text(size = 12, colour = "darkgrey"),
          axis.title = element_text(size = 12, face = "bold", colour = "black"),
          legend.text = element_text(size = 12, colour = "darkgrey"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))
  
    
  return(plot)
  rm(list(tmp, plot))
}

get_results <- function(dds, sig_log2FC, sig_pval, contrast = NULL, name = NULL){
  tmp <- results(dds, contrast = contrast, name = name,
                 independentFiltering = T, pAdjustMethod = "BH", alpha = 0.05)
  fdr <- fdrtool(tmp$stat, statistic = "normal")
  tmp$padj <- p.adjust(fdr$pval, method = "BH")
  
  df <- data.frame(tmp, geneID = rownames(assay(dds)), row.names = rownames(assay(dds)))
  
  significance <- rep('NS', nrow(df))
  significance[which(abs(df$log2FoldChange) > sig_log2FC & df$pvalue > sig_pval)] <- 'log2FoldChange'
  significance[which(abs(df$log2FoldChange) < sig_log2FC & df$pvalue < sig_pval)] <- '-Log10P'
  significance[which(df$log2FoldChange < -(sig_log2FC) & df$pvalue < sig_pval)] <- 'Signif. down-regulated'
  significance[which(df$log2FoldChange > sig_log2FC & df$pvalue < sig_pval)] <- 'Signif. up-regulated'
  df$significance <- as.factor(significance)
  
  sig_df <- df[which(df$significance == 'Signif. down-regulated' | 
                       df$significance == 'Signif. up-regulated'),]
  
  return(list(results = tmp, df = df, sig_df = sig_df))
  rm(list(tmp, fdr, df, significance, sig_df))
}

MA_plotting <- function(x){
  ggplot(x, aes(baseMean, log2FoldChange, colour = padj)) + geom_point(size=1) + 
    scale_y_continuous(limits = c(-3,3), oob = squish) + scale_x_log10() + 
    geom_hline(yintercept=0, colour = "darkorchid", size = 1, linetype="longdash") + 
    labs(x="Mean of normalized count", y="log2 fold change") + 
    scale_colour_viridis(direction = -1, trans = "sqrt") + theme_bw() + 
    geom_density_2d(colour = "black", size = 1)
}

vulcan_plotting <- function(x){return(ggplot(data = na.omit(x), 
                                          aes(x = log2FoldChange, y = -log10(pvalue), colour = significance)) + 
                                        geom_point(mapping = aes(), inherit.aes = T, size = 2.5) + 
                                        scale_color_manual(
                                          values = c(
                                            'NS'='#c1c1c1',
                                            '-Log10P'='#767676',
                                            'log2FoldChange'='#363636', 
                                            'Signif. down-regulated'='#000f64',
                                            'Signif. up-regulated'='#841f27')) +
                                        scale_x_continuous(expand = expansion(0.1)) + 
                                        labs(
                                          x = expression(paste(log[2], 'FoldChange')),
                                          y = expression(paste(-log[10], italic('P')))) + 
                                        theme(axis.title = element_text(size = 14),
                                              axis.text = element_text(size = 14),
                                              legend.position = 'none') + 
                                        geom_vline(xintercept = c(-1.5, 1.5), linetype = 'dotted', size = 1) + 
                                        geom_hline(yintercept = -log10(0.05), linetype = 'dotted', size = 1) +    
                                        geom_text(
                                          data = subset.data.frame(na.omit(x), pvalue < 0.05 & (log2FoldChange >= 1.5 | log2FoldChange <= -1.5)), 
                                          label = subset.data.frame(na.omit(x), pvalue < 0.05 & (log2FoldChange >= 1.5 | log2FoldChange <= -1.5))$geneID,
                                          hjust = 0, vjust = 1.5, colour = 'black',
                                          position = 'identity', show.legend = F, check_overlap = T))
}


KEGG_plotting <- function(x, y){
  return(enrichKEGG(x,
                    organism = 'hsa',
                    keyType = 'kegg',
                    universe = names(y),
                    pvalueCutoff = 0.05))
}
GO_plotting <- function(x, y, type){
  return(enrichGO(x, 
                  'org.Hs.eg.db', 
                  pvalueCutoff = 0.05,
                  universe = names(y),
                  minGSSize = 5,
                  readable = T,
                  ont = type,
                  pool = T))
}

library(ggplot2)
library(dplyr)
library(forcats)
library(numbers)
library(plyr)
library(stringr)
make_matrix <- function(mat, filter){
  mat <- mat - rowMeans(mat)
  mat <- mat[which(row.names(mat) %in% filter),]
  
  return(mat)
  rm(mat)
}

expression_heatmap <- function(mat, clusters, columns, coldata, annot, method1, method2){
  dend <- dendsort(hclust(dist(mat, method = method1), method = method2))
  dend <- color_branches(dend, k = clusters) # `color_branches()` returns a dendrogram object
  
  col_fun <- colorRamp2(c(min(mat), 0, max(mat)), 
                        c("purple","white","green"))

  # panel_fun = function(index, nm) {
  #   pushViewport(viewport(xscale = range(barplot), yscale = c(0, 2)))
  #   grid.rect()
  #   grid.xaxis(gp = gpar(fontsize = 8))
  #   anno_barplot(x = barplot[,"GO"], baseline = 0, )
  #   popViewport()
  # }
  # anno <- anno_link(align_to = clusters, which("row"), panel_fun = panel_fun,
  #                   size = unit(2, "cm"), gap = unit(1, "cm"), width = unit(4, "cm"))
  hm <- Heatmap(mat, cluster_rows = dend, 
                row_split = clusters, 
                row_title = "Cluster %s",
                row_title_side = "left", row_title_rot = 0,
                show_column_dend = T, show_row_names = F,
                row_dend_reorder = T,
                column_dend_reorder = T,
                column_title = "Conditional Cluster %s",
                name = "normalized expression",
                col = col_fun, border = T, 
                column_split = columns, 
                column_labels = coldata$sample,
                bottom_annotation = annot)
  
  
  return(list(heatmap = hm, dendrogram = dend))
  rm(list(dend, colors, annot, hm))
}

logFC_heatmap <- function(mat, dend){
  row_order <- c(unlist(dend[[1]]), unlist(dend[[2]]))
  col_fun <- colorRamp2(c(min(mat), 0, max(mat)), 
                        c("blue","white","red"))
  
  hm <- Heatmap(mat, name = "log2FC",
                na_col = "white", border = T,
                row_order = row_order,
                column_title = "Expression change",
                row_names_side = "left")
  return(hm)
}

merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

make_simMatrix <- function(list){
  list(BP = calculateSimMatrix(list$BP$ID,
                               orgdb="org.Sc.sgd.db",
                               ont="BP",
                               method="Wang"),
       MF = calculateSimMatrix(list$MF$ID,
                               orgdb="org.Sc.sgd.db",
                               ont="MF",
                               method="Wang"),
       CC = calculateSimMatrix(list$CC$ID,
                               orgdb="org.Sc.sgd.db",
                               ont="CC",
                               method="Wang"))
}

get_reducedTerms <- function(simm, scores, treshold){
  tmp <- reduceSimMatrix(simm,
                         scores,
                         threshold= treshold,
                         orgdb="org.Hs.eg.db")
  pca <- prcomp(simm)
  pca <- as.data.frame(pca$x[,1:2])
  pca <- pca %>% 
    setNames(c('x','y')) %>%
    rownames_to_column(., var = "go")
  
  tmp <- merge(tmp, pca, by = 'go')
  tmp$parentTerm <- as.factor(tmp$parentTerm)
  subset <- tmp [ tmp$parentSimScore > treshold, ]
  
  return(list(reduced = tmp, subset = subset))
}

make_GO_simplot <- function(reduced, subset){
  return(ggplot( data = reduced )
         + geom_point( aes( x = x, y = y, colour = parentTerm, size = score), alpha = I(0.6))
         + geom_point( aes( x = x, y = y, size = score), shape = 21, 
                       fill = "transparent", colour = I (alpha ("black", 0.6) ))
         + stat_ellipse(geom = "polygon", aes( x = x, y = y, colour = parentTerm, fill = parentTerm), type = "norm", level = 0.68, 
                      alpha = .3, linetype = 2, show.legend = F)
         + scale_size( range=c(5, 30))
         + theme_bw()
         + geom_label_repel( data = subset, aes(x = x, y = y, label = term, colour = parentTerm), size = 3 )
         + labs (y = "semantic space x", x = "semantic space y", colour = "Cluster")
         + theme(legend.key = element_blank()))
}

make_FGplot <- function(mydata, type=c("KEGG","GO")){
  attach(mydata)
  
  limits = c(min(Count),max(Count))
  modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
  
  breaks = if (max(Count) > 100) {
    unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
  } else if (max(Count) > 10 & max(Count) < 20){
    unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
  } else if (max(Count) > 20){
    unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
  } else {
    unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
  }
  detach(mydata)
 
  return(
    if(type == "GO"){
      mydata %>%
        mutate(Description = fct_reorder(Description, GeneRatio)) %>%
        group_by(ONTOLOGY) %>%
        dplyr::slice_head(n = 10) %>%
        ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
        geom_point(aes(color = pvalue)) + 
        facet_grid(ONTOLOGY ~ ., scales = "free",
                   labeller = as_labeller(c(
                     'BP' = 'Biological processes',
                     'MF' = 'Molecular functions',
                     'CC' = 'Cellular components'
                   ))) +
        scale_size(range = c(5,15), limits = limits, breaks =  breaks) +
        theme_classic() +
        guides(size = guide_legend(title = "Gene count", order = 2),
               color = guide_colorbar(title = "p-value", order = 1)) + 
        theme(axis.text.x = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              legend.title = element_text(size = 14),
              panel.border = element_rect(colour = "black", fill="transparent",
                                          linewidth = 1),
              panel.grid.major = element_line(colour = "grey80", linewidth = .5),
              legend.text = element_text(size = 14))}
    else{
      mydata %>%
        mutate(Description = fct_reorder(Description, GeneRatio)) %>%
        dplyr::slice_head(n = 10) %>%
        ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
        geom_point(aes(color = pvalue)) + 
        scale_size(range = c(5,15), limits = limits, breaks =  breaks) +
        theme_classic() + 
        guides(size = guide_legend(title = "Gene count", order = 2),
               color = guide_colorbar(title = "p-value", order = 1)) + 
        theme(axis.text.x = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              legend.title = element_text(size = 14),
              legend.text = element_text(size = 14))}
    )
}

get_CategoryExpressionPlot <- function(results, genes){
  subset <- subset.data.frame(results, subset = geneID %in% genes)

  return(ggplot(data = na.omit(results)) 
         + geom_point(mapping = aes(x = log2FoldChange, y = -log10(pvalue), 
                                    colour = significance), size = 2.5) 
         + geom_point(data = subset, mapping = aes(x = log2FoldChange, y = -log10(pvalue)),
                      size = 2.5, fill = "transparent", colour = I (alpha ("yellow", 0.6) )) 
         + scale_color_manual(values = colourPalette) 
         + labs( x = expression(paste(log[2], 'FoldChange')),
                 y = expression(paste(log[10], italic('P')))) 
         + theme(axis.title = element_text(size = 14), 
                 axis.text = element_text(size = 14), 
                 legend.position = 'none') 
         + coord_cartesian(ylim = c(0, -log10(min(results$pvalue))))
         + geom_vline(xintercept = c(-(1.5), 1.5), linetype = 'dotted', size = 1) 
         + geom_hline(yintercept = -log10(0.05), linetype = 'dotted', size = 1) 
         + geom_label_repel(data = subset, aes(x = log2FoldChange, y = -log10(pvalue)),
                            colour = 'black', position = 'identity', max.overlaps = 20,
                            show.legend = F, label.padding = .5, direction = "both",
                            label = paste(subset$geneID)))
  
}

make_gsea <- function(x, GOset){
  set.seed(123) 
  
  return(GSEA(
    geneList = x,
    minGSSize = 25,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    eps = 0,
    seed = TRUE,
    by = "fgsea",
    pAdjustMethod = "BH",
    TERM2GENE = dplyr::select(
      GOset,
      gs_name,
      entrez_gene) %>%
      dplyr::mutate(entrez_gene = as.character(entrez_gene))
    ))
}

make_gseaplot <- function(data){
  return(data %>%
           as.data.frame(.) %>%
           mutate(
             Direction = ifelse(NES > 0, '+','-'),
             Direction = as.factor(Direction),
             Description = fct_reorder(Description, abs(NES))
           ) %>%
           group_by(ONTOLOGY) %>%
           dplyr::slice_head(n = 10) %>%
           ggplot(., aes(x = abs(NES), y = Description, size = abs(NES))) + 
           geom_point(aes(fill = Direction, shape = Direction), color = "black") +
           scale_fill_manual(
             values = c(
               "+" = "brown1",
               "-" = "cornflowerblue"
             ),
             labels = c("+" = "Up-regulated",
                        "-" = "Down-regulated"),
             guide = guide_legend(
               title = "Direction of enrichment",
               order = 1
             ),
             aesthetics = "fill") + 
           scale_shape_manual(
             values = c(
               "+" = 24,
               "-" = 25),
             labels = c("+" = "Up-regulated",
                        "-" = "Down-regulated"),
             guide = guide_legend(
               title = "Direction of enrichment",
               override.aes = list(size = 8),
               order = 1
             )) + 
           scale_size(range = c(6,10),
                      guide = "none") + 
           facet_grid(ONTOLOGY~., scales = "free", #space = "free_x", 
                      labeller = as_labeller(c("BP"="Biological processes",
                                               "MF"="Molecular functions",
                                               "CC"="Cell components"))) +
           theme(axis.text = element_text(size = 14),
                 axis.title = element_text(size = 14),
                 strip.background = element_rect(fill = "white",colour = "grey25",
                                                 linewidth = 1),
                 strip.text = element_text(size = 14, face = "bold"),
                 legend.position = "right",
                 legend.key.size = unit(2,'cm'),
                 legend.key = element_rect(fill = "white"),
                 legend.title = element_text(size = 14),
                 legend.text = element_text(size = 14)) + 
           labs(x = "Absolute normalized enrichment score", y=""))
}

make_GObase <- function (terms, genes) 
{
  colnames(terms) <- tolower(colnames(terms))
  terms$genes <- toupper(terms$genes)
  genes$geneID <- toupper(genes$geneID)
  tgenes <- strsplit(as.vector(terms$genes), ", ")
  if (length(tgenes[[1]]) == 1) 
    tgenes <- strsplit(as.vector(terms$genes), ",")
  count <- sapply(1:length(tgenes), function(x) length(tgenes[[x]]))
  logFC <- sapply(unlist(tgenes), function(x) genes$log2FoldChange[match(x, 
                                                                genes$geneID)])
  if (class(logFC) == "factor") {
    logFC <- gsub(",", ".", gsub("\\.", "", logFC))
    logFC <- as.numeric(logFC)
  }
  s <- 1
  zsc <- c()
  for (c in 1:length(count)) {
    value <- 0
    e <- s + count[c] - 1
    value <- sapply(logFC[s:e], function(x) ifelse(x > 0, 
                                                   1, -1))
    zsc <- c(zsc, sum(value)/sqrt(count[c]))
    s <- e + 1
  }
  terms$zscore <- zsc
  return(terms)
}

compound_GOplot <-function (data) {
  colnames(data) <- tolower(colnames(data))
  subset <- subset.data.frame(data, padj < 0.05)
  subset2 <- subset.data.frame(subset, slim == TRUE)
  
  dummy_col <- data.frame(category = factor(c("BP", "MF", "CC", "KEGG"),
                                            levels = c("BP", "MF", "CC", "KEGG")),
                          padj = data$padj[1:4], zscore = data$zscore[1:4],
                          size = 1:4,
                          count = 1:4)
  data <- data %>% 
    dplyr::mutate(category = factor(category, levels=c("BP", "MF", "CC", "KEGG")))
  
  (plot = ggplot(data) + 
      geom_point(aes(x = zscore, y = -log10(padj), size = count),
                 shape = 21, fill = "grey50", col = "black", alpha = .2) +
      geom_point(data = subset, 
                 aes(x = zscore, y = -log10(padj), fill = category, size = count),
                 shape = 21, col = "black", alpha = .5) +
      geom_rect(data = dummy_col, aes(fill = category),
                xmin = -Inf, xmax = Inf,
                ymin = -Inf, ymax = Inf,
                alpha = 0.1, show.legend = F) + 
      facet_grid(cols = vars(category), #space = "free_x", scales = "free_x",
                 labeller = as_labeller(c("BP"="Biological processes",
                                          "MF"="Molecular functions",
                                          "CC"="Cell components",
                                          "KEGG"="KEGG pathways"))) +
      scale_size(range = c(5, 30), 
                 guide = guide_legend(
                   title = "Gene count",
                   position = "right",
                   override.aes = list(shape = 16, fill = "grey"),
                   order = 2)) + 
      scale_x_continuous(expand = expansion(0.1), n.breaks = 10) + 
      labs(x = "z-score", 
           y = "-log (adj p-value)") + 
      scale_fill_manual(
        name = "Categories",
        labels = c("Biological processes","Molecular functions",
                   "Cell components","KEGG metabolic pathways"),
        values = c("BP" = "chartreuse4", 
                   "MF" = "brown2",
                   "CC" = "cornflowerblue",
                   "KEGG" = "orange"),
        aesthetics = c("fill","color"),
        guide = "none") + 
      geom_hline(yintercept = 1.3, col = "red") + 
      geom_label_repel(data = subset2, 
                       aes(x = zscore, y = -log10(padj), 
                           label = description, fill = category),
                       box.padding = 5, max.overlaps = Inf,
                       arrow = arrow(length = unit(0.015, "npc")),
                       label.padding = .5, show.legend = F) +
      theme_bw(base_size = 14) + 
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 14),
            axis.line = element_line(colour = "grey80"),
            axis.ticks = element_line(colour = "grey80"),
            strip.text.x = element_text(size = 14),
            panel.border = element_rect(fill = "transparent", colour = "grey80"),
            plot.margin = margin(.5,1,.5,1, 'cm'),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 12),
            legend.position = "bottom",
            panel.grid.major = element_line(colour = "grey80"),
            plot.background = element_rect(color = "white")))
}