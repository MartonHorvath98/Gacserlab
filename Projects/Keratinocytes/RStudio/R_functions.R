################################################################################
# 1. Utility functions                                                         #
################################################################################

# Pause the script until the user provides a valid answer
user_answer <- function(section){
  answer <- ""
  # Loop until the user provides a valid answer
  while(!(answer %in% c("yes", "no"))) {
    # Prompt the user
    answer <- readline(prompt="Do you wish to proceed? (yes/no): ")
    
    # Check the user's answer
    if(answer %in% c("yes", "no")) {
      break # Exit the loop if answer is valid
    } else {
      cat("Please answer with 'yes' or 'no'.\n") # Prompt again if answer is invalid
    }
  }
  
  # Act based on the user's answer
  if(answer == "yes") {
    message(str_glue("Continuing with {section}...\n"))
    return(TRUE)
  } else {
    message(str_glue("You chose to skip {section} ...\n"))
    return(FALSE)
  }
}

# Recursively merges a list of data frames into a single data frame
merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

################################################################################
# 1. Functions for differential expression analysis                            #
################################################################################

# Function to perform differential expression analysis using DESeq2
make_deseq <- function(matrix, coldata, design){
  # Create a formula for the design matrix
  .design <- reformulate(design)
  # Create a DESeqDataSet object
  deseq <- DESeqDataSetFromMatrix(matrix, colData = coldata, design = .design)
  # Perform differential expression analysis
  dds <- DESeq(deseq)
  # Remove features with low counts, <1 read per million in at least 3 samples
  subset <- rowSums(cpm(counts(dds), normalized = T) >1 ) >= 3
  dds <- dds[subset,]
  # Perform variance stabilizing transformation
  rld <- rlog(dds, blind = F)
  
  return(list(deseq = deseq, dds = dds, rld = rld))
  rm(list(deseq, dds, subset, rld))
}

# Function to extract a result table from the DESeq analysis
get_results <- function(dds, sig_log2FC, sig_pval, contrast = NULL, name = NULL){
  # Extract results table
  tmp <- results(dds, contrast = contrast, name = name,
                 independentFiltering = T, pAdjustMethod = "BH", alpha = 0.05)
  # Perform multiple testing correction
  fdr <- fdrtool(tmp$stat, statistic = "normal")
  tmp$padj <- p.adjust(fdr$pval, method = "BH")
  # Create a data frame from the results table
  df <- data.frame(tmp, geneID = rownames(assay(dds)), 
                   row.names = rownames(assay(dds)))
  
  # Add a column for significance using the specified thresholds
  significance <- rep('NS', nrow(df))
  significance[which(abs(df$log2FoldChange) > sig_log2FC & df$pvalue > sig_pval)] <- 'log2FoldChange'
  significance[which(abs(df$log2FoldChange) < sig_log2FC & df$pvalue < sig_pval)] <- '-Log10P'
  significance[which(df$log2FoldChange < -(sig_log2FC) & df$pvalue < sig_pval)] <- 'Signif. down-regulated'
  significance[which(df$log2FoldChange > sig_log2FC & df$pvalue < sig_pval)] <- 'Signif. up-regulated'
  df$significance <- as.factor(significance)
  # Create a data frame with only significant results
  sig_df <- df[which(df$significance == 'Signif. down-regulated' | 
                       df$significance == 'Signif. up-regulated'),]
  
  return(list(results = tmp, df = df, sig_df = sig_df))
  rm(list(tmp, fdr, df, significance, sig_df))
}

# Function to annotate the results table with gene names and entrez IDs
annotate_results <- function(x){
  return(x %>%
           dplyr::mutate(
             # map the ensembl ID to the gene name and entrezID
             ensembl = geneID,
             geneName = mapIds(org.Hs.eg.db, x[["geneID"]],
                               keytype = "ENSEMBL", column = "SYMBOL",
                               multiVals = "first"),
             entrezID = mapIds(org.Hs.eg.db, x[["geneID"]],
                               keytype = "ENSEMBL", column = "ENTREZID",
                               multiVals = "first")) %>%
           dplyr::select(!geneID) %>%
           # move the ID columns to the front
           dplyr::relocate(where(is.character), .before = where(is.numeric)) %>%
           dplyr::relocate(significance, .after = everything()))
}

#
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
  #rm(list(tmp, plot))
}

make_matrix <- function(mat, filter){
  mat <- mat - rowMeans(mat)
  mat <- mat[which(row.names(mat) %in% filter),]
  
  return(mat)
  rm(mat)
}

choose_method <- function(
    matrix, n_clust, distance_options = c(), clustering_options = c()
){
  
  best_avg_sil_width <- -Inf
  best_method1 <- ""
  best_method2 <- ""
  
  for (distance in distance_options) {
    for (cluster in clustering_options) {
      # Compute the distance matrix
      dist_mat <- dist(matrix, method = distance)
      
      # Perform hierarchical clustering
      hc <- hclust(dist_mat, method = cluster)
      
      # Cut the tree to obtain a reasonable number of clusters, e.g., k
      # The choice of k depends on your specific dataset and goals
      cluster_assignment <- cutree(hc, k = n_clust)
      
      # Calculate silhouette widths for the clustering
      silhouette_info <- silhouette(cluster_assignment, dist_mat)
      avg_sil_width <- mean(silhouette_info[, "sil_width"])
      
      # Check if this is the best combination so far
      if (avg_sil_width > best_avg_sil_width) {
        best_avg_sil_width <- avg_sil_width
        best_dist <- distance
        best_clust <- cluster
      }
    }
  }
  
  # Print the best combination
  cat("Best combination:\n")
  cat("Distance metric:", best_dist, "\n")
  cat("Clustering method:", best_clust, "\n")
  cat("With an average silhouette width of:", best_avg_sil_width, "\n")
  return(list(best_dist, best_clust))
}

expression_heatmap <- function(mat, clusters, columns, coldata, annot, method1, method2){
  dend <- dendsort(hclust(dist(mat, method = method1), method = method2))
  dend <- color_branches(dend, k = clusters) 
  
  col_fun <- colorRamp2(c(min(mat), 0, max(mat)), 
                        c("purple","white","green"))
  
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

arrange_regions <- function(set1, set2, set3, set4, names, filename){
  # Select meaningful columns from the input data
  tmp <- list(
    set1[,c("ensembl","log2FoldChange")],
    set2[,c("ensembl","log2FoldChange")],
    set3[,c("ensembl","log2FoldChange")],
    set4[,c("ensembl","log2FoldChange")]
  )
  names(tmp) <- names
  
  # Merge the data frames
  table <- merge.rec(tmp, by = "ensembl", all = T) %>%
    setNames(.,c("geneID", names))
  
  # Extract the allocation of each gene
  regions <- table %>%
    # Pivot the data frame to long format along the expression values
    tidyr::pivot_longer(cols = -geneID, names_to = "Region", values_to = "Value") %>%
    # Filter out missing values
    dplyr::filter(!is.na(Value)) %>%
    # Group by geneID and concatenate the regions for genes in intersecting sets
    dplyr::group_by(geneID) %>%
    dplyr::summarise(Regions = paste(Region, collapse = "")) %>%
    dplyr::ungroup() 
  
  # Merge the regions with the table
  table <- merge(table, regions, by = "geneID")
  
  # Calculate the number of genes in each subset
  regions <- regions %>% 
    dplyr::group_by(Regions) %>% 
    dplyr::summarise(Size = n()) %>% 
    dplyr::pull(Size, name = Regions)
  
  return(list(table = table, regions = regions)) 
}

make_venn <- function(expression, size, names = list()){
  # Create a data frame for plotting the ellipses on the Venn diagram
  region.matrix <- data.frame(
    labels = names(names),
    # Define the coordinates and dimensions of the ellipses: A, B, C, D
    x0=c(0.35, 0.5, 0.5, 0.65),
    y0=c(0.47, 0.57, 0.57, 0.47),
    a=c(0.35, 0.35, 0.35, 0.35),
    b=c(0.2, 0.15, 0.15, 0.2),
    angle=c(-10, -10, 10, 10),
    stringsAsFactors = T)
    
  # Create a data frame for plotting the labels on the Venn diagram
  label.matrix <- matrix(byrow = T, nrow = 15, ncol = 3,
                         # The labels equal the number of genes in each region
                         c(size[['B']], 0.35, 0.75,
                           size[['BC']], 0.5, 0.69,
                           size[['C']], 0.65, 0.75,
                           size[['AB']], 0.27, 0.65,
                           size[['ABC']], 0.4, 0.6,
                           size[['ABCD']], 0.5, 0.5,
                           size[['BCD']], 0.6, 0.6,
                           size[['CD']], 0.73, 0.65,
                           size[['A']], 0.15, 0.6,
                           size[['AC']], 0.3, 0.45,
                           size[['ACD']], 0.4, 0.4,
                           size[['ABD']], 0.6, 0.4,
                           size[['BD']], 0.7, 0.45,
                           size[['D']], 0.85, 0.6,
                           size[['AD']], 0.5, 0.28))
  colnames(label.matrix) <- c("label", "x", "y")
  row.names(label.matrix) <- c('B','BC','C','AB','ABC','ABCD','BCD','CD','A',
                               'AC','ACD','ABD','BD','D','AD')
  
  set.seed(0)
  # Create a data frame for plotting the genes on the Venn diagram
  point.matrix <- merge(human.venn$table, label.matrix %>%
                          as.data.frame() %>% 
                          tibble::rownames_to_column("Regions"), 
                        by = "Regions")
  point.matrix <- point.matrix %>% 
    rowwise(.) %>% 
    # Add jitter to the coordinates for better visualization
    dplyr::mutate(
      x = rnorm(1, x, sd = 0.015),
      y = rnorm(1, y, sd = 0.015),
      # Determine the trend of gene expression
      trend = dplyr::case_when(
        all(dplyr::c_across(A:D) > 0, na.rm = T) ~ 'UP',
        all(dplyr::c_across(A:D) < 0, na.rm = T) ~ 'DOWN',
        TRUE ~ 'CHANGE'
      ),
      trend = as.factor(trend)
    ) %>% 
    ungroup(.)
  
  # Create the Venn diagram
  plot <- ggplot() + 
           # Plot the ellipses
           ggforce::geom_ellipse(
             data = region.matrix, # ellipse coordinates
             aes(x0=x0,y0=y0,a=a,b=b, angle=angle, fill = labels),
             color = "grey25", linewidth = .8, alpha = .3) + 
           # Paint the ellipses
           scale_fill_manual(
             values = c("A" = "steelblue",
                        "B" = "salmon",
                        "C" = "darkred",
                        "D" = "darkblue"),
             labels = c(paste("A", names[["A"]], sep = ") "),
                        paste("B", names[["B"]], sep = ") "),
                        paste("C", names[["C"]], sep = ") "),
                        paste("D", names[["D"]], sep = ") ")),
             guide = guide_legend(title = "Regions",
                                  order = 1,
                                  override.aes = list(
                                    linetype = c(0,0,0,0),
                                    alpha = 0))) +
           geom_point(
             data = point.matrix, # Gene coordinates
             aes(x = x, y = y, colour = trend),
             size = 3) +
           # Paint the genes
           scale_color_manual(
             name = "Trend",
             values = c("UP" = "red",
                        "DOWN" = "blue",
                        "CHANGE" = "green"),
             labels = c("Up-regulated", 
                        "Down-regulated",
                        "Opposing changes")) +
           # Plot the number of genes in each region
           geom_label(
             data = label.matrix, # Label coordinates
             aes(x = x, y = y, label = label),
             size = 6, alpha = 0.5, fill = 'white', color = 'grey15') + 
           # Clear the background and axes
           theme_dendro() +
           # Increase font size and adjust legend placement
           theme(text = element_text(size = 16),
                 legend.key.spacing = unit(.5, 'cm'),
                 legend.spacing = unit(2, 'cm')) + 
           # Expand the plot limits
           expand_limits(x = 0, y = .9) + 
           # Add region labels
           annotate('label', x = 0, y = 0.6, size = 8, label = 'A') +
           annotate('label', x = 0.25, y = 0.85, size = 8, label = 'B') +
           annotate('label', x = 0.75, y = 0.85, size = 8, label = 'C') + 
           annotate('label', x = 1, y = 0.6, size = 8, label = 'D') +
           # Add a title
           labs(title = "Venn diagram of gene expression regions")
  return(list(plot = plot, point_matrix = point.matrix))
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
                                          label = subset.data.frame(na.omit(x), pvalue < 0.05 & (log2FoldChange >= 1.5 | log2FoldChange <= -1.5))$geneName,
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
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) + 
        theme_classic() +
        guides(size = guide_legend(title = "Gene count", order = 2),
               color = guide_colorbar(title = "p-value", order = 1)) + 
        theme(axis.text.x = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              legend.text = element_text(size = 14),
              panel.border = element_rect(colour = "black", fill="transparent",
                                          linewidth = 1),
              panel.grid.major = element_line(colour = "grey80", linewidth = .5),
              strip.text = element_text(size = 14))+
        labs(y = "", x = "Gene ratio")}
    else{
      mydata %>%
        mutate(Description = fct_reorder(Description, GeneRatio)) %>%
        dplyr::slice_head(n = 10) %>%
        ggplot(aes(x = GeneRatio, y = Description,
                   size = as.numeric(Count))) + 
        geom_point(aes(color = pvalue)) + 
        scale_size(range = c(5,15), limits = limits, breaks =  breaks) +
        scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) + 
        theme_classic() + 
        guides(size = guide_legend(title = "Gene count", order = 2),
               color = guide_colorbar(title = "p-value", order = 1)) + 
        theme(axis.text.x = element_text(size = 14),
              axis.title.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              legend.text = element_text(size = 14)) +
        labs(y = "", x = "Gene ratio")}
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
      geom_label_repel(data = subset2, force = 0.5,
                       aes(x = zscore, y = -log10(padj), 
                           label = stringr::str_wrap(description, 35),
                           fill = category), alpha = .8,
                       fontface = 'bold',
                       box.padding = 5, max.overlaps = Inf,
                       arrow = arrow(length = unit(0.015, "npc")),
                       label.padding = .5, show.legend = F) +
      theme_bw(base_size = 14) + 
      theme(axis.title = element_text(size = 16),
            axis.text = element_text(size = 16),
            axis.line = element_line(colour = "grey80"),
            axis.ticks = element_line(colour = "grey80"),
            strip.text.x = element_text(size = 16),
            panel.border = element_rect(fill = "transparent", colour = "grey80"),
            plot.margin = margin(.5,1,.5,1, 'cm'),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            legend.text = element_text(size = 14),
            legend.position = "bottom",
            panel.grid.major = element_line(colour = "grey80"),
            plot.background = element_rect(color = "white")))
}

# Function to build a phylogenetic tree
buildTree <- function(df) {
  # Create a root node
  root <- Node$new("47608")
  
  # Recursive function to add children
  addChildren <- function(node) {
    children <- df[df$from == node$name, "to"]
    for (child in children) {
      # Add child node with additional information
      childNode <- node$AddChild(child)
      addChildren(childNode) # Recursion for any children of the current child
    }
  }
  # Start the tree building process
  addChildren(root)
  return(root)
}

make_spiderbase <- function(df){
  
  colnames(df) <- tolower(colnames(df))
  bg.genes = df$geneid[which(df$panc == T)]
  bg.count = length(bg.genes)
  
  genes = list(
    ca = df$geneid[which(df$ca11 == T & df$geneid %in% bg.genes)],
    cp = df$geneid[which(df$cp51 == T & df$geneid %in% bg.genes)]  
  )
  
  counts = lapply(genes, length)
  
  zsc = list()
  sub = list()
  for(i in names(genes)){
    sub[[i]] = subset.data.frame(df, df$geneid %in% genes[[i]])
    
    value = 0
    value = sum(ifelse(sub[[i]]$trend == "Signif. up-regulated", 1, -1))
    
    zsc[[i]] = value/sqrt(counts[[i]])
  }
  
  ratios = lapply(counts, function(x) round((x/bg.count)*100, digits = 1))
  
  return(
    data.frame(
      condition = c("CA11","CP51"),
      genes = sapply(genes, paste, collapse = ", "),
      count = unlist(counts),
      geneRation = unlist(ratios),
      zscore = unlist(zsc)
    ))
}

make_circPlot <- function(df){
  df <- df %>% 
    dplyr::mutate(GeneRatio = round(GeneRatio, digits = 3),
                  term = as.factor(term)) %>% 
    dplyr::group_by(condition, term) %>% 
    dplyr::slice_head(n = 1) 
  
  return(ggplot(df,
          aes(
            x = reorder(str_wrap(term, 20), Count),
            y = Count,
            fill = zscore,
            color = zscore
          )) +
     geom_bar(
       position = "dodge2",
       stat='identity',
       show.legend = TRUE,
       alpha = .9
     ) + 
     geom_label_repel(
       data = subset.data.frame(df, Count > 0),
       aes(label = round(zscore, 2), y = max(Count)/2), direction = "y",
       fill="white",label.padding = unit(5, "pt"),
       position=position_dodge(width=0.9)) +
     scale_fill_gradient(
       "Activation (z-score)",
       low = "grey50", high = "darkred", 
       limits = c(0, 5), breaks = c(0, 1, 2, 3, 4, 5),
       aesthetics = c("fill","color")
     ) + 
     facet_grid(.~condition, drop = T, scales = "free",
                labeller = as_labeller(
                  c("CTRL_vs_C.albi_HaCat" = "C.albicans (MOI 1:1) vs HaCat",
                    "CTRL_vs_C.para_HaCat" = "C.parapsilosis (MOI 5:1) vs HaCat",
                    "CTRL_vs_C.albi_HPV-KER" = "C.albicans (MOI 1:1) vs HPV-KER",
                    "CTRL_vs_C.para_HPV-KER" = "C.parapsilosis (MOI 5:1) vs HPV-KER")
                )) +
     labs(x = "Gene count") +
     theme_minimal() + 
     theme(
       # Remove axis ticks and text
       axis.ticks = element_blank(),
       # Use grey text for the region names
       line = element_line(color = "grey80", linewidth = 1),
       axis.title.x =  element_text(color = "grey12", size = 14),
       axis.title.y = element_blank(),
       axis.text.y = element_text(color = "grey12", size = 14),
       axis.text.x = element_text(color = "grey12", size = 12),
       panel.background = element_rect(fill = "white", color = "white"),
       #strip parameters
       strip.clip = "off", 
       strip.text = element_text(color = "grey12", size = 16),
       # Move the legend to the bottom
       legend.position = "bottom",
       legend.title = element_text(size = 14),
     ) +
     scale_y_continuous(expand = expansion(0.1), n.breaks = 5) +
     coord_radial(clip = "off", expand = T, r_axis_inside = T, direction = -1,
                  start = -0.95 * pi, end = 0.95 * pi, rotate_angle = T, inner.radius = 0.1))
}

