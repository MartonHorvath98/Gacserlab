###############################
# Load files for the analysis #
###############################
require_file <- function(file,...){
  if (file.exists(file)) {
    return(read.csv(file = file,...))
  }}

####################################
# Differential expression analysis #
####################################
calc_DiffExp <- function(matrix, coldata, design){
  .design <- reformulate(design)
  deseq <- DESeqDataSetFromMatrix(matrix, colData = coldata, design = .design)
  dds <- DESeq(deseq)
  
  dds_norm <- vst(dds)
  
  return(list(deseq = deseq, dds = dds, dds_norm = dds_norm))
}

get_results <- function(dds, contrast = NULL, name = NULL){
  res <- results(dds, contrast = contrast, name = name,
                 independentFiltering = T, pAdjustMethod = "BH", alpha = 0.05)
  fdr <- fdrtool(res$stat, statistic = "normal")
  res$padj <- p.adjust(fdr$pval, method = "BH")
  
  tmp <- stringr::str_split(contrast,'_')
  resLFC <- lfcShrink(dds, contrast = contrast, 
                      res = res, type = "ashr")
  
  return(list(res = res, resLFC = resLFC))
}

sig_results <- function(results, sig_log2FC, sig_pval){
  df <- data.frame(results, geneID = rownames(results), row.names = rownames(results))
  df <- df %>%
    dplyr::mutate(geneID = mapIds(org.Hs.eg.db, 
                                  row.names(.), 
                                  keytype = "ENSEMBL", 
                                  column = "SYMBOL",
                                  multiVals = "first"),
                  entrezID = mapIds(org.Hs.eg.db, 
                                    row.names(.), 
                                    keytype = "ENSEMBL", 
                                    column = "ENTREZID",
                                    multiVals = "first")) %>%
    dplyr::mutate(significance = dplyr::case_when(abs(log2FoldChange) > sig_log2FC & 
                                             pvalue > sig_pval ~ 'log2FoldChange',
                                           abs(log2FoldChange) < sig_log2FC & 
                                             pvalue < sig_pval ~ '-Log10P',
                                           log2FoldChange < (-1)*sig_log2FC & 
                                             pvalue < sig_pval ~ 'Signif. down-regulated',
                                           log2FoldChange > sig_log2FC & 
                                             pvalue < sig_pval ~ 'Signif. up-regulated',
                                           T ~ 'NS')) %>%
    dplyr::relocate(c("geneID","entrezID"), .after = everything())
 
  sig_df <- df %>%
    dplyr::filter(significance == 'Signif. down-regulated' | 
                  significance == 'Signif. up-regulated') %>%
    dplyr::filter(!is.na(geneID))

  return(list(df = df, sig_df = sig_df))
}

merge.rec <- function(.list, ...){
  if(length(.list)==1) return(.list[[1]])
  Recall(c(list(merge(.list[[1]], .list[[2]], ...)), .list[-(1:2)]), ...)
}

## calculate median expression level
diffexp_GEO <- function(expr, design){
  cutoff <- median(expr)
  is_expressed <- expr > cutoff
  keep <- rowSums(is_expressed) > 2
  table(keep)
  
  expr <- expr[keep,]
  aw <- arrayWeights(expr, design)
  
  fit <- lmFit(expr, design, weights = aw)
  
  levels <- colnames(design)
  contrasts <- makeContrasts(contrasts = paste(levels[2],"-", levels[1]), levels=design)
  
  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  
  table <- topTable(fit2, number=Inf)
  
  return(list(fit = summary(fit2), df = table))
}


get_CategoryExpressionPlot <- function(results,sig_log2FC, sig_pval, names){
  subset <- subset.data.frame(results, subset = geneID %in% names) 
  
  return(ggplot(data = na.omit(results)) 
         + geom_point(mapping = aes(x = log2FoldChange, y = -log10(pvalue), 
                                    colour = significance), size = 2.5) 
         + geom_point(data = subset, mapping = aes(x = log2FoldChange, y = -log10(pvalue)),
                      size = 2.5, fill = "transparent", colour = I (alpha ("yellow", 0.6) )) 
         + scale_color_manual(values = c("Signif. down-regulated" = "#000f64", "Signif. up-regulated" = "#841f27")) 
         + labs( x = expression(paste(log[2], 'FoldChange')),
                 y = expression(paste(log[10], italic('P')))) 
         + theme(axis.title = element_text(size = 14), 
                 axis.text = element_text(size = 14), 
                 legend.position = 'none') 
         + coord_cartesian(ylim = c(0, -log10(min(results$pvalue))))
         + geom_vline(xintercept = c(-(sig_log2FC), sig_log2FC), linetype = 'dotted', size = 1) 
         + geom_hline(yintercept = -log10(sig_pval), linetype = 'dotted', size = 1) 
         + geom_label_repel(data = subset, aes(x = log2FoldChange, y = -log10(pvalue)),
                            colour = 'black', position = 'identity', 
                            show.legend = F, label.padding = .5, direction = "both",
                            label = paste(subset$geneID)))
  
}
###############################
# Overrepresentation analysis #
###############################
get_entrez <- function(list){
  background <- list$df %>%
    dplyr::pull("entrezID")
  background <- unique(as.character(background))
  
  interest <- list$sig_df %>%
    dplyr::pull("entrezID")
  interest <- unique(as.character(interest))
  
  return(list(background = background, interest = interest))
}

kegg_results <- function(genes_of_interest, background_set, kegg_df){
  kegg <- enricher( gene = genes_of_interest, 
                     pvalueCutoff = 0.1, 
                     pAdjustMethod = "BH", 
                     universe = background_set,
                     TERM2GENE = dplyr::select(
                       kegg_df,
                       gs_name,
                       human_entrez_gene)
                     )
  kegg <- setReadable(kegg, org.Hs.eg.db, keyType = "ENTREZID")
  
  df <- as.data.frame(kegg) 
  df <- df[order(df$pvalue, decreasing = F),]
  df <- df %>% 
    mutate(GeneRatio = sapply(stringr::str_split(df$GeneRatio, "/"), 
                              function(x) as.numeric(x[1])/as.numeric(x[2])))
  return(list(kegg = kegg, df = df))
}

gsea_results <- function(df, hallmark_set){
  dup_genes <- df %>%
    dplyr::filter(!is.na(geneID)) %>%
    dplyr::filter(duplicated(geneID)) %>%
    dplyr::pull(geneID)
  
  dup_df <- df %>%
    dplyr::filter(geneID %in% dup_genes) %>%
    dplyr::arrange(geneID) %>%
    dplyr::group_by(geneID) %>%
    arrange(dplyr::desc(abs(log2FoldChange)))
  
  df <- df %>%
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::distinct(geneID, .keep_all = TRUE)
  
  genelist <- df$log2FoldChange
  names(genelist) <- df$geneID
  genelist <- sort(genelist, decreasing = TRUE)
  
  set.seed(42)
  gsea <- GSEA(
    geneList = genelist, 
    minGSSize = 25, 
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    eps = 0, 
    seed = TRUE, 
    pAdjustMethod = "BH",
    TERM2GENE = dplyr::select(
      hallmark_set,
      gs_name,
      gene_symbol
    )
  )
  gsea_df <- as.data.frame(gsea@result)
  gsea_df <- gsea_df %>%
    dplyr::slice_max(abs(NES), n = 10) %>%
    dplyr::mutate(Direction = dplyr::case_when(NES > 0 ~ '+',
                                        NES < 0 ~ '-'))
  return(list(dup_df = dup_df, gsea = gsea, gsea_df = gsea_df))
}


####################################
# Visualisation of gene expression #
####################################
make_pca <- function(vst, group){
  tmp <- prcomp(t(assay(vst)), center =T, scale. = TRUE)
  plot <- ggbiplot(pcobj = tmp, choices = 1:2, scale = 1,
                   groups = group, var.axes = F, circle = T) +
    geom_point(size = 3, aes(color = group)) +
    geom_label_repel(aes(label = vst@colData@key, color = group), size = 5)
  
  
  return(plot)
}

make_dotplot <- function(mydata, Count){
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
  
  mydata %>%
    dplyr::arrange(desc(GeneRatio)) %>%
    dplyr::mutate(Description = fct_reorder(Description, GeneRatio)) %>%
    dplyr::slice(1:10) %>%
    ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
    geom_point(aes(color = pvalue)) + 
    scale_color_gradient(low = "blue",
                         high = "deeppink3",
                         space = "Lab",
                         na.value = "grey50",
                         guide = "colourbar",
                         aesthetics = "colour") +
    scale_size(range = c(4,10), limits = limits, breaks =  breaks) + 
    theme_classic() + 
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)) + 
    labs(main = names(mydata), size = "Count", y = "")
}

make_vulcanplot <- function(total, sig){
  colourPalette <- c('#c1c1c1', '#363636', '#000f64','#841f27')
  
  return((ggplot(data = na.omit(total), 
                 aes(x = log2FoldChange, 
                     y = -log10(pvalue), 
                     colour = significance)) 
          + geom_point(mapping = aes(), inherit.aes = T, size = 2.5) 
          + scale_color_manual(values = colourPalette)
          + labs(x = expression(paste(log[2], 'FoldChange')),
                 y = expression(paste(log[10], italic('P')))) 
          + theme(axis.title = element_text(size = 14), 
                  axis.text = element_text(size = 14), 
                  legend.position = 'none') 
          + geom_vline(xintercept = c(-1.5, 1.5), 
                       linetype = 'dotted', size = 1) 
          + geom_hline(yintercept = -log10(0.05), 
                       linetype = 'dotted', size = 1) 
          + geom_text(data = sig, hjust = 0, vjust = 1.5, 
                      colour = 'black', position = 'identity', 
                      show.legend = F, check_overlap = T,
                      label = sig$geneID)))
}

make_gseaplot <- function(data){
  return((ggplot(data, aes(x = abs(enrichmentScore), y = ID, size = abs(NES))) 
          + geom_point(aes(fill = pvalue, shape = Direction), color = "white") 
          + scale_color_gradient(low = "blue",
                         high = "deeppink3",
                         space = "Lab",
                         na.value = "grey50",
                         guide = "colourbar",
                         aesthetics = "fill") 
          + scale_shape_manual(values = c(
            "+" = 24,
            "-" = 25
          ))
          + scale_size(range = c(6,10))
          + theme_classic() 
          + theme(axis.text.x = element_text(size = 14),
                  axis.title.x = element_text(size = 14),
                  legend.title = element_text(size = 14),
                  legend.text = element_text(size = 14))
          + labs(main = names(data), 
                 size = "Normalized enrichment score",
                 shape = "Direction of enrichment",
                 x = "Enrichment score", y="")))
}

make_simMatrix <- function(df, type){
  df <- as.data.frame(df)
  simM <- calculateSimMatrix(df$ID,
                             orgdb="org.Hs.eg.db",
                             ont=type,
                             method="Wang")
  
  scores <- setNames(-log10(df$pvalue), df$ID)
  reduced <- reduceSimMatrix(simM,
                             scores,
                             threshold=0.9,
                             orgdb="org.Hs.eg.db")
  
  pca <- prcomp(simM)
  pca <- as.data.frame(pca$x[,1:2])
  pca <- pca %>% 
    setNames(c('x','y')) %>%
    rownames_to_column(., var = "go")
  
  reduced <- merge(reduced, pca, by = 'go')
  reduced$parentTerm <- as.factor(reduced$parentTerm)
  
  subset <- subset.data.frame(reduced, reduced$term == reduced$parentTerm)
  
  return(list(simM = simM, reduced = reduced, subset = subset))
}

make_GO_simplot <- function(reduced, subset){
  return(ggplot( data = reduced) +
           stat_ellipse(geom = "polygon", aes( x = x, y = y, colour = parentTerm, fill = parentTerm), type = "norm", level = 0.68, 
                        alpha = .15, linetype = 2, show.legend = F) +
           geom_point( aes( x = x, y = y, colour = parentTerm, size = score), alpha = I(0.6), show.legend = F) + 
           geom_point( aes( x = x, y = y, size = score), shape = 21, 
                       fill = "transparent", colour = I (alpha ("black", 0.6) ), show.legend = F) + 
           geom_label_repel(data = subset, aes(x = x, y = y, label = term, colour = parentTerm), show.legend = F) + 
           labs (y = "semantic space x", x = "semantic space y") + 
           theme(legend.key = element_blank()))
}


##############################
# Venn diagram & upset plots #
##############################
make_venn <- function(set1, set2, set3, names, filename){
  tmp <- list(
    set1[,c("geneID","log2FoldChange")],
    set2[,c("geneID","log2FoldChange")],
    set3[,c("geneID","log2FoldChange")]
  )
  names(tmp) <- names
  venn <- venn.diagram(x = list(tmp[[1]]$geneID, tmp[[2]]$geneID, tmp[[3]]$geneID),
                       category.names = names(tmp),
                       filename = paste0(plots_dir, "/", filename,".png"), 
                       output = T,
                       
                       # Output features
                       imagetype="png" ,
                       height = 960 , 
                       width = 960 , 
                       resolution = 300,
                       compression = "lzw",
                       # Circles
                       lwd = 2,
                       lty = 'blank',
                       fill = brewer.pal(3, "Pastel2"),
                       
                       # Numbers
                       cex = .6,
                       fontface = "bold",
                       fontfamily = "sans",
                       
                       # Set names
                       cat.cex = 0.6,
                       cat.fontface = "bold",
                       cat.default.pos = "outer",
                       cat.pos = c(-27, 27, 135),
                       cat.dist = c(0.055, 0.055, 0.085),
                       cat.fontfamily = "sans",
                       rotation = 1)
  
  list <- GOVenn(tmp[[1]], tmp[[2]], tmp[[3]], plot = F)$table
  list <- lapply(list, rownames_to_column, var = "geneID")
  list[["A_only"]] <- list[["A_only"]] %>% dplyr::rename(logFC_A = logFC)
  list[["B_only"]] <- list[["B_only"]] %>% dplyr::rename(logFC_B = logFC)
  list[["C_only"]] <- list[["C_only"]] %>% dplyr::rename(logFC_C = logFC)
   
  table <- do.call(rbind.fill,rev(list))
  table <- table %>% dplyr::distinct(geneID, .keep_all = T)
   
  return(list(venn = venn, table = table))
}

make_upsetbase <- function(table, vars){
  table <- table %>%
    dplyr::mutate(var1 = ifelse(is.na(logFC_A), FALSE, TRUE)) %>%
    dplyr::mutate(var2 = ifelse(is.na(logFC_B), FALSE, TRUE)) %>%
    dplyr::mutate(var3 = ifelse(is.na(logFC_C), FALSE, TRUE)) %>%
    relocate(where(is.logical), .before = where(is.character))
  
  colnames(table)[1:3] <- vars
  
  return(table)
}



