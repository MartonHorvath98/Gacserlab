condition_palette <- c("ctrl" = "#869B5B",
                       "PANC"="green",
                       "CA11"="#FF7751",
                       "CP11"="#81EDF7",
                       "CP51"="#00A4C0")
regulation_palette <- c("Signif. up-regulated" = "#841f27",
                        "Signif. down-regulated" = "#000f64")

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
  
  subset <- rowSums(cpm(counts(dds), normalized = T) >1 ) >= 3
  dds <- dds[subset,]
  
  dds_norm <- rlog(dds)
  
  return(list(deseq = deseq, dds = dds, dds_norm = dds_norm))
}

get_results <- function(dds, contrast = NULL, name = NULL, lfc_treshold, pval_treshold){
  # Extract results table
  tmp <- results(dds, contrast = contrast, name = name,
                 independentFiltering = T, pAdjustMethod = "BH", alpha = 0.05)
  # Perform multiple testing correction
  fdr <- fdrtool(tmp$stat, statistic = "normal", plot = F)
  tmp$padj <- p.adjust(fdr$pval, method = "BH")
  # Create a data frame from the results table
  df <- data.frame(tmp, geneID = rownames(tmp), row.names = rownames(tmp))
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
    dplyr::mutate(significance = dplyr::case_when(abs(log2FoldChange) > lfc_treshold & 
                                                    padj > pval_treshold ~ 'log2FoldChange',
                                                  abs(log2FoldChange) < lfc_treshold & 
                                                    padj < pval_treshold ~ 'Log10P',
                                                  log2FoldChange < (-1)*lfc_treshold & 
                                                    padj < pval_treshold ~ 'Signif. down-regulated',
                                                  log2FoldChange > lfc_treshold & 
                                                    padj < pval_treshold ~ 'Signif. up-regulated',
                                                  T ~ 'NS')) %>%
    dplyr::relocate(c("geneID","entrezID"), .after = everything())
  
  sig_df <- df %>%
    dplyr::filter(significance == 'Signif. down-regulated' | 
                    significance == 'Signif. up-regulated') %>%
    dplyr::filter(!is.na(geneID))
  
  return(list(df = df, sig_df = sig_df))
}

miR_results <- function(dds, contrast = NULL, name = NULL, lfc_treshold, p_treshold, tissue){
  # Extract results table
  tmp <- results(dds, contrast = contrast, name = name,
                 independentFiltering = T, pAdjustMethod = "BH", alpha = 0.05)
  # Perform multiple testing correction
  fdr <- fdrtool(tmp$stat, statistic = "normal", plot = F)
  tmp$padj <- p.adjust(fdr$pval, method = "BH")
  # Create a data frame from the results table
  df <- data.frame(tmp, miRname = rownames(tmp), row.names = rownames(tmp))
  df <- df %>%
    dplyr::mutate(significance = dplyr::case_when(abs(log2FoldChange) > lfc_treshold & 
                                                    padj > p_treshold ~ 'log2FoldChange',
                                                  abs(log2FoldChange) < lfc_treshold & 
                                                    padj < p_treshold ~ 'Log10P',
                                                  log2FoldChange < (-1)*lfc_treshold & 
                                                    padj < p_treshold ~ 'Signif. down-regulated',
                                                  log2FoldChange > lfc_treshold & 
                                                    padj < p_treshold ~ 'Signif. up-regulated',
                                                  T ~ 'NS')) %>% 
  dplyr::left_join(tissue, by = c("miRname" = "miRname"))
  
  sig_df <- df %>%
    dplyr::filter(significance == 'Signif. down-regulated' | 
                    significance == 'Signif. up-regulated') %>%
    dplyr::filter(complete.cases(.))
  
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


get_CategoryExpressionPlot <- function(results, sig_log2FC, sig_pval, targets,
                                       labels = c("predicted", "validated")){
  subsets <- lapply(targets, function(x) {
    x %>% 
      dplyr::group_by(genesymbol) %>%
      dplyr::arrange(desc(baseMean.miR)) %>% 
      dplyr::distinct(genesymbol, .keep_all = T) %>%
      dplyr::rowwise() %>% 
      dplyr::mutate(pairs = paste0(genesymbol, " [", miRname, "]")) %>% 
      dplyr::select(genesymbol, pairs)
    })
  
  point.data <- dplyr::inner_join(results, subsets[["predicted"]], 
                                  by = c("geneID" = "genesymbol"),
                                  multiple = "first")
  
  label.data <- dplyr::inner_join(results, subsets[["validated"]], 
                                  by = c("geneID" = "genesymbol"),
                                  multiple = "first")
  
  plot <- ggplot(data = na.omit(results), 
         aes(x = log2FoldChange, 
             y = -log10(pvalue), 
             colour = significance)) +
    geom_point(mapping = aes(), inherit.aes = T, size = 2.5) +
    geom_point(data = point.data,
               position = "identity",
               size = 2.5, fill = "transparent",
               colour = I (alpha ("yellow", 0.6))) + 
    scale_color_manual(values = c(
      "NS" = "#c1c1c1",
      "Log10P" = '#363636',
      "Log2FoldChange" = '#767676',
      "Signif. up-regulated" = '#841f27',
      "Signif. down-regulated" = '#000f64'
    )) + 
    labs(x = expression(paste(log[2], 'FoldChange')),
         y = expression(paste(log[10], italic('P')))) + 
    theme(axis.title = element_text(size = 14), 
          axis.text = element_text(size = 14), 
          legend.position = 'none')  + 
    scale_x_continuous(expand = expansion(0.2)) + 
    geom_vline(xintercept = c(-1.5, 1.5), 
               linetype = 'dotted', size = 1) #+ 
    #geom_hline(yintercept = -log10(0.05), 
     #          linetype = 'dotted', size = 1)
  
  if (labels == "predicted") {
    plot <- plot +
      geom_label_repel(
        data = point.data, hjust = .5, vjust = .5, 
        colour = 'black', position = 'identity', 
        show.legend = F, label = point.data[,"pairs"],
        max.overlaps = 100, min.segment.length = 0, size = 3,
        arrow = arrow(type = "closed", angle = 15, length = unit(0.1,"in")),
        box.padding = unit(0.5,"in"), point.padding = unit(0.1,"in"),
        force = 1, direction = "both")
  } else {
    plot <- plot +
      geom_label_repel(
        data = label.data, hjust = .5, vjust = .5, 
        colour = 'black', position = 'identity', 
        show.legend = F, label = label.data[,"pairs"],
        max.overlaps = 100, min.segment.length = 0, size = 3,
        arrow = arrow(type = "closed", angle = 15, length = unit(0.1,"in")),
        box.padding = unit(0.5,"in"), point.padding = unit(0.1,"in"),
        force = 1, direction = "both")
  }
  
  return(plot)
}


###############################
# Overrepresentation analysis #
###############################
# I. EXTRACT GENE LISTS
get_genelist <- function(.df, .filter, .value, .name){
  # Extract the background gene list of every expressed gene
  background <- .df %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull(.value, name = .name) %>% 
    sort(., decreasing = T)
  # Extract the gene list of interest of DEGs
  interest <- .df %>%
    dplyr::filter(.filter) %>%
    dplyr::distinct(entrezID, .keep_all = T) %>%
    dplyr::pull(.value, name = .name) %>% 
    sort(., decreasing = T) 
  interest <- interest[!is.na(names(interest))]
  
  return(list(background = background, interest = interest))
}

# II. OVER-REPRESENTATION ANALYSIS
run_ora <- function(.interest, .background, .pathways){
  ora <- enricher(gene = names(.interest), # gene set of interest 
                  pvalueCutoff = 1, 
                  qvalueCutoff = 1,
                  pAdjustMethod = "BH", 
                  universe = names(.background), # background gene set
                  TERM2GENE = dplyr::select(
                    .pathways,
                    gs_name,
                    human_entrez_gene))
  ora <- setReadable(ora, org.Hs.eg.db, keyType = "ENTREZID")
  return(list("ora" = ora))
}

extract_ora_results <- function(.ora, .db){
  .db <- .db %>% 
    dplyr::select(gs_name, gs_exact_source, gs_description) %>% 
    dplyr::distinct()
  # extract data frames
  df <- as.data.frame(.ora@result)
  # order on p-value
  df <- df[order(df$`p.adjust`, decreasing = F),]
  # change gene ratio to numeric values
  df <- df %>% 
    mutate(GeneRatio = sapply(stringr::str_split(df$GeneRatio, "/"), 
                              function(x) 100*(as.numeric(x[1])/as.numeric(x[2]))),
           BgRatio = sapply(stringr::str_split(df$BgRatio, "/"), 
                            function(x) 100*(as.numeric(x[1])/as.numeric(x[2]))))
  # extract pathway IDs and description from database
  
  # add the pathway IDs and descriptions to the data frame
  ids <- df$ID
  database <- sapply(stringr::str_split(ids, "_"), 
                     function(x) return(x[1]))
  
  df$ID = .db[match(ids, .db$gs_name),][["gs_exact_source"]]
  df$Description = .db[match(ids, .db$gs_name),][["gs_description"]]
  df$Database = database
  
  # extract significant results: adjusted p-value < 0.05
  sig_df <- df %>% 
    dplyr::filter(p.adjust < 0.05)
  
  #return data frames
  return(list("df" = df, "sig_df" = sig_df))
}

# III. GENE SET ENRICHMENT ANALYSIS
run_gsea <- function(.geneset, .terms){
  set.seed(42)
  ## run the GSEA analysis
  res <- GSEA(
    geneList = .geneset, # gene set of interest (ordered on effect size)
    minGSSize = 10, # minimum size of the gene set
    maxGSSize = 500, # maximum size of the gene set
    pvalueCutoff = 1, # adjusted p-value cutoff
    eps =  0, # p-value cutoff (minimum)
    seed = TRUE, # seed for reproducibility
    pAdjustMethod = "BH", # p-value adjustment method
    TERM2GENE = dplyr::select(
      .terms,
      gs_name,
      human_entrez_gene
    ))
  res <- setReadable(res, org.Hs.eg.db, keyType = "ENTREZID")
  
  # extract data frame
  return(list("gsea" = res))
}

extract_gsea_results <- function(.gsea, .db){
  .db <- .db %>%
    dplyr::select(gs_subcat, gs_name, gs_exact_source, gs_description) %>% 
    dplyr::distinct()
  # extract data frames
  df <- as.data.frame(.gsea@result)
  # order on p-value
  df <- df[order(df$`p.adjust`, decreasing = F),]
  
  df <- df %>%
    dplyr::mutate(Direction = dplyr::case_when(NES > 0 ~ '+',
                                               NES < 0 ~ '-'))
  df <- df %>%
    dplyr::mutate(
      Name = sapply(stringr::str_split(ID, "_"), 
                    function(x) return(paste(x[-1], collapse = "_")))) %>%
    dplyr::rowwise(.) %>% 
    # Calculate background ratio
    dplyr::mutate(
      Count = length(unlist(strsplit(core_enrichment,"\\/"))),
      geneRatio = length(unlist(strsplit(core_enrichment,"\\/")))/setSize
    ) %>% 
    dplyr::ungroup(.)
  
  # add the pathway IDs and descriptions to the data frame
  ids <- df$ID
  df$ID = .db[match(ids, .db$gs_name),][["gs_exact_source"]]
  df$Description = .db[match(ids, .db$gs_name),][["gs_description"]]
  df$Database = .db[match(ids, .db$gs_name),][["gs_subcat"]]
  
  # extract significant results: adjusted p-value < 0.05
  sig_df <- df %>% 
    dplyr::filter(p.adjust < 0.05)
  
  #return data frames
  return(list("df" = df, "sig_df" = sig_df))
}


####################################
# Visualisation of gene expression #
####################################
make_pca <- function(vst, key, group){
  
  tmp <- prcomp(t(assay(vst)), center =T, scale. = TRUE)
  plot <- ggbiplot(pcobj = tmp, choices = 1:2, scale = 1,
                   groups = group, var.axes = F, circle = T) +
    geom_point(size = 3, aes(color = group)) +
    stat_ellipse(geom = "polygon", 
                 aes(color = group, fill = group), 
                 type = "norm", level = 0.68, 
                 alpha = .3, linetype = 2, show.legend = F) + 
    geom_label_repel(aes(label = vst@colData[,key], color = group), 
                     size = 5, show.legend = F) + 
    scale_color_manual(
      values = condition_palette, 
      labels = c("Control", "C.albicans (MOI 1:1)",
                 "C.parapsilosis (MOI 1:1)","C.parapsilosis (MOI 5:1)"),
      name = "Treatment groups",
      aesthetics = c("colour", "fill")
      )
  
  
  return(plot)
}

make_heatmap <- function(dds, key, contrast = "Intercept"){
  cols <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(255)
  #get annotations
  anno <- data.frame(colData(dds)[key])
  ids <-  na.omit(mapIds(org.Hs.eg.db, row.names(dds),
                 keytype = "ENSEMBL",column = "SYMBOL",
                 multiVals = "first"))
  
  res <- results(dds, name = contrast)
  res <- res[match(names(ids), row.names(res)),]
  #prepare matrix
  mat <- assay(rlog(dds))[match(names(ids), row.names(assay(dds))),]
  mat <- mat[head(order(res$stat, decreasing = T),20),]
  mat <- mat - rowMeans(mat)
  
  return(pheatmap(mat = mat,
                  color = cols,
                  labels_row = ids[match(row.names(mat), names(ids))], 
                  annotation_col = anno))
}

make_dotplot <- function(mydata, Count, type=c("GO","KEGG")){
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
  
  if(type == "KEGG"){
    data <- mydata %>% 
      #dplyr::mutate(
      #GeneRatio = as.numeric(str_split_i(GeneRatio,pattern = "/",1))/as.numeric(str_split_i(GeneRatio,pattern = "/", 2))) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::mutate(Description = fct_reorder(Description, GeneRatio)) %>%
      .[1:min(10, (dim(mydata)[1])),]
  } else {
    data <- mydata %>%
      #dplyr::mutate(
      #  GeneRatio = as.numeric(str_split_i(GeneRatio,pattern = "/",1))/as.numeric(str_split_i(GeneRatio,pattern = "/", 2))) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::mutate(Name = fct_reorder(Name, geneRatio)) %>%
      group_by(., Database) %>%
      dplyr::slice_head(n = 10)
  }
  
  
  
  if(type == "KEGG"){
    plot <-         
      ggplot(data, aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
      geom_point(aes(color = pvalue)) + 
      guides(color = guide_colorbar(title = "p-value", order = 1),
             size = guide_legend(title = "Gene count", order = 2)) +
      scale_size(range = c(2,7), limits = limits, breaks =  breaks) + 
      theme(axis.text.x = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            strip.text.y = element_text(size = 14)) + 
      labs(size = "Count", y = "")
  } else {
    plot <-         
      ggplot(data, aes(x = geneRatio, y = Name, size = as.numeric(Count))) + 
      geom_point(aes(color = pvalue)) + 
      guides(color = guide_colorbar(title = "p-value", order = 1),
             size = guide_legend(title = "Gene count", order = 2)) +
      scale_size(range = c(2,7), limits = limits, breaks =  breaks) + 
      theme(axis.text.x = element_text(size = 14),
            axis.title.x = element_text(size = 14),
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            strip.text.y = element_text(size = 14)) + 
      labs(size = "Count", y = "") + 
             facet_grid(Database ~ ., scales = "free", 
                        labeller = as_labeller(
                          c("GO:BP" = "Biological processes",
                            "GO:MF" = "Molecular functions",
                            "GO:CC" = "Cellular components"))) +
             theme(strip.text.y = element_text(size = 14))
    
  } 
  
      return(plot)
}


make_vulcanplot <- function(total, sig){
  colourPalette <- c('#c1c1c1', '#363636','#222222', '#000f64','#841f27')
  
  return((ggplot(data = na.omit(total), 
                 aes(x = log2FoldChange, 
                     y = -log10(pvalue), 
                     colour = significance)) 
          + geom_point(mapping = aes(), alpha = .5,
                       inherit.aes = T, size = 2.5) 
          + scale_color_manual(values = c(
            "NS" = "#c1c1c1",
            "Log10P" = '#363636',
            "Log2FoldChange" = '#767676',
            "Signif. up-regulated" = '#841f27',
            "Signif. down-regulated" = '#000f64'
          ))
          + labs(x = expression(paste(log[2], 'FoldChange')),
                 y = expression(paste(log[10], italic('P'))))
          + scale_x_continuous(expand = expansion(0.2))
          + geom_vline(xintercept = c(-1.5, 1.5), 
                       linetype = 'dotted', size = 1) 
          # + geom_hline(yintercept = -log10(0.05), 
          #              linetype = 'dotted', size = 1) 
          + geom_text(data = na.omit(sig), hjust = 0, vjust = 1.5, 
                      colour = 'black', position = 'identity', 
                      show.legend = F, check_overlap = T,
                      label = na.omit(sig)[,"geneID"])
          + theme_minimal()
          + theme(axis.title = element_text(size = 14), 
                  axis.text = element_text(size = 14), 
                  legend.position = 'none') ))
}

make_gseaplot <- function(data, type = NULL){
  plot <- data %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::mutate(Name = fct_reorder(Name, abs(NES))) %>%
    dplyr::slice_head(n=10) %>%
    ggplot(aes(x = abs(NES), 
               y = Name,
               size = abs(NES))) + 
    geom_point(aes(fill = Direction, shape = Direction)) + 
    guides(size = FALSE) +
    scale_fill_manual(
      values = c(
        "+" = "red",
        "-" = "blue"),
      labels = c("+" = "Positive", 
                 "-" = "Negative"),
      na.value = "grey50",
      guide = "legend",
      name = "Direction of enrichment") + 
    scale_shape_manual(
      values = c(
        "+" = 24,
        "-" = 25
      ),
      labels = c("+" = "Positive", 
                 "-" = "Negative"),
      guide = "legend",
      name = "Direction of enrichment"
    ) +  
    theme_classic() + 
    theme(axis.text.x = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14)) + 
    labs(x = "Normalized Enrichment score (NES)", y="")
  
  if(type == "GO"){
    return(plot + 
             facet_grid(Database ~ ., scales = "free", 
                        labeller = as_labeller(
                          c("GO:BP" = "Biological processes",
                            "GO:MF" = "Molecular functions",
                            "GO:CC" = "Cellular components"))) +
             theme(strip.text.y = element_text(size = 14),
                   strip.background = element_rect(fill = "white"),
                   legend.key = element_blank(),
                   panel.background = element_rect(color = "black",
                                                   linewidth = 0.5))
    )
  } else {
    return(plot)
  }
}

make_simMatrix <- function(df, type, treshold){
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
    tibble::rownames_to_column(., var = "go")
  
  reduced <- merge(reduced, pca, by = 'go')
  reduced$parentTerm <- as.factor(reduced$parentTerm)

  return(list(simM = simM, reduced = reduced))
}

# ---------------------------------------------------------------------- #
# Function - Cluster analysis                                            #
# ---------------------------------------------------------------------- #
# a function calculating the Cohen's Kappa coefficient between a list of gene sets
cohen_kappa <- function(.list){
  N <- length(.list)
  kappa_mat <- matrix(0, nrow = N, ncol = N,
                      dimnames = list(names(.list), names(.list)))
  diag(kappa_mat) <- 1
  
  total <- length(unique(unlist(.list)))
  
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      genes_i <- .list[[i]]
      genes_j <- .list[[j]]
      both <- length(intersect(genes_i, genes_j))
      term_i <- length(base::setdiff(genes_i, genes_j))
      term_j <- length(base::setdiff(genes_j, genes_i))
      no_terms <- total - sum(both, term_i, term_j)
      observed <- (both + no_terms)/total
      chance <- (both + term_i) * (both + term_j)
      chance <- chance + (term_j + no_terms) * (term_i + 
                                                  no_terms)
      chance <- chance/total^2
      kappa_mat[j, i] <- kappa_mat[i, j] <- (observed - 
                                               chance)/(1 - chance)
    }}
  return(kappa_mat)
}

# A function to create a network graph of gene sets based on the Cohen's Kappa coefficient
# then cluster the network based on interconnectivity using Louvain algorithm
get_cluster <- function(.df, .matrix, .column, .threshold = 0.25){
  # extract the list of enriched gene sets 
  genes <- split(strsplit(.df[["core_enrichment"]],"/"), .df[["ID"]])
  
  # subset the cohen's matrix for the gene sets of interest
  similarity.matrix <- .matrix[.df[['ID']], .df[['ID']]]
  
  # create an adjacency matrix based on the similarity matrix
  adjacency.matrix <- similarity.matrix > .threshold
  
  set.seed(42)
  # create a network graph from the adjacency matrix
  graph <- graph_from_adjacency_matrix(adjacency.matrix, mode = "undirected", diag = F)
  
  # annotate the nodes of the graph with the enrichment results
  V(graph)$Name <- .df$Name
  V(graph)$NES <- .df$NES
  V(graph)$geneRatio <- .df$geneRatio
  
  # Community detection (Louvain algorithm)
  clusters <- cluster_louvain(graph)
  # annotate nodes with cluster membership
  V(graph)$cluster <- membership(clusters)
  
  # extract cluster information
  cluster_summary <- data.frame(
    ID = names(V(graph)),
    cluster = as.factor(V(graph)$cluster))
  
  cluster_summary <- distinct(cluster_summary, .keep_all = T)
  
  df = .df %>% 
    dplyr::inner_join(., cluster_summary, by = "ID")
  
  return(list(graph = graph, df = df))
}

get_cluster_representative <- function(.cluster, .degs){
  # add gene expression values to the clusters
  linkage <- .cluster %>% 
    dplyr::select(ID, Name, core_enrichment, cluster) %>%
    tidyr::separate_rows(core_enrichment, sep = "/")
  
  linkage$logFC = .degs$log2FoldChange[match(linkage$core_enrichment, .degs$geneID)]
  linkage$padj = .degs$padj[match(linkage$core_enrichment, .degs$geneID)]
  
  # calculate normalized term weight
  linkage$weight = abs(linkage$logFC)*(-log10(linkage$padj))
  linkage$weight = linkage$weight/max(linkage$weight, na.rm = T)
  
  linkage = linkage %>% 
    dplyr::select(c("core_enrichment", "ID", "weight", "cluster")) %>%
    setNames(.,c("node1", "node2", "weight", "cluster")) %>% 
    as.data.frame(.)
  
  # create a network visualization of gene and GO-term relationships
  net <- graph_from_data_frame(linkage)
  net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)
  
  # calculate hub score of each gene
  hs <-  igraph::hub_score(net, scale = T, weights = linkage$weight)$vector
  
  # summarize gene hub scores, to determine the final weight of each GO term
  linkage$geneRatio = .cluster$geneRatio[match(linkage$node2, .cluster$ID)]
  
  linkage <- linkage %>%
    dplyr::mutate(hub_score = geneRatio * hs[match(node1, names(hs))]) %>% 
    dplyr::rename(geneID = node1, ID = node2) %>% 
    dplyr::group_by(ID, cluster) %>%
    dplyr::summarise(core_enrichment = paste(geneID, collapse = "/"),
                     hub_score = sum(hub_score, na.rm = T))
  
  out_df <- inner_join(.cluster, linkage, by = c("ID","core_enrichment","cluster")) %>%
    dplyr::mutate(cluster = as.factor(cluster))
  
  # select cluster representatives according to calculated weight
  representative.terms <- out_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(desc(hub_score)) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::pull(ID)
  
  out_df <- out_df %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(Representative = ifelse(ID %in% representative.terms, T, F))
  
  return(out_df)
}

filter_graph <- function(.graph, .threshold){
  cl <- table(V(.graph)$cluster)
  
  valid_cl <- which(cl >= .threshold)
  
  filtered_vertices <- which(V(.graph)$cluster %in% valid_cl)
  
  subgraph <- induced_subgraph(.graph, vids = filtered_vertices)
  
  return(subgraph)
}

plot_network <- function(.net, .layout, .labels, .df){
  edges = as.data.frame(as_edgelist(.net)) %>% 
    setNames(c("from", "to"))
  
  df = data.frame(x = as.data.frame(.layout)[,1],
                  y = as.data.frame(.layout)[,2],
                  ID = names(V(.net)),
                  Name = V(.net)$Description,
                  cluster = as.factor(V(.net)$cluster),
                  Representative = V(.net)$Representative,
                  weight = abs(.df[["hub_score"]][match(names(V(.net)), .df[["ID"]])]),
                  regulation = ifelse(.df[["NES"]][match(names(V(.net)), .df[["ID"]])] > 0, "UP", "DOWN")
  )
  
  lines = edges %>% 
    dplyr::mutate(
      from.x = df$x[match(edges$from, df$ID)],
      from.y = df$y[match(edges$from, df$ID)],
      to.x = df$x[match(edges$to, df$ID)],
      to.y = df$y[match(edges$to, df$ID)]
    )
  
  plot = ggplot() +
    stat_ellipse(data = df, geom = "polygon", 
                 aes(x = x, y = y, group = cluster), fill = "grey75", color = "grey15", #color = cluster, fill = cluster), 
                 type = "norm", level = 0.9,
                 alpha = .1, linetype = 2, show.legend = F) +
    geom_segment(data = lines, aes(x = from.x, y = from.y, xend = to.x, yend = to.y),
                 color = "grey25") +
    geom_point(data = df, aes(x = x, y = y, size = weight, fill = regulation), #fill = cluster, size = weight),
               shape = 21, colour = "black") +
    # geom_point(data = subset(df, Representative), aes(x = x, y = y), #fill = cluster, size = 5),
    #            shape = 21, colour = "orange", fill = "transparent", stroke = 2, size = 15) +
    geom_label_repel(data = subset(df, Name %in% .labels), 
                     aes(x = x, y = y,
                         label = ifelse(Representative, str_wrap(gsub("_"," ",Name), 30), "")),
                     max.overlaps = 100, min.segment.length = 0, size = 6,
                     color = "grey15", fill = "white", show.legend = F,
                     arrow = arrow(type = "closed", angle = 15, length = unit(0.1,"in")),
                     box.padding = unit(0.5,"in"), point.padding = unit(0.1,"in"),
                     force = 1, direction = "both") +
    scale_fill_manual(values = c("UP" = "darkred", "DOWN" = "blue"),
                      labels = c("UP" = "Activated", "DOWN" = "Inhibited"),
                      name = c("Regulation"),
                      guide = guide_legend(override.aes = c(size = 10))) +
    scale_size_continuous(guide = "none", range=c(5,15)) +
    theme_void() +
    theme(legend.title=element_text(size=28), 
          legend.text=element_text(size=28),
          legend.spacing=unit(1.5,"lines"),
          legend.position = "right")
  
  return(plot)
}

getCircplotData <- function(.cluster, .deg, .interest_cluster, .interest_cluster_genes, .palette){
  linkage <- .cluster %>% 
    dplyr::filter(ID %in% names(.interest_cluster)) %>% 
    dplyr::select(ID, Name, core_enrichment, cluster) %>%
    tidyr::separate_rows(core_enrichment, sep = "/") %>% 
    dplyr::left_join(.deg[,c("geneID","log2FoldChange","pvalue")],
                     by = c("core_enrichment" = "geneID")) %>% 
    dplyr::rowwise() %>%
    # ...as a function of the log2FoldChange and adjusted p-value
    dplyr::mutate(weight = abs(log2FoldChange)*(-log10(pvalue))) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(complete.cases(.)) %>% 
    dplyr::mutate(weight = as.numeric(weight/max(weight))) %>%
    dplyr::select(c("core_enrichment", "Name", "weight", "cluster")) %>%
    setNames(.,c("node1", "node2", "weight", "cluster")) %>%
    as.data.frame(.)
  
  # create a network visualization of gene and GO-term relationships
  net <- graph_from_data_frame(linkage)
  net <- igraph::simplify(net, remove.multiple = F, remove.loops = T)
  # calculate hub score of each gene
  hs <-  igraph::hub_score(net, scale = T, weights = linkage$weight)$vector
  
  # summarize gene hub scores, to determine the final weight of each GO term
  linkage <- linkage %>%
    dplyr::rowwise(.) %>%
    dplyr::mutate(hub_score = hs[which(names(hs) == node1)])
  
  linkage_wide <- linkage %>% 
    dplyr::select(node1, node2, hub_score) %>%
    tidyr::pivot_wider(names_from = node2, values_from = hub_score) %>%
    dplyr::mutate(across(where(is.numeric), ~ifelse(is.na(.), 0, .))) %>% 
    tibble::column_to_rownames("node1") %>% 
    as.matrix()
  
  background = unique(linkage$node1[-c(which(linkage$node1 %in% .interest_cluster_genes))])
  focus = .interest_cluster_genes                   
  grid.col = c(.palette,
               setNames(rep("grey", length(background)), background),
               setNames(rep("red", length(focus)), focus))
  
  # Extract sector names
  from_sectors <- rownames(linkage_wide)
  to_sectors <- colnames(linkage_wide)
  
  border_mat = matrix(NA, nrow = nrow(linkage_wide), ncol = ncol(linkage_wide))
  border_mat[row.names(linkage_wide) %in% .interest_cluster_genes, ] = "red"
  dimnames(border_mat) = dimnames(linkage_wide)
  
  return(list(data.mat = linkage_wide, border.mat = border_mat, grid.col = grid.col))
}

plotCircplot <- function(.path, .data, .color, .links, .labels){
  png(.path, 
      width = 20, height = 12, res = 300, units = "in")
  par(mfrow=c(1,1),mar=c(2,4,4,2), xpd = TRUE)
  
  circos.par(start.degree = 90)
  circlize::chordDiagram(.data, annotationTrack = "grid",
                         grid.col = .color, transparency = 0.5,
                         link.border = .links,
                         preAllocateTracks = list(track.height = 0.05))
  
  # Add labels only to the selected sectors
  circos.track(track.index = 1, panel.fun = function(x, y) {
    sector_name <- get.cell.meta.data("sector.index")
    
    # Label all 'to' sectors and only selected 'from' sectors
    if (sector_name %in% .labels) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1],
                  str_wrap(gsub("_", " ", CELL_META$sector.index), 25), 
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
    }
  }, bg.border = NA)
  dev.off()
  circos.clear()
}

make_GO_simplot <- function(.reduced){
  return(ggplot( data = .reduced )
         + stat_ellipse(geom = "polygon", aes( x = x, y = y, colour = cluster, fill = cluster), 
                        type = "norm", level = 0.68, 
                        alpha = .15, linetype = 2, show.legend = F)
         + geom_point( aes( x = x, y = y, colour = cluster, size = hub_score), 
                       alpha = I(0.6), show.legend = F)
         + geom_point( aes( x = x, y = y, size = hub_score), shape = 21, 
                       fill = "transparent", colour = I (alpha ("black", 0.3) ), show.legend = F)
         + scale_size( range=c(5, 30))
         + theme_bw()
         + geom_label_repel( data = subset.data.frame(.reduced, Representative),
                             aes(x = x, y = y, label = Name, colour = cluster), 
                             size = 3, show.legend = F )
         + guides( color = "none")
         + labs (y = "semantic space x", x = "semantic space y")
         + theme(legend.key = element_blank()))
}


##############################
# Venn diagram & upset plots #
##############################
make_upsetbase <- function(table, vars){
  table <- table %>%
    dplyr::mutate(var1 = ifelse(is.na(logFC_A), FALSE, TRUE)) %>%
    dplyr::mutate(var2 = ifelse(is.na(logFC_B), FALSE, TRUE)) %>%
    dplyr::mutate(var3 = ifelse(is.na(logFC_C), FALSE, TRUE)) %>%
    relocate(where(is.logical), .before = where(is.character))
  
  colnames(table)[1:3] <- vars
  
  return(table)
}

make_upsetvennbase <- function(arranged, df){
  tmp = data.frame(
    region = c("CA11","CA11-CP51","CA11-CP51-PANC","CA11-PANC","CP51","CP51-PANC","PANC"),
    x = c(-1.3, 0, 0, -0.65, 1.3, 0.65, 0),
    y = c(1.0392, 1.0392, 0.2887, -0.0866, 1.0392, -0.0866, -1.2124)
  )
  tmp = merge(tmp, arranged[,c("region","size")], by = "region", all = T)
  
  set.seed(42)
  
  xy = rbind(
    calculateCircle(x = tmp$x[1], y = tmp$y[1], r = 0.2, 
                    noiseFun = function(x) (x + rnorm(1,0,0.2)),
                    steps = ifelse(is.na(tmp$size[1]),0, tmp$size[1]), 
                    randomDist = T, randomFun = rnorm),
    calculateEllipse(x = tmp$x[2], y = tmp$y[2], a = 0.1, b = .1, angle = 0, 
                     noiseFun = function(x) (x + rnorm(1,0,0.1)),
                     steps = ifelse(is.na(tmp$size[2]),0, tmp$size[2]), 
                     randomDist = T, randomFun = rnorm),
    calculateCircle(x = tmp$x[3], y = tmp$y[3], r = .1, 
                    noiseFun = function(x) (x + rnorm(1,0,0.1)),
                    steps = ifelse(is.na(tmp$size[3]),0, tmp$size[3]), 
                    randomDist = T, randomFun = rnorm),
    calculateEllipse(x = tmp$x[4], y = tmp$y[4], a = 0.1, b = .1, angle = -120, 
                     noiseFun = function(x) (x + rnorm(1,0,0.1)),
                     steps = ifelse(is.na(tmp$size[4]),0, tmp$size[4]), 
                     randomDist = T, randomFun = rnorm),
    calculateCircle(x = tmp$x[5], y = tmp$y[5], r = .2,  
                    noiseFun = function(x) (x + rnorm(1,0,0.2)),
                    steps = ifelse(is.na(tmp$size[5]),0, tmp$size[5]), 
                    randomDist = T, randomFun = rnorm),
    calculateEllipse(x = tmp$x[6], y = tmp$y[6], a = 0.1, b = .1, angle = -240, 
                     noiseFun = function(x) (x + rnorm(1,0,0.1)),
                     steps = ifelse(is.na(tmp$size[6]),0, tmp$size[6]), 
                     randomDist = T, randomFun = rnorm),
    calculateCircle(x = tmp$x[7], y = tmp$y[7], r = .2, 
                    noiseFun = function(x) (x + rnorm(1,0,0.2)),
                    steps = ifelse(is.na(tmp$size[7]),0, tmp$size[7]), 
                    randomDist = T, randomFun = rnorm)
  )
  
  base <- df %>%
    dplyr::mutate(
      region = dplyr::case_when(
        CA11 & CP51 & PANC ~ "CA11-CP51-PANC",
        CA11 & !CP51 & PANC ~ "CA11-PANC",
        !CA11 & CP51 & PANC ~ "CP51-PANC",
        CA11 & CP51 & !PANC ~ "CA11-CP51",
        !CA11 & !CP51 & PANC ~ "PANC",
        !CA11 & CP51 & !PANC ~ "CP51",
        CA11 & !CP51 & !PANC ~ "CA11"),
      region = as.factor(region)) %>%
    dplyr::arrange(region)
  
  base <- base %>%
    dplyr::mutate(
      x = xy[,1],
      y = xy[,2]) %>%
    dplyr::select(c("CA11","CP51","PANC","region","x","y","Trend","geneID"))
  return(base)
}

make_upsetvenn <- function(data, sets, labels){
  panc_label = labels
  
  plot <- ggplot(data) + 
    coord_fixed(clip = "off") +
    theme_void() + 
    scale_color_manual(
      name = "Regulation",
      values = regulation_palette,
      na.value = "grey",
      aesthetics = c("colour")
    )  + 
    geom_venn_region(data, sets = sets,
                     alpha = 0.1, show.legend = F) + 
    geom_venn_circle(data, sets = sets,
                     size = .5, color = "black") + 
    geom_point(aes(x = x, y = y, colour = Trend), size = 2)  + 
    scale_fill_venn_mix(data, sets = sets, colors = condition_palette) + 
    geom_venn_label_set(data, outwards_adjust = 2.5,
                        sets = sets, size = 8, label.padding = unit(0.7, "lines"),
                        fill = alpha("white", .9), 
                        aes(label = panc_label)) + 
    geom_venn_label_region(
      data, sets = sets, size = 6,
      aes(label=size), 
      outwards_adjust=1.5, 
      position=position_nudge(y=0)) + 
    theme(
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14))
  
  return(plot)
}

make_quad_venn <- function(set1, set2, names, title, dir){
  list <- list(S1 = set1 %>%
                 dplyr::filter(region %in% c("CA11-PANC","CA11-CP51-PANC")) %>%
                 dplyr::pull(geneID) %>%
                 as.character(.),
               S2 = set1 %>%
                 dplyr::filter(region %in% c("CP51-PANC","CA11-CP51-PANC")) %>%
                 dplyr::pull(geneID) %>%
                 as.character(.),
               S3 = set2 %>%
                 dplyr::filter(region %in% c("CA11-PANC","CA11-CP51-PANC")) %>%
                 dplyr::pull(geneID),
               S4 = set2 %>%
                 dplyr::filter(region %in% c("CP51-PANC","CA11-CP51-PANC")) %>%
                 dplyr::pull(geneID))
  names(list) = names
  
  (plot <- venn_plot(list,
                     show_elements = TRUE,
                     label_sep = "\n",
                     text_size = 3))
  
  venn <- nVennR::plotVenn(list,
                           setColors = c("#8B0000","#045D5D","#FFA500","#00FFFF"),
                           borderWidth=2, opacity=0.2,
                           outFile = paste0(dir, "/",title, ".svg"))
  
  genes <- listVennRegions(venn)
  
  genes <- lapply(names(genes), function(x){
    tmp <- genes[[x]]
    df <- tmp %>% data.frame(genes = tmp, locate = x)
  })
  
  genes <- do.call(rbind, genes) %>%
    dplyr::mutate(region = dplyr::case_when(
      stringr::str_detect(locate, "1, 0, 0, 0") ~ names[1],
      stringr::str_detect(locate, "0, 1, 0, 0") ~ names[2],
      stringr::str_detect(locate, "0, 0, 1, 0") ~ names[3],
      stringr::str_detect(locate, "0, 0, 0, 1") ~ names[4],
      stringr::str_detect(locate, "1, 1, 0, 0") ~ paste(names[c(1,2)],collapse = ", "),
      stringr::str_detect(locate, "1, 0, 1, 0") ~ paste(names[c(1,3)],collapse = ", "),
      stringr::str_detect(locate, "1, 0, 0, 1") ~ paste(names[c(1,4)],collapse = ", "),
      stringr::str_detect(locate, "0, 1, 1, 0") ~ paste(names[c(2,3)],collapse = ", "),
      stringr::str_detect(locate, "0, 1, 0, 1") ~ paste(names[c(2,4)],collapse = ", "),
      stringr::str_detect(locate, "0, 0, 1, 1") ~ paste(names[c(3,4)],collapse = ", "),
      stringr::str_detect(locate, "1, 1, 1, 0") ~ paste(names[c(1,2,3)],collapse = ", "),
      stringr::str_detect(locate, "1, 1, 0, 1") ~ paste(names[c(1,2,4)],collapse = ", "),
      stringr::str_detect(locate, "1, 0, 1, 1") ~ paste(names[c(1,3,4)],collapse = ", "),
      stringr::str_detect(locate, "0, 1, 1, 1") ~ paste(names[c(2,3,4)],collapse = ", "),
      stringr::str_detect(locate, "1, 1, 1, 1") ~ paste(names,collapse = ", "))) %>%
    dplyr::select(c("genes","region"))
  
  return(list(plot = plot, df = genes))
}

make_venn_df <- function(df, names){
  df %>%
    merge(.,Ca11_6h_df$sig_df[c("geneID","log2FoldChange")], 
          by.x = "genes", by.y = "geneID", all.x = T) %>%
    merge(.,Cp51_6h_df$sig_df[c("geneID","log2FoldChange")], 
          by.x = "genes", by.y = "geneID", all.x = T, 
          suffixes=c(".1", ".2")) %>%
    merge(.,OKF6_Ca11_6h[c("geneID","log2FoldChange")], 
          by.x = "genes", by.y = "geneID", all.x = T) %>%
    merge(.,OKF6_Cp51_6h[c("geneID","log2FoldChange")], 
          by.x = "genes", by.y = "geneID", all.x = T, 
          suffixes=c(".3", ".4")) %>%
    setNames(c("geneID","region",names))
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
  (ggplot(df,
          aes(
            x = reorder(str_wrap(category, 5), geneRation),
            y = geneRation,
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
       data = subset.data.frame(df, geneRation > 0),
       aes(label = str_wrap(paste(paste0("zsc:",round(zscore,2)),
                                  paste0("(",geneRation,"%)")), 10),
           y = max(geneRation)-5), direction = "y",
       fill="white",label.padding = unit(5, "pt"),
       position=position_dodge(width=0.9), vjust=-0.25) +
     scale_fill_gradient2(
       "Activation (z-score)",
       low = "#000f64", mid = muted("#000f64"), high = "#841f27", 
       midpoint = 0,
       aesthetics = c("fill","color")
     ) + 
     facet_grid(.~condition, 
                labeller = as_labeller(
                  c("CA11" = "C.albicans (MOI 1:1)",
                    "CP51" = "C.parapsilosis (MOI 5:1)")
                )) +
     theme_bw() + 
     theme(
       # Remove axis ticks and text
       line = element_line(color = "gray80", linewidth = 1),
       axis.title = element_blank(),
       axis.ticks = element_blank(),
       axis.text.y = element_blank(),
       # Use gray text for the region names
       axis.text.x = element_text(color = "gray12", size = 10),
       panel.background = element_rect(fill = "white",
                                       color = "white"),
       #strip parameters
       strip.clip = "off", 
       strip.text = element_text(color = "gray12", size = 14),
       # Move the legend to the bottom
       legend.position = "bottom",
     ) + 
     coord_polar())
}

###########################
## microRNA visualisation #
###########################
miRNA_plot <- function(data){
  data %>%
    dplyr::arrange(desc(abs(log2FoldChange))) %>%
    dplyr::mutate(miRname = fct_reorder(miRname, log2FoldChange)) %>%
    ggplot(aes(x = log2FoldChange, y = miRname)) + 
    geom_col(width = .5, aes(fill = significance)) +
    scale_fill_manual(values = regulation_palette,
    name = "Expression change:",
    guide = "legend") +
    theme_classic() +
    geom_errorbar(aes(y = miRname, xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), 
                  colour = 'black', position = position_dodge(.1), width = .3) + 
    labs( y = '',
          x = expression(paste(log[2], 'FoldChange'))) +
    geom_vline(xintercept = c(0), linetype = 'solid', size = 1) + 
    geom_vline(xintercept = c(-1, 1), linetype = 'dotted', size = 1) + 
    theme(axis.text = element_text(size = 16, face = 'bold', colour = 'grey25'),
          axis.title = element_text(size = 16, face = 'bold', colour = 'black'),
          legend.title = element_text(size = 14, face = 'bold', colour = 'black'),
          legend.text = element_text(size = 14, colour = 'black'),
          panel.spacing = unit(2, "cm"),
          legend.position = 'bottom',
          plot.margin = margin(1,1.5,0,0, "cm"))
}


miRNA_target_plot <- function(data, targets){
  subset <- subset.data.frame(data, subset = mirbaseId %in% targets)
  
  return(data %>%
           dplyr::arrange(desc(abs(log2FoldChange))) %>%
           dplyr::mutate(mirbaseId = fct_reorder(mirbaseId, log2FoldChange)) %>%
           ggplot(aes(x = log2FoldChange, y = mirbaseId)) + 
           geom_col(width = .5, aes(fill = regulation)) +
           geom_col(data = subset, aes(x = log2FoldChange, y = mirbaseId),
                    width = .6, fill = "transparent", linewidth = 2, colour = I (alpha ("yellow", 0.6))) + 
           scale_fill_manual(values = c(
             "Signif. up-regulated" = '#841f27',
             "Signif. down-regulated" = '#000f64'
           ),
           name = "Expression change:",
           guide = "legend") +
           theme_classic() +
           geom_errorbar(aes(y = mirbaseId, xmin = log2FoldChange - lfcSE, xmax = log2FoldChange + lfcSE), 
                         colour = 'black', position = position_dodge(.1), width = .3) + 
           labs( y = '',
                 x = expression(paste(log[2], 'FoldChange'))) +
           geom_vline(xintercept = c(0), linetype = 'solid', size = 1) + 
           geom_vline(xintercept = c(-1, 1), linetype = 'dotted', size = 1) + 
           theme(axis.text = element_text(size = 16, face = 'bold', colour = 'grey25'),
                 axis.title = element_text(size = 16, face = 'bold', colour = 'black'),
                 legend.title = element_text(size = 14, face = 'bold', colour = 'black'),
                 legend.text = element_text(size = 14, colour = 'black'),
                 legend.position = 'bottom'))
}

plot_miRNA_targets <- function(miRNA_data) {
  ggplot(miRNA_data, aes(x = genesymbol, y = log2FC)) +
    geom_bar(aes(fill = Type), show.legend = F,
             stat = "identity", position = position_dodge(width = 0.7), 
             width = 0.6, color = "black") +
    scale_fill_manual(values = c("miRNA" = "purple", "mRNA" = "steelblue")) +
    facet_wrap(~miRname,scales = "free_x") +
    labs(title = "",
         x = "", y = expression(paste(log[2], 'FoldChange'))) +
    geom_hline(yintercept = c(0), 
               linetype = 'solid', size = .5) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 12),
          strip.background = element_blank(),
          strip.text = element_text(size = 12, colour = "black"),
          panel.border = element_rect(colour = "black", fill = NA, size = 1))
}