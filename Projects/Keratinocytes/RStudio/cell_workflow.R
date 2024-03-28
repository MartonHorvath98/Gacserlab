################################################################################
# 1.) Set up the work environment and the directory structure                  #
################################################################################

# Set working directory
setwd(dir = choose.dir(getwd()))
# Set downstream path
human.folder <- "Cells"
if (!dir.exists(human.folder)) {
  dir.create(human.folder) # create the main results folder
}
date <- format(Sys.Date(), "%Y-%m-%d") # get the current date
if (!dir.exists(file.path(human.folder,date))) {
  dir.create(file.path(human.folder, date)) # create the dated results human.folder
}

# Create the sub human.folders for: results, data, and pictures
data_dir <- "data" # general for every results (only created once)
data_needed <- FALSE # flag to check if the data human.folder is needed
if (!dir.exists(file.path(human.folder,data_dir))) {
  data_needed <- TRUE
  dir.create(file.path(human.folder,data_dir)) # create the data human.folder
}

#plots directory
plots_dir <- "plots"
if (!dir.exists(file.path(human.folder, date, plots_dir))) {
  dir.create(file.path(human.folder, date, plots_dir)) # create the plots human.folder
}

#results directory
results_dir <- "tables" 
if (!dir.exists(file.path(human.folder, date, results_dir))) {
  dir.create(file.path(human.folder, date, results_dir)) # create the results human.folder
}

################################################################################
# 2.) Load the data and the required files                                     #
################################################################################
# Load the user-defined functions
source("./human.packages.R")
source("./R_functions.R")
cat(crayon::white$bold("Loading the data and the required files\n"))

if (!file.exists(
  paste(file.path(human.folder, data_dir), "total_readcounts.xlsx", sep = "/")
)) {
  # EXTRACT ANNOTATIONS
  cat(crayon::white("-> Extracting annotations from the GFF3 file...\n"))
  # Extracting transcript annotations from the GFF3 file into a TxDb object
  human.txdb <- makeTxDbFromGFF("../reference/Homo_sapiens.GRCh38.96.gff3", 
                                dataSource = "Ensembl", 
                                organism = "Homo sapiens")
  # Extracting gene features from the TxDb object
  human.genes <- exonsBy(human.txdb, by = "gene")
  
  # READ IN BAM ALIGNMENT FILES
  cat(crayon::white("-> Extracting read counts from the BAM files...\n"))
  human.path  <- choose.dir(getwd(), "Select the directory containing the '.bam' files")
  human.files <- list.files(human.path, pattern = ".bam$", full.names = T)
  human.list <- BamFileList(human.files, yieldSize = 2000000)
  
  # CALCULATE READCOUNTS
  cat(crayon::white("-> Calculating read counts..."))
  # Set up parallel computing
  register(SnowParam())
  # Calculate read counts
  human.reads <- summarizeOverlaps(features = human.genes, 
                                   reads = human.list, 
                                   # (default) count reads overlapping single exon
                                   mode = "Union",  
                                   # use strand-specificity
                                   ignore.strand = F, 
                                   # paired-end reads
                                   singleEnd = F, 
                                   fragments = T, 
                                   # strand-specificity => first strand, RF
                                   preprocess.reads = invertStrand)

  # CREATE DATA FRAME WITH READ COUNTS
  cat(crayon::white("-> creating data frame with readcounts...\n"))
  # Extract count matrix
  human.counts <- assay(human.reads)
  human.counts <- human.counts %>%
    as.data.frame() %>%
    # rename columns
    setNames(str_remove(colnames(.), ".bam"))
  # Save the read counts to an excel file
  write.xlsx(human.counts, 
             paste(file.path(human.folder, data_dir),"total_readcounts.xlsx", sep = "/"),
             keepNA = T, na.string = "NA", colNames = T, rowNames = T)
}

################################################################################
# 3.) Perform differential expression analysis                                 #
################################################################################
# Read in the (newly) created count table for further analysis
human.countTable <- read.xlsx(paste(file.path(human.folder, data_dir),"total_readcounts.xlsx", sep = "/"),
                             colNames = T, rowNames = T)

cat(crayon::white$bold("Differential expression analysis\n"))
run_DE <- user_answer("differential expression analysis")

if (run_DE){
  # CREATE THE DESIGN MATRIX
  cat(crayon::white("-> creating the design matrix...\n"))
  # get sample names
  human.samples <- colnames(human.countTable)
  human.coldata <- data.frame("samplenames" = human.samples) %>%
    # extract the cell line from the sample names
    dplyr::mutate(cells = dplyr::case_when(
      stringr::str_detect(human.samples, "HC_") ~ "HaCaT",
      stringr::str_detect(human.samples, "HK_") ~ "HPVKER"
    ),
    cells = factor(cells, levels = c("HaCaT","HPVKER"))) %>%
    # extract the treatment from the sample names
    dplyr::mutate(samplenames = as.factor(samplenames),
                  treatment = dplyr::case_when(
                    stringr::str_detect(human.samples, "_Ca_") ~ "C.albi",
                    stringr::str_detect(human.samples, "_Cp_") ~ "C.para",
                    TRUE ~ "ctrl"),
                  treatment = factor(treatment, levels = c("ctrl","C.albi","C.para"))) %>% 
    # extract the replicates from the sample names
    dplyr::mutate(rep = dplyr::case_when(
      stringr::str_detect(human.samples, "_1$") ~ "run1",
      stringr::str_detect(human.samples, "_2$") ~ "run2",
      stringr::str_detect(human.samples, "_3$") ~ "run3"
    ),
    rep = factor(rep, levels = c("run1","run2","run3"))) %>%
    # create the experiment and full setup
    dplyr::mutate(condition = as.factor(str_glue("{cells}_{treatment}",
                                                  cells = cells,
                                                  treatment = treatment)))

  # PERFORM THE DIFFERENTIAL EXPRESSION ANALYSIS
  cat(crayon::white("-> performing the differential expression analysis...\n"))
  # Effects of C.albi and C.para infection on HaCaT
  hacat.deseq <- make_deseq(matrix = human.countTable[,1:9],
                            coldata = human.coldata[1:9,],
                            design = "treatment")
  # Effects of C.albi and C.para infection on HPV-KER
  hpvker.deseq <- make_deseq(matrix = human.countTable[,10:18],
                            coldata = human.coldata[10:18,],
                            design = "treatment")
  
  # GET RESULTS
  cat(crayon::white("-> Calculating significant the results...\n"))
  # Significance thresholds: log2FC > 1.5, p-value < 0.05
  # 1. HaCaT response to Candida albicans
  hacat_vs_ca.res <- get_results(hacat.deseq$dds, 1.5, 0.05, 
                                  contrast = c("treatment","C.albi","ctrl"), 
                                  name = c(''))
  # 2. HaCaT response to Candida parapsilosis
  hacat_vs_cp.res <- get_results(hacat.deseq$dds, 1.5, 0.05, 
                                  contrast = c("treatment","C.para","ctrl"), 
                                  name = c(''))
  # 3. HPV-KER response to Candida albicans
  hpvker_vs_ca.res <- get_results(hpvker.deseq$dds, 1.5, 0.05, 
                                   contrast = c("treatment","C.albi","ctrl"), 
                                   name = c(''))
  # 4. HPV-KER response to Candida parapsilosis
  hpvker_vs_cp.res <- get_results(hpvker.deseq$dds, 1.5, 0.05, 
                                   contrast = c("treatment","C.para","ctrl"), 
                                   name = c(''))
  
  # EXPORT THE RESULTS
  cat(crayon::white("-> exporting the results...\n"))
  # Save the list of names for the results
  names <- c("CTRL_vs_C.albi_HaCat", "CTRL_vs_C.para_HaCat",
             "CTRL_vs_C.albi_HPV-KER", "CTRL_vs_C.para_HPV-KER")
  # Save the results to a list
  human.results <- list(hacat_vs_ca.res$df, hacat_vs_cp.res$df,
                    hpvker_vs_ca.res$df, hpvker_vs_cp.res$df)
  # Save the significant results to a list
  human.sig.results <- list(hacat_vs_ca.res$sig_df, hacat_vs_cp.res$sig_df,
                            hpvker_vs_ca.res$sig_df, hpvker_vs_cp.res$sig_df)
  # Add gene symbol and entrezID annotation to the results
  human.results <- lapply(human.results, annotate_results)
  human.sig.results <- lapply(human.sig.results, annotate_results)
  # names the results
  names(human.results) <- names
  names(human.sig.results) <- names
  
  # Export the results to an excel file
  sapply(names(human.results), function(x){
    openxlsx::write.xlsx(
      human.res[[x]],
      paste(file.path(human.folder,date,results_dir),
            paste0(x,".xlsx"), sep = "/"))})
  sapply(names(human.sig.results), function(x){
    openxlsx::write.xlsx(
      human.res[[x]],
      paste(file.path(human.folder,date,results_dir),
            paste0(x,"_sigGene.xlsx"), sep = "/"))})
}


################################################################################
# 4.) Perform functional analyses                                              #
################################################################################
cat(crayon::white$bold("Functional analysis\n"))
run_HM <- user_answer("heatmap of DEGs")

### 4.1 Create a heatmap of the differentially expressed genes
if (run_HM){
  cat(crayon::white("-> Performing cluster analysis on the differentially expressed genes...\n"))
  # Create a matrix of the differentially expressed genes
  human.log.matrix <- rlog(assay(human.reads))
  human.lognorm.matrix <- make_matrix(human.log.matrix,
                                      unique(
                                        unlist(
                                          lapply(human.sig.results, rownames))))
  
  # Find optimal number of clusters
  (fviz_nbclust(human.lognorm.matrix, kmeans, method = "silhouette"))
  
  # CREATE THE HEATMAP
  cat(crayon::white("-> Creating the heatmap...\n"))
  # Create the heatmap annotation
  human.annot <- HeatmapAnnotation(
    Cells = as.factor(human.coldata$cells),
    Treament = as.factor(human.coldata$treatment),
    show_annotation_name = F) 
  
  # Create an expression heatmap
  human.heatmap <- expression_heatmap(mat = human.lognorm.matrix, clusters = 4, columns = 2,
                              coldata = human.coldata, annot = human.annot,
                              method1 = "euclidean", method2 = "ward.D")
  # Create a fold-change heatmap
  human.expr <- lapply(human.sig.results, 
                       function(x){x %>% select(ensembl, log2FoldChange)}) %>% 
    merge.rec(., by="ensembl", all=T, suffixes=c("", "")) %>% 
    column_to_rownames("ensembl") %>%
    setNames(names)
  human.expr.matrix <- logFC_heatmap(mat = as.matrix(human.expr), 
                                     dend = human.heatmap$dendrogram)

  # SAVE THE HEATMAP
  cat(crayon::white("-> Saving the heatmap...\n"))
  # Save the heatmap to a file
  png(file.path(human.folder, date, plots_dir, "full_heatmap.png"), 
      width = 14, height = 10, units = 'in', res = 300)
  draw(human.heatmap$heatmap + human.expr.matrix, merge_legend = T)
  dev.off()
}

# PLOT MEAN-ABUNDANCE PLOTS
cat(crayon::white("-> Creating the mean-abundance plots...\n"))
# Create the mean-abundance plots
human.MA <- list()
for (i in names(human.results)) {
  human.MA[[i]] <- MA_plotting(human.results[[i]])
  # Save the plot to a file
  ggsave(file.path(human.folder, date, plots_dir, paste0(i,"_MAplot.png")), 
                plot = human.MA[[i]], width = 10, height = 6, units = 'in')
}
# PLOT VOLCANO PLOTS
cat(crayon::white("-> Creating the volcano plots...\n"))
# Create the volcano plots
human.expr_plot <- list()
for (i in names(human.results)) {
  human.expr_plot[[i]] <- vulcan_plotting(human.results[[i]])
  # Save the plot to a file
  ggsave(file.path(human.folder, date, plots_dir, paste0(i,"vulcanplot.png")), 
         plot = human.expr_plot[[i]], width = 12, height = 8, units = 'in')
}

# Create Venn diagram of shared differentially expressed genes 
human.venn <- make_venn(human.sig.results$CTRL_vs_C.albi_HaCat, human.sig.results$CTRL_vs_C.para_HaCat,
                        human.sig.results$`CTRL_vs_C.albi_HPV-KER`, human.sig.results$`CTRL_vs_C.para_HPV-KER`,
                        c("HaCaT vs C.albi", "HaCaT vs C.para", 
                          "HPV-KER vs C.albi", "HPV-KER vs C.para"), "venn_diagram") 

venn.diagram(x = list(human.sig.results$CTRL_vs_C.albi_HaCat$ensembl,
                      human.sig.results$CTRL_vs_C.para_HaCat$ensembl,
                      human.sig.results$`CTRL_vs_C.albi_HPV-KER`$ensembl,
                      human.sig.results$`CTRL_vs_C.para_HPV-KER`$ensembl),
             category.names = c("HaCaT vs C.albi", "HaCaT vs C.para", 
                                "HPV-KER vs C.albi", "HPV-KER vs C.para"),
             filename = paste0(file.path(human.folder, date, plots_dir), "/","venn.png"), 
             output = T,
             
             # Output features
             imagetype="png" ,
             height = 1960 , 
             width = 1960 , 
             resolution = 300,
             compression = "lzw",
             # Circles
             lwd = 2,
             lty = 'solid',
             col = "grey",
             fill = brewer.pal(4, "Dark2"),
             
             # Numbers
             cex = .6,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 0.6,
             cat.fontface = "bold",
             cat.fontfamily = "sans")










library(clusterProfiler)
library(KEGGgraph)
library(Rgraphviz)

genelists <- lapply(human.results, "[", c("entrezID","log2FoldChange"))

for (i in 1:length(genelists)) {
  tmp <- genelists[[i]] %>% dplyr::pull(log2FoldChange,name = NULL)
  names <- genelists[[i]] %>% dplyr::pull(entrezID, name = NULL)
  
  tmp <- setNames(tmp, names)
  tmp <- sort(tmp, decreasing = T)
  tmp <- tmp[-which(is.na(names(tmp)))]
  tmp <- tmp[-which(duplicated(names(tmp)))]
  
  genelists[[i]] <- tmp
  rm(tmp)
}

genes <- lapply(human.sig_results,"[", c("entrezID","log2FoldChange"))

for (i in 1:length(genes)) {
  tmp <- genes[[i]] %>% dplyr::pull(log2FoldChange,name = NULL)
  names <- genes[[i]] %>% dplyr::pull(entrezID, name = NULL)
  
  tmp <- setNames(tmp, names)
  tmp <- sort(tmp, decreasing = T)
  tmp <- tmp[-which(is.na(names(tmp)))]
  #tmp <- tmp[-which(duplicated(names(tmp)))]
  
  
  genes[[i]] <- tmp
  rm(tmp)
}

KEGGs <- list()
for (i in 1:length(genelists)) {
  KEGGs[[length(KEGGs) + 1]] <- KEGG_plotting(names(genes[[i]]), names(genelists[[i]]))
  KEGGs[[length(KEGGs)]] <- setReadable(KEGGs[[length(KEGGs)]], 'org.Hs.eg.db', keyType = 'ENTREZID')
  
}
KEGGs_df <- lapply(KEGGs, function(x) as.data.frame(x@result))
KEGGs_df <- lapply(KEGGs_df, function(x) x[order(x$pvalue, decreasing = F),])
KEGGs_df <- lapply(KEGGs_df, function(x){
  return(
    x %>% mutate(GeneRatio = sapply(stringr::str_split(x$GeneRatio, "/"), function(y) as.numeric(y[1])/as.numeric(y[2])))
  )
})
names(KEGGs_df) <- names(genelists)

KEGG_plots <- lapply(KEGGs_df, make_FGplot)

for (i in names(KEGG_plots)) {
  ggsave(paste0("./Results/Rstudio/pictures/",i,"KEGG.png"), 
         plot = KEGG_plots[[i]], width = 12, height = 8, units = 'in')
}
sapply(names(KEGGs_df), function(x){
  openxlsx::write.xlsx(KEGGs_df[[x]], paste0("./Results/RStudio/tables/",x," KEGG.xlsx"))})





human.GOs <- list()
for (i in 1:length(genes)) {
  human.GOs[[length(human.GOs) + 1]] <- GO_plotting(names(genes[[i]]), genelists[[i]], "ALL")
}

human.GOs_df <- lapply(human.GOs, function(x) as.data.frame(x@result))
human.GOs_df <- lapply(human.GOs_df, function(x) x[order(x$pvalue, decreasing = F),])
human.GOs_df <- lapply(human.GOs_df, function(x){
  return(
    x %>% mutate(GeneRatio = sapply(stringr::str_split(x$GeneRatio, "/"), function(y) as.numeric(y[1])/as.numeric(y[2])))
  )
})
names(human.GOs_df) <- set_names

human.GOs_plots <- lapply(human.GOs_df, function(x) make_FGplot(x, "GO"))
for (i in names(human.GOs_plots)) {
  ggsave(paste0("./Results/Rstudio/pictures/",i," GO.png"), 
         plot = human.GOs_plots[[i]], width = 12, height = 12, units = 'in')
}
sapply(names(human.GOs_df), function(x){
  openxlsx::write.xlsx(human.GOs_df[[x]], paste0("./Results/RStudio/tables/",x," GO.xlsx"))})


############################################################
#GSEA
library(msigdbr)
msigdbr_df <- msigdbr(species = "Homo sapiens")
GO_sets <- msigdbr_df %>%
  dplyr::filter(gs_cat == "C5" & gs_subcat != "HPO")

human_GSEA <- list()
for (i in 1:length(genelists)) {
  
  gene_list = genelists[[i]][!duplicated(names(genelists[[i]]))]
  gene_list = sort(gene_list, decreasing = TRUE)
  
  human_GSEA[[length(human_GSEA) + 1]] <- make_gsea(gene_list, GO_sets)
  
}
names(human_GSEA) <- names(genelists)

human_GSEA.df <- lapply(human_GSEA, function(x){ 
  x %>%  
    as.data.frame(.) %>% 
    dplyr::mutate(ONTOLOGY = str_split_i(ID,"_",1),
                  ONTOLOGY = as.factor(str_remove(ONTOLOGY,"GO")),
                  Description = str_replace_all(str_remove(ID,"GOMF_|GOBP_|GOCC_"),"_"," ")) %>%
    dplyr::arrange(desc(abs(NES))) %>%
    dplyr::filter(p.adjust < 0.05) %>%
    rowwise(.) %>%
    dplyr::mutate(Count = length(unlist(strsplit(core_enrichment,"/"))))
  })

for(x in names(human_GSEA.df)){
  genes <- unlist(str_split(human_GSEA.df[[x]]$core_enrichment,"/")) 
  genes <- mapIds(org.Hs.eg.db, genes,
                 keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
  
  counts <- human_GSEA.df[[x]]$Count
  nrows <- dim(human_GSEA.df[[x]])[1]
  n <- 1
  core_enrichment <- list()
  for (i in 1:nrows) {
    last <- n + counts[i] - 1
    core_enrichment[i] <- paste(genes[n:last],collapse = "/")
    n <- last
  }
  
  core_enrichment <- unlist(core_enrichment)
  human_GSEA.df[[x]]$core_enrichment <- core_enrichment
}

human_gseaplots <- lapply(human_GSEA.df, make_gseaplot)
for (i in names(human_gseaplots)) {
  ggsave(paste0("./Results/Rstudio/pictures/",i," GSEA.png"), 
         plot = human_gseaplots[[i]], width = 18, height = 12, units = 'in')
}
sapply(names(human_GSEA.df), function(x){
  openxlsx::write.xlsx(human_GSEA.df[[x]], paste0("./Results/RStudio/tables/",x," GO (GSEA).xlsx"))})


sapply(names(human_GSEA.df), function(x){
  openxlsx::write.xlsx(human_GSEA.df[[x]], paste0("./Results/RStudio/tables/",x," GO (GSEA).xlsx"))})

 
# human_gsea_merged <- lapply(human_GSEA.df, function(x){
#   x %>% 
#     select(c("ID","ONTOLOGY","Description","Count","NES")) %>%
#     group_by(ONTOLOGY) %>%
#     dplyr::slice_head(n = 10)
# }) %>% merge.rec(., by = c("ID","Description","ONTOLOGY"), all = T, suffixes=c("",""))




#########################################################
library(ggplot2)
library(ComplexHeatmap)
library(ComplexUpset)
library(GO.db)

con <- GO.db::GO_dbconn()
go_term <- tbl(con, "go_term")
go_mf_offspring <- tbl(con, "go_mf_offspring")
go_bp_offspring <- tbl(con, "go_bp_offspring")


human.upset_base <- list(
  data.frame("geneID" = hacat.c_albi_res$sig_df$geneID,
             `CTRL vs C.albi (HaCat)` = hacat.c_albi_res$sig_df$log2FoldChange,
             row.names = row.names(hacat.c_albi_res$sig_df)),
  data.frame("geneID" = hacat.c_para_res$sig_df$geneID,
             `CTRL vs C.para (HaCat)` = hacat.c_para_res$sig_df$log2FoldChange,
             row.names = row.names(hacat.c_para_res$sig_df)),
  data.frame("geneID" = hpvker.c_albi_res$sig_df$geneID,
             `CTRL vs C.albi (HPV-KER)` = hpvker.c_albi_res$sig_df$log2FoldChange,
             row.names = row.names(hpvker.c_albi_res$sig_df)),
  data.frame("geneID" = hpvker.c_para_res$sig_df$geneID,
             `CTRL vs C.para (HPV-KER)` = hpvker.c_para_res$sig_df$log2FoldChange,
             row.names = row.names(hpvker.c_para_res$sig_df)))

human.upset_base <- merge.rec(human.upset_base, by="geneID", all=T, suffixes=c("", ""))
human.upset_base <- human.upset_base %>% 
  mutate(geneID = mapIds(org.Hs.eg.db,human.upset_base$geneID, 
                         keytype = "ENSEMBL", column = "SYMBOL")) %>%
  mutate(CTRL.vs.C.albi..HaCat. = ifelse(is.na(CTRL.vs.C.albi..HaCat.), 
                                         FALSE, TRUE)) %>%
  mutate(CTRL.vs.C.para..HaCat. = ifelse(is.na(CTRL.vs.C.para..HaCat.), 
                                         FALSE, TRUE)) %>%
  mutate(CTRL.vs.C.albi..HPV.KER. = ifelse(is.na(CTRL.vs.C.albi..HPV.KER.), 
                                           FALSE, TRUE)) %>%
  mutate(CTRL.vs.C.para..HPV.KER. = ifelse(is.na(CTRL.vs.C.para..HPV.KER.), 
                                           FALSE, TRUE)) %>%
  relocate(where(is.logical), .before = where(is.character))
human.upset_levels  <- as.factor(colnames(human.upset_base)[1:4])



###### Molecular functions
MF_gene_annot <- AnnotationDbi::select(org.Hs.eg.db,human.upset_base$geneID, "GO","SYMBOL") %>%
  dplyr::filter(ONTOLOGY =="MF" & EVIDENCE != "ND") %>%
  dplyr::select(-c(EVIDENCE,ONTOLOGY)) %>%
  distinct(SYMBOL, .keep_all = T)

MF_gene_annot <- MF_gene_annot %>% 
  dplyr::mutate(
    select(GO.db, GO, c("TERM"), "GOID")
  ) %>%
  dplyr::select(!GOID)


human.upset_MF <- merge(na.omit(human.upset_base), MF_gene_annot, 
                          by.x = "geneID", by.y="SYMBOL", all.x = T) %>%
  relocate(where(is.logical), .before = where(is.character))
human.upset_MF <- human.upset_MF %>%
  dplyr::group_by(TERM) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup()

MF_keys <- human.upset_MF %>%
  dplyr::filter(count < 10) %>%
  dplyr::pull(GO) %>%
  unique()

MF_key2id <- go_term %>%
  dplyr::filter(go_id %in% MF_keys) %>%
  dplyr::select(go_id, `_id`)

MF_id2ancestor <- dplyr::left_join(MF_key2id, go_mf_offspring, by = c(`_id` = "_offspring_id"))

MF_ancestor2go_id <- dplyr::left_join(MF_id2ancestor, go_term, by = c(`_id.y` = "_id")) %>%
  dplyr::rename(go_id = `go_id.x`, ancestor_go_id = `go_id.y`) %>%
  dplyr::select(-c(`_id`, `_id.y`, ontology, definition)) %>%
  dplyr::filter(!ancestor_go_id %in% c("GO:0005488","GO:0003674","all")) %>%
  dplyr::collect()

MF_gene_annot_full <- merge(MF_gene_annot, ancestor2go_id[,-c(4)], by.x = "GO", by.y="go_id", all = T)
colnames(MF_gene_annot_full) <- c("ID","Symbol","term","parent_ID","parent_term")

MF_gene_annot_full <- MF_gene_annot_full %>%
  dplyr::group_by(ID) %>%
  dplyr::arrange(parent_ID) %>%
  dplyr::distinct(Symbol, .keep_all = T) %>%
  ungroup()

MF_gene_annot_full <- MF_gene_annot_full %>%
  dplyr::mutate(TERM = case_when(
    ID %in% c('GO:0000976','GO:0003676','GO:0000166','GO:0004857') ~ ID,
    parent_ID %in% c('GO:000493') ~ 'GO:0004888',
    parent_ID == NA ~ ID,
    T ~ parent_ID
  )) %>%
  dplyr::select(-c(ID, term, parent_ID, parent_term)) %>%
  dplyr::rename(ID = TERM) %>% 
  dplyr::mutate(
    select(GO.db, ID, c("TERM"))
  )

human.upset_MF <- merge(human.upset_MF, MF_gene_annot_full[,c(1,2,4)], 
                          by.x = c("geneID"), by.y=c("Symbol"), all = T) %>%
  dplyr::mutate(GO = ifelse(is.na(ID), GO, ID),
                TERM = ifelse(is.na(ID), TERM.x, TERM.y)) %>%
  dplyr::select(-c(TERM.x,ID,TERM.y)) %>%
  relocate(where(is.logical), .before = where(is.character)) %>%
  relocate(where(is.numeric), .after = everything()) %>%
  dplyr::group_by(TERM) %>%
  dplyr::mutate(count = n(),
                TERM = ifelse(n() > 10, TERM, "other")) %>%
  dplyr::ungroup()

colnames(human.upset_MF)[1:4] <- names(human.results)
openxlsx::write.xlsx(human.upset_MF, "./Results/RStudio/tables/upset_MF.xlsx")



upset_MF_plot <- upset(data = human.upset_MF,
                       intersect = human.upset_levels, 
                       name='Conditions',
                       mode='inclusive_intersection',
                       min_degree = 1,
                       max_degree = 2,
                       min_size = 5,
                       set_sizes = (
                         upset_set_size(position = "right") + 
                           geom_label(aes(label=..count..), stat='count', 
                                      position = position_stack(vjust = .5))
                       ),
                       base_annotations = list(
                         'Intersection size'=(
                           intersection_size(
                             mapping=aes(
                               stat = "count",
                               fill=TERM),
                             mode='inclusive_intersection',
                             bar_number_threshold = 1,
                             counts=T,
                             color = "black")
                           + theme(plot.background=element_rect(fill='#E5D3B3', color = "transparent"),
                                   legend.background = element_blank(),
                                   legend.position = "top",
                                   legend.spacing.y = unit(5,'mm'),
                                   legend.direction = "horizontal",
                                   legend.key.size = unit(5,'mm'),
                                   legend.text = element_text(size=10))
                           + labs(y='# of genes', fill = "Molecular functions")
                         )
                       ),
                       themes=upset_modify_themes(
                         list(
                           'intersections_matrix'=theme(text=element_text(size=16)),
                           'overall_sizes'=theme(text=element_text(size=16))
                         )),
                       stripes=c('cornsilk1', 'deepskyblue1'))

ggsave("./Results/RStudio/pictures/MF_upset_plot.png",plot = upset_MF_plot,
       width = 22, height = 10, units = 'in')

####### Biological processes
BP_gene_annot <- AnnotationDbi::select(org.Hs.eg.db,human.upset_base$geneID, "GO","SYMBOL") %>%
  dplyr::filter(ONTOLOGY =="BP" & EVIDENCE != "ND") %>%
  dplyr::select(-c(EVIDENCE,ONTOLOGY)) %>%
  distinct(SYMBOL, .keep_all = T)

BP_gene_annot <- BP_gene_annot %>% 
  dplyr::mutate(
    select(GO.db, GO, c("TERM"), "GOID")
  ) %>%
  dplyr::select(!GOID)


human.upset_BP <- merge(na.omit(human.upset_base), BP_gene_annot, 
                        by.x = "geneID", by.y="SYMBOL", all.x = T) %>%
  relocate(where(is.logical), .before = where(is.character)) %>%
  dplyr::group_by(TERM) %>%
  dplyr::mutate(count = n()) %>%
  dplyr::ungroup()

BP_keys <- human.upset_BP %>%
  dplyr::filter(count < 10) %>%
  dplyr::pull(GO) %>%
  unique()

BP_key2id <- go_term %>%
  dplyr::filter(go_id %in% BP_keys) %>%
  dplyr::select(go_id, `_id`)

BP_id2ancestor <- dplyr::left_join(BP_key2id, go_bp_offspring, by = c(`_id` = "_offspring_id"))

BP_ancestor2go_id <- dplyr::left_join(BP_id2ancestor, go_term, by = c(`_id.y` = "_id")) %>%
  dplyr::rename(go_id = `go_id.x`, ancestor_go_id = `go_id.y`) %>%
  dplyr::select(-c(`_id`, `_id.y`, ontology, definition)) %>%
  dplyr::filter(!ancestor_go_id %in% c("GO:0008150","all")) %>%
  dplyr::collect()

BP_gene_annot_full <- merge(BP_gene_annot, BP_ancestor2go_id, by.x = "GO", by.y="go_id", all = T)
colnames(BP_gene_annot_full) <- c("ID","Symbol","term","parent_ID","parent_term")

BP_gene_annot_full <- BP_gene_annot_full %>%
  dplyr::group_by(ID) %>%
  dplyr::arrange(parent_ID) %>%
  dplyr::distinct(Symbol, .keep_all = T) %>%
  ungroup()

BP_gene_annot_full <- BP_gene_annot_full %>%
  dplyr::mutate(TERM = case_when(
    ID %in% c('GO:0002376','GO:0007155','GO:0001568','GO:0000165') ~ ID,
    parent_ID == NA ~ ID,
    T ~ parent_ID
  )) %>%
  dplyr::select(-c(ID, term, parent_ID, parent_term)) %>%
  dplyr::rename(ID = TERM) %>% 
  dplyr::mutate(
    select(GO.db, ID, c("TERM"))
  ) %>%
  dplyr::select(-c(GOID))

human.upset_BP <- merge(human.upset_BP, BP_gene_annot_full, 
                        by.x = c("geneID"), by.y=c("Symbol"), all = T) %>%
  dplyr::mutate(GO = ifelse(is.na(ID), GO, ID),
                TERM = ifelse(is.na(ID), TERM.x, TERM.y)) %>%
  dplyr::select(-c(TERM.x,ID,TERM.y)) %>%
  relocate(where(is.logical), .before = where(is.character)) %>%
  relocate(where(is.numeric), .after = everything()) %>%
  dplyr::group_by(TERM) %>%
  dplyr::mutate(count = n(),
                TERM = ifelse(n() > 5, TERM, "other")) %>%
  dplyr::ungroup()

colnames(human.upset_BP)[1:4] <- names(human.results)
openxlsx::write.xlsx(human.upset_BP, "./Results/RStudio/tables/upset_BP.xlsx")



upset_BP_plot <- upset(data = human.upset_BP,
                       intersect = human.upset_levels, 
                       name='Conditions',
                       mode='inclusive_intersection',
                       min_degree = 1,
                       max_degree = 2,
                       min_size = 5,
                       set_sizes = (
                         upset_set_size(position = "right") + 
                           geom_label(aes(label=..count..), stat='count', 
                                      position = position_stack(vjust = .5))
                       ),
                       base_annotations = list(
                         'Intersection size'=(
                           intersection_size(
                             mapping=aes(
                               stat = "count",
                               fill=TERM),
                             mode='inclusive_intersection',
                             bar_number_threshold = 1,
                             counts=T,
                             color = "black")
                           + theme(#plot.background=element_rect(fill='#E5D3B3', color = "transparent"),
                                   legend.background = element_blank(),
                                   legend.location = "panel",
                                   legend.position = "top",
                                   legend.title.position = "top",
                                   legend.spacing.y = unit(5,'mm'),
                                   legend.direction = "horizontal",
                                   legend.key.size = unit(5,'mm'),
                                   legend.text = element_text(size=10))
                           + labs(y='# of genes', fill = "Biological proces")
                         )
                       ),
                       themes=upset_modify_themes(
                         list(
                           'intersections_matrix'=theme(text=element_text(size=16)),
                           'overall_sizes'=theme(text=element_text(size=16))
                         )),
                       stripes=c('cornsilk1', 'deepskyblue1'))

ggsave("./Results/RStudio/pictures/BP_upset_plot.png",plot = upset_BP_plot,
       width = 22, height = 10, units = 'in')


human.gene_expr <- merge.rec(lapply(human.sig_results, select, c("geneName","log2FoldChange")),
          by = "geneName", all = T, suffixes = c("",""))
colnames(human.gene_expr)[2:5] <- set_names


human.gene_term <- merge(human.upset_MF[,-8], human.upset_BP[,-8], 
      by = c("CTRL vs C.albi (HaCat)","CTRL vs C.para (HaCat)",
             "CTRL vs C.albi (HPV-KER)","CTRL vs C.para (HPV-KER)", "geneID"),
      all = T, suffixes = c(" (BP)", " (MF)"))

human.gene_term <- merge(human.gene_expr, human.gene_term[,-c(1:4)], 
                         by.x = "geneName", by.y = "geneID",
                         all.y = T)
openxlsx::write.xlsx(human.gene_term, "./Results/RStudio/tables/human.gene_terms.xlsx")
