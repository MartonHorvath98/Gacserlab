setwd("D:/Adam/Cells")
source("../R_functions.R")


library(GenomicFeatures)
library(Rsamtools)
library(GenomeInfoDb)

human.txdb <- vulcankeTxDbFromGFF("../../grch38/Homo_sapiens.GRCh38.96.gff3", 
                              dataSource = "Ensembl", 
                              organism = "Homo sapiens")
human.genes <- exonsBy(human.txdb, by = "gene")

human.files <- list.files("./Align/", pattern = ".bam$", full.names = T)
human.list <- BamFileList(human.files, yieldSize = 2000000)

library(BiocParallel)
register(SnowParam())

library(GenomicAlignments)
human.reads <- summarizeOverlaps(features = human.genes, 
                                 reads = human.list, 
                                 mode = "Union", 
                                 ignore.strand = F, 
                                 singleEnd = F, 
                                 fragments = T, 
                                 preprocess.reads = invertStrand)
head(assay(human.reads))
readcounts <- assay(human.reads)


library(stringr)
bam.names <- as.list(lapply(strsplit(human.files, "/"),"[[",3)) #get samplenames
bam.treatment <- sapply(sapply(bam.names,strsplit, "_"), "[[", 2)
bam.treatment <- factor(bam.treatment, levels = c("C","Ca","Cp"), 
                        labels = c("ctrl","C.albi","C.para"))
bam.cell_line <- sapply(sapply(bam.names,strsplit, "_"), "[[", 1)
bam.cell_line <- factor(bam.cell_line, levels = c("HC","HK"), 
                        labels = c("HaCat","HPVKER"))
bam.condition <-  as.factor(paste(bam.cell_line, bam.treatment))

#separate treatment
bam.run <- as.factor(sapply(bam.names, function(x) paste("run", (substr(x, str_locate(x, ".bam") - 1, str_locate(x, ".bam") - 1))))) #separate run

human.coldata <- data.frame(samplenames = as.character(bam.names),
                      cells = bam.cell_line,
                      treatment = bam.treatment,
                      condition = bam.condition,
                      run = bam.run) #vulcanke metadate table

colnames(readcounts) <- as.character(human.coldata$condition)

write.table(readcounts, "./Results/RStudio/tables/total_readcounts.csv", 
            sep = ",", na = "NA", dec = ".", row.names = T, col.names = T)


library(DESeq2)
library(edgeR)
human.deseq <- make_deseq(matrix = assay(human.reads),
                          coldata = human.coldata,
                          design = "cells + treatment")
resultsNames(human.deseq$dds)

hacat.deseq <- make_deseq(matrix = assay(human.reads)[,1:9],
                          coldata = human.coldata[1:9,],
                          design = "treatment")
resultsNames(hacat.deseq$dds)

hpvker.deseq <- make_deseq(matrix = assay(human.reads)[,10:18],
                          coldata = human.coldata[10:18,],
                          design = "treatment")
resultsNames(hpvker.deseq$dds)

library(ggplot2)
library(ggbiplot)
library(ggrepel)
hacat.pca <- make_pca(rld = hacat.deseq$rld, 
                      group = hacat.deseq$rld$treatment,
                      labs = c("Control", "C.albicans (SC5314)",
                               "C.parapsilosis (CLIB214)"),
                      cols = c("C.para" = "salmon",
                               "C.albi" = "darkred",
                               "ctrl" = "turquoise"))
hacat.pca <- hacat.pca + 
  theme(legend.position = "bottom")
ggsave("./Results/Rstudio/pictures/hacat_pca.png", plot = hacat.pca,
       width = 10, height = 10, units = 'in')

hpvker.pca <- make_pca(rld = hpvker.deseq$rld, 
                      group = hpvker.deseq$rld$treatment,
                      labs = c("Control", "C.albicans (SC5314)",
                               "C.parapsilosis (CLIB214)"),
                      cols = c("C.para" = "salmon",
                               "C.albi" = "darkred",
                               "ctrl" = "turquoise"))
hpvker.pca <- hpvker.pca + 
  theme(legend.position = "bottom")
ggsave("./Results/Rstudio/pictures/hpvker_pca.png", plot = hpvker.pca,
       width = 10, height = 10, units = 'in')

cells.pca <- make_pca(rld = human.deseq$rld, 
                      group = human.deseq$rld$condition,
                      labs = c("C.albicans (HaCat)", "C.parapsilosis (HaCat)", "Control (HaCat)",
                               "C.albicans (HPV-KER)", "C.parapsilosis (HPV-KER)", "Control (HPV-KER)"),
                      cols = c("HaCat C.para" = "salmon",
                               "HaCat C.albi" = "darkred",
                               "HaCat ctrl" = "orange",
                               "HPVKER C.para" = "blue",
                               "HPVKER C.albi" = "purple",
                               "HPVKER ctrl" = "turquoise"))
ggsave("./Results/Rstudio/pictures/total_pca.png", plot = cells.pca,
       width = 12, height = 12, units = 'in')

library(fdrtool)
library(openxlsx)
library(AnnotationDbi)
library(org.Hs.eg.db)
# cells.res <- get_results(human.deseq$dds, 1.5, 0.05, 
#                          contrast = list("cells_HPVKER_vs_HaCat"), 
#                          name = c(''))

hacat.c_albi_res <- get_results(hacat.deseq$dds, 1.5, 0.05, 
                                contrast = list("treatment_C.albi_vs_ctrl"), 
                                name = c(''))

hacat.c_para_res <- get_results(hacat.deseq$dds, 1.5, 0.05, 
                                contrast = list("treatment_C.para_vs_ctrl"), 
                                name = c(''))

hpvker.c_albi_res <- get_results(hpvker.deseq$dds, 1.5, 0.05, 
                                 contrast = list("treatment_C.albi_vs_ctrl"), 
                                 name = c(''))

hpvker.c_para_res <- get_results(hpvker.deseq$dds, 1.5, 0.05, 
                                 contrast = list("treatment_C.para_vs_ctrl"), 
                                 name = c(''))

human.sig_results <- list(
  "CTRL vs C.albi (HaCat)" = hacat.c_albi_res$sig_df,
  "CTRL vs C.para (HaCat)" = hacat.c_para_res$sig_df,
  "CTRL vs C.albi (HPV-KER)" = hpvker.c_albi_res$sig_df,
  "CTRL vs C.para (HPV-KER)" = hpvker.c_para_res$sig_df
)
set_names <- names(human.sig_results)

human.sig_results <- lapply(names(human.sig_results), function(x){
  return(human.sig_results[[x]] %>%
           dplyr::mutate(
             esnembl = geneID,
             geneName = mapIds(org.Hs.eg.db, human.sig_results[[x]]$geneID,
                               keytype = "ENSEMBL", column = "SYMBOL",
                               multiVals = "first"),
             entrezID = mapIds(org.Hs.eg.db, human.sig_results[[x]]$geneID,
                               keytype = "ENSEMBL", column = "ENTREZID",
                               multiVals = "first")) %>%
           dplyr::select(!geneID) %>%
           dplyr::relocate(where(is.character), .before = where(is.numeric)) %>%
           dplyr::relocate(significance, .after = everything()))
  
})
names(human.sig_results) <- set_names
sapply(names(human.sig_results), function(x){
  openxlsx::write.xlsx(human.sig_results[[x]], paste0("./Results/RStudio/tables/",x,".xlsx"))})


library(genefilter)
library(pheatmap)
library(RColorBrewer)
library(devtools)
library(ComplexHeatmap)
library(dendextend)
library(dendsort)
library(cluster)
library(factoextra)
library(circlize)
human.mat <- rlog(assay(human.reads))
human.count.mat <- make_matrix(human.mat, 
                               unique(c(#rownames(cells.res$sig_df),
                                        rownames(hacat.c_albi_res$sig_df),
                                        rownames(hacat.c_para_res$sig_df),
                                        rownames(hpvker.c_albi_res$sig_df),
                                        rownames(hpvker.c_para_res$sig_df))))

fviz_nbclust(human.count.mat, kmeans, method = "wss")
human.annot <- HeatmapAnnotation(Cells = as.factor(human.coldata$cells),
                           Treament = as.factor(human.coldata$treatment),
                           show_annotation_name = F)                                        

human.h1 <- expression_heatmap(mat = human.count.mat, clusters = 4, columns = 2,
                              coldata = human.coldata, annot = human.annot,
                              method1 = "euclidean", method2 = "ward.D")

library(dplyr)
library(plyr)
library(tidyr)
library(tibble)
human.expr <- list(
  # data.frame("geneID" = cells.res$sig_df$geneID,
  #            `HaCat vs HPV-KER` = cells.res$sig_df$log2FoldChange, 
  #            row.names = row.names(cells.res$sig_df)),
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

human.expr.df <- merge.rec(human.expr, by="geneID", all=T, suffixes=c("", ""))
human.expr.df <- human.expr.df %>% column_to_rownames("geneID") %>%
  setNames(., c("CTRL vs C.albi (HaCat)", "CTRL vs C.para (HaCat)",#"HaCat vs HPV-KER", 
                "CTRL vs C.albi (HPV-KER)","CTRL vs C.para (HPV-KER)"))


human.expr.mat <- logFC_heatmap(mat = as.matrix(human.expr.df), dend = human.h1$dendrogram)

png("./Results/Rstudio/pictures/total_heatmap.png", 
    width = 14, height = 10, units = 'in', res = 300)
draw(human.h1$heatmap + human.expr.mat, merge_legend = T)
dev.off()

library(scales)
library(viridis)
human.results <- list(
  #"HaCat vs HPV-KER" = cells.res$df,
  "CTRL vs C.albi (HaCat)" = hacat.c_albi_res$df,
  "CTRL vs C.para (HaCat)" = hacat.c_para_res$df,
  "CTRL vs C.albi (HPV-KER)" = hpvker.c_albi_res$df,
  "CTRL vs C.para (HPV-KER)" = hpvker.c_para_res$df
)
human.MA <- list()
for (i in names(human.results)) {
  human.MA[[i]] <- MA_plotting(human.results[[i]])
  ggsave(paste("./Results/Rstudio/pictures/",i,"MAplot.png"), 
                plot = human.MA[[i]], width = 10, height = 6, units = 'in')
}

library(org.Hs.eg.db)
library(ggplot2)
#colourPalette <- 
human.results <- lapply(human.results, function(x){
  return(x %>%
           dplyr::rename(., ensembl = geneID) %>%
           dplyr::mutate(geneID = mapIds(org.Hs.eg.db, x$ensembl, 
                                         keytype = "ENSEMBL", column = "SYMBOL",multiVals = "first"),
                         entrezID = mapIds(org.Hs.eg.db, x$ensembl, 
                                           keytype = "ENSEMBL", column = "ENTREZID",multiVals = "first")) %>%
           dplyr::relocate(where(is.character), .before = where(is.numeric)) %>%
           dplyr::relocate(significance, .after = everything()))
    
})


human.expr_plot <- list()
for (i in names(human.results)) {
  human.expr_plot[[i]] <- vulcan_plotting(human.results[[i]])
  ggsave(paste("./Results/Rstudio/pictures/",i,"vulcanplot.png"), 
         plot = human.expr_plot[[i]], width = 12, height = 8, units = 'in')
}

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
                           + theme(plot.background=element_rect(fill='#E5D3B3', color = "transparent"),
                                   legend.background = element_blank(),
                                   legend.position = "top",
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
