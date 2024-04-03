############################################## 
# 1.) Load analysis data                     #
##############################################
mir_data_dir <- "mir_data"
if (!dir.exists(mir_data_dir)) {
  dir.create(mir_data_dir)
}

#plots directory
mir_plots_dir <- "mir_plots"
if (!dir.exists(mir_plots_dir)) {
  dir.create(mir_plots_dir)
}

#results directory
mir_results_dir <- "mir_results" 
if (!dir.exists(mir_results_dir)) {
  dir.create(mir_results_dir)
}

miRcount.paths <- list.files("E:/HSC2_Reni/miRNA/mirdeep2/counts/", 
                          pattern = ".xlsx$", full.names = T)

miRcount.files <- lapply(miRcount.paths, read.xlsx, colNames = TRUE)

library(dplyr)
miRcount.list <- lapply(miRcount.files, function(x) {
  x  %>% 
    select(c(1:2)) %>%
    setNames(c("miRname","Counts")) %>%
    dplyr::group_by(miRname) %>%
    dplyr::summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
})
miRcount.table <- merge.rec(miRcount.list, by = "miRname", suffixes = c("",""))

colnames(miRcount.table)[-1] <- file.names
write.table(miRcount.table, paste(mir_data_dir,"miR_readcounts.csv", sep = "/"), 
            sep =",", na = "NA", dec = ".", row.names = F, col.names = T)

##############################################
# 2.) Make differential analyses with DESeq2 #
##############################################
miRreadcounts <- read.csv(file = paste(mir_data_dir,"miR_readcounts.csv", sep = "/"),
                       sep = ",", header = T, na.strings = NA, row.names = 1)
miRreadcounts <- as.matrix(miRreadcounts)
miRreadcounts <- list("1h" = miRreadcounts[,1:6], "6h" = miRreadcounts[,7:18])

miR_DESeq.1h <- calc_DiffExp(matrix = miRreadcounts[["1h"]], 
                         coldata = coldata.list[["1h"]],
                         design = "condition")
miR_DESeq.6h <- calc_DiffExp(matrix = miRreadcounts[["6h"]],
                         coldata = coldata.list[["6h"]],
                         design = "condition")

(miR_pca.1h <- make_pca(miR_DESeq.1h$dds_norm, "samples", coldata.list[["1h"]]$condition))
ggsave(paste(plots_dir,"miR_1h_pca.png",sep = "/"), plot = miR_pca.1h,
       width = 8, height = 8, units = 'in')

(miR_pca.6h <- make_pca(miR_DESeq.6h$dds_norm, "samples", coldata.list[["6h"]]$condition))
ggsave(paste(plots_dir,"miR_6h_pca.png",sep = "/"), plot = miR_pca.6h,
       width = 8, height = 8, units = 'in')

miR_Ca11_1h <- miR_results(miR_DESeq.1h$dds,
                           contrast = list("condition_CA11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1,
                           fdr_treshold = 0.1)
sapply(names(miR_Ca11_1h), function(x){
  openxlsx::write.xlsx(miR_Ca11_1h[[x]], paste0(mir_results_dir, "/Ca11_1h_miR_",x,".xlsx"), rowNames = T)
})

miR_Ca11_6h <- miR_results(miR_DESeq.6h$dds,
                            contrast = list("condition_CA11_vs_ctrl"),
                            name = c(''),
                            lfc_treshold = 1,
                            fdr_treshold = 0.1)
sapply(names(miR_Ca11_6h), function(x){
  openxlsx::write.xlsx(miR_Ca11_6h[[x]], paste0(mir_results_dir, "/Ca11_6h_miR_",x,".xlsx"), rowNames = T)
})

miR_Cp11_6h <- miR_results(miR_DESeq.6h$dds,
                           contrast = list("condition_CP11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1,
                           fdr_treshold = 0.1)
sapply(names(miR_Cp11_6h), function(x){
  openxlsx::write.xlsx(miR_Cp11_6h[[x]], paste0(mir_results_dir, "/Cp11_6h_miR_",x,".xlsx"), rowNames = T)
})

miR_Cp51_6h <- miR_results(miR_DESeq.6h$dds,
                           contrast = list("condition_CP51_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1,
                           fdr_treshold = 0.1)
sapply(names(miR_Cp51_6h), function(x){
  openxlsx::write.xlsx(miR_Cp51_6h[[x]], paste0(mir_results_dir, "/Cp51_6h_miR_",x,".xlsx"), rowNames = T)
})


(miR_Ca11_1h_plot <- miRNA_plot(miR_Ca11_1h$sig_df))
ggsave(paste(mir_plots_dir,"miR_Ca11_1h_plot.png",sep = "/"), plot = miR_Ca11_1h_plot,
       width = 10, height = 8, units = 'in')

(miR_Ca11_6h_plot <- miRNA_plot(miR_Ca11_6h$sig_df))
ggsave(paste(mir_plots_dir,"miR_Ca11_6h_plot.png",sep = "/"), plot = miR_Ca11_6h_plot,
       width = 10, height = 8, units = 'in')

(miR_Cp11_6h_plot <- miRNA_plot(miR_Cp11_6h$sig_df))
ggsave(paste(mir_plots_dir,"miR_Cp11_6h_plot.png",sep = "/"), plot = miR_Cp11_6h_plot,
       width = 10, height = 8, units = 'in')

(miR_Cp51_6h_plot <- miRNA_plot(miR_Cp51_6h$sig_df))
ggsave(paste(mir_plots_dir,"miR_Cp51_6h_plot.png",sep = "/"), plot = miR_Cp51_6h_plot,
       width = 10, height = 3, units = 'in')

##############################################
# 3.) Check miRNA target genes (predicted)   #
##############################################

# C.albicans (1h)
Ca11_1h_miR_unfiltered.file <- require_file(paste(mir_data_dir, "Ca11_1h_miRWalk_unfiltered.csv", sep = "/"))
Ca11_1h_miR_unfiltered.df <- Ca11_1h_miR_unfiltered.file %>%
  dplyr::group_by(mirnaid) %>%
  dplyr::distinct(genesymbol,.keep_all = T) %>%
  dplyr::select(mirnaid:genesymbol) %>%
  dplyr::ungroup(.)

Ca11_1h_miR_unfiltered.df <- merge(miR_Ca11_1h$sig_df[,c(7,2,5)], Ca11_1h_miR_unfiltered.df,
                                    by.x = "miRname", by.y = "mirnaid", all = T)
Ca11_1h_miR_unfiltered.filter <- merge(Ca11_1h.list$sig_df[,c(8,2,5)], Ca11_1h_miR_unfiltered.df,
                                    by.x = "geneID", by.y = "genesymbol", suffixes = c("",".miR"))
Ca11_1h_miR_unfiltered.filter <- Ca11_1h_miR_unfiltered.filter %>%
  dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                  log2FoldChange < 0 & log2FoldChange.miR > 0)

Ca11_1h_miR_unfiltered.genes <- Ca11_1h_miR_unfiltered.filter %>%
  dplyr::mutate(pairs = as.character(stringr::str_glue("{geneID} [{miRname}]"))) %>%
  dplyr::pull(pairs, name = "geneID")
# visualisation
(Ca11_1h_miR_unfiltered.plot <- get_CategoryExpressionPlot(results = Ca11_1h.list$df, 
                                                        sig_log2FC = 1.5, 
                                                        sig_pval = 0.05, 
                                                        targets = Ca11_1h_miR_unfiltered.genes))
ggsave(paste(mir_plots_dir,"Ca11_1h_miR_predicted_targets.png",sep = "/"), plot = Ca11_1h_miR_unfiltered.plot,
       width = 12, height = 10, units = 'in')

# C.albicans (6h)
Ca11_6h_miR_unfiltered.file <- require_file(paste(mir_data_dir, "Ca11_6h_miRWalk_unfiltered.csv", sep = "/"))
Ca11_6h_miR_unfiltered.df <- Ca11_6h_miR_unfiltered.file %>%
  dplyr::group_by(mirnaid) %>%
  dplyr::distinct(genesymbol,.keep_all = T) %>%
  dplyr::select(mirnaid:genesymbol) %>%
  dplyr::ungroup(.)

Ca11_6h_miR_unfiltered.df <- merge(miR_Ca11_6h$sig_df[,c(7,2,5)], Ca11_6h_miR_unfiltered.df,
                                by.x = "miRname", by.y = "mirnaid", all = T)
Ca11_6h_miR_unfiltered.filter <- merge(Ca11_6h.list$sig_df[,c(8,2,5)], Ca11_6h_miR_unfiltered.df,
                                    by.x = "geneID", by.y = "genesymbol", suffixes = c("",".miR"))

Ca11_6h_miR_unfiltered.filter <- Ca11_6h_miR_unfiltered.filter %>%
  dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                  log2FoldChange < 0 & log2FoldChange.miR > 0)

Ca11_6h_miR_unfiltered.genes <- Ca11_6h_miR_unfiltered.filter %>%
  dplyr::select(c('geneID','miRname')) %>%
  
  dplyr::mutate(pairs = as.character(stringr::str_glue("{geneID} [{miRname}]"))) %>%
  dplyr::pull(pairs, name = "geneID")
# visualisation
(Ca11_6h_miR_unfiltered.plot <- get_CategoryExpressionPlot(results = Ca11_6h.list$df, 
                                                           sig_log2FC = 1.5, 
                                                           sig_pval = 0.05, 
                                                           targets = Ca11_6h_miR_unfiltered.genes))
ggsave(paste(mir_plots_dir,"Ca11_1h_miR_predicted_targets.png",sep = "/"), plot = Ca11_1h_miR_unfiltered.plot,
       width = 12, height = 10, units = 'in')


# C.parapsilosis MOI 1:1 (6h)
Cp11_6h_miR_unfiltered.file <- require_file(paste(mir_data_dir, "Cp11_6h_miRWalk_unfiltered.csv", sep = "/"))
Cp11_6h_miR_unfiltered.df <- Cp11_6h_miR_unfiltered.file %>%
  dplyr::group_by(mirnaid) %>%
  dplyr::distinct(genesymbol,.keep_all = T) %>%
  dplyr::select(mirnaid:genesymbol) %>%
  dplyr::ungroup(.)

Cp11_6h_miR_unfiltered.df <- merge(miR_Cp11_6h$sig_df[,c(7,2,5)], Cp11_6h_miR_unfiltered.df,
                                by.x = "miRname", by.y = "mirnaid", all = T)
Cp11_6h_miR_unfiltered.filter <- merge(Cp11_6h.list$sig_df[,c(8,2,5)], Cp11_6h_miR_unfiltered.df,
                                    by.x = "geneID", by.y = "genesymbol", suffixes = c("",".miR"))
Cp11_6h_miR_unfiltered.filter <- Cp11_6h_miR_unfiltered.filter %>%
  dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                log2FoldChange < 0 & log2FoldChange.miR > 0)

Cp11_6h_miR_unfiltered.genes <- Cp11_6h_miR_unfiltered.filter %>%
  dplyr::mutate(pairs = as.character(stringr::str_glue("{geneID} [{miRname}]"))) %>%
  dplyr::pull(pairs, name = "geneID")
# visualisation
(Cp11_6h_miR_unfiltered.plot <- get_CategoryExpressionPlot(results = Cp11_6h.list$df, 
                                                           sig_log2FC = 1.5, 
                                                           sig_pval = 0.05, 
                                                           targets = Cp11_6h_miR_unfiltered.genes))
ggsave(paste(mir_plots_dir,"Cp11_6h_miR_predicted_targets.png",sep = "/"), plot = Cp11_6h_miR_unfiltered.plot,
       width = 12, height = 10, units = 'in')


# C.parapsilosis MOI 5:1 (6h)
Cp51_6h_miR_unfiltered.file <- require_file(paste(mir_data_dir, "Cp51_6h_miRWalk_unfiltered.csv", sep = "/"))
Cp51_6h_miR_unfiltered.df <- Cp51_6h_miR_unfiltered.file %>%
  dplyr::group_by(mirnaid) %>%
  dplyr::distinct(genesymbol,.keep_all = T) %>%
  dplyr::select(mirnaid:genesymbol) %>%
  dplyr::ungroup(.)

Cp51_6h_miR_unfiltered.df <- merge(miR_Cp51_6h$sig_df[,c(7,2,5)], Cp51_6h_miR_unfiltered.df,
                                by.x = "miRname", by.y = "mirnaid", all = T)
Cp51_6h_miR_unfiltered.filter <- merge(Cp51_6h.list$sig_df[,c(8,2,5)], Cp51_6h_miR_unfiltered.df,
                                    by.x = "geneID", by.y = "genesymbol", suffixes = c("",".miR"))
Cp51_6h_miR_unfiltered.filter <- Cp51_6h_miR_unfiltered.filter %>%
 dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
               log2FoldChange < 0 & log2FoldChange.miR > 0)

Cp51_6h_miR_unfiltered.genes <- Cp51_6h_miR_unfiltered.filter %>%
  dplyr::mutate(pairs = as.character(stringr::str_glue("{geneID} [{miRname}]"))) %>%
  dplyr::pull(pairs, name = "geneID")
# visualisation
(Cp51_6h_miR_unfiltered.plot <- get_CategoryExpressionPlot(results = Cp51_6h.list$df, 
                                                           sig_log2FC = 1.5, 
                                                           sig_pval = 0.05, 
                                                           targets = Cp51_6h_miR_unfiltered.genes))
ggsave(paste(mir_plots_dir,"Cp51_6h_miR_predicted_targets.png",sep = "/"), plot = Cp51_6h_miR_unfiltered.plot,
       width = 12, height = 10, units = 'in')


##############################################
# 4.) Check miRNA target genes (validated)   #
##############################################

# C.albicans (1h)
Ca11_1h_miR_targets.file <- require_file(paste(mir_data_dir, "Ca11_1h_miRWalk_targets.csv", sep = "/"))
Ca11_1h_miR_targets.df <- Ca11_1h_miR_targets %>%
  dplyr::group_by(mirnaid) %>%
  dplyr::distinct(genesymbol,.keep_all = T) %>%
  dplyr::select(mirnaid:genesymbol) %>%
  dplyr::ungroup(.)

Ca11_1h_miR_targets.df <- merge(miR_Ca11_1h$sig_df[,c(7,2,5)], Ca11_1h_miR_targets.df,
                                by.x = "miRname", by.y = "mirnaid", all = T)
Ca11_1h_miR_targets.filter <- merge(Ca11_1h.list$sig_df[,c(8,2,5)], Ca11_1h_miR_targets.df,
                                    by.x = "geneID", by.y = "genesymbol", suffixes = c("",".miR"))

# C.albicans (6h)
Ca11_6h_miR_targets.file <- require_file(paste(mir_data_dir, "Ca11_6h_miRWalk_targets.csv", sep = "/"))
Ca11_6h_miR_targets.df <- Ca11_6h_miR_targets.file %>%
  dplyr::group_by(mirnaid) %>%
  dplyr::distinct(genesymbol,.keep_all = T) %>%
  dplyr::select(mirnaid:genesymbol) %>%
  dplyr::ungroup(.)

Ca11_6h_miR_targets.df <- merge(miR_Ca11_6h$sig_df[,c(7,2,5)], Ca11_6h_miR_targets.df,
                                by.x = "miRname", by.y = "mirnaid", all = T)
Ca11_6h_miR_targets.filter <- merge(Ca11_6h.list$sig_df[,c(8,2,5)], Ca11_6h_miR_targets.df,
                                    by.x = "geneID", by.y = "genesymbol", suffixes = c("",".miR"))

Ca11_6h_miR_targets.filter <- Ca11_6h_miR_targets.filter %>%
  dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                  log2FoldChange < 0 & log2FoldChange.miR > 0)

Ca11_6h_miR_targets.genes <- Ca11_6h_miR_targets.filter %>%
  dplyr::mutate(pairs = as.character(stringr::str_glue("{geneID} [{miRname}]"))) %>%
  dplyr::pull(pairs, name = "geneID")

# visualisation
(Ca11_6h_miR_targets.plot <- get_CategoryExpressionPlot(results = Ca11_6h.list$df, 
                                                       sig_log2FC = 1.5, 
                                                       sig_pval = 0.05, 
                                                       targets = Ca11_6h_miR_targets.genes))
ggsave(paste(mir_plots_dir,"Ca11_6h_miR_validated_targets.png",sep = "/"), plot = Ca11_6h_miR_targets.plot,
       width = 12, height = 10, units = 'in')
# C.parapsilosis MOI 1:1 (6h)
Cp11_6h_miR_targets.file <- require_file(paste(mir_data_dir, "Cp11_6h_miRWalk_targets.csv", sep = "/"))
Cp11_6h_miR_targets.df <- Cp11_6h_miR_targets.file %>%
  dplyr::group_by(mirnaid) %>%
  dplyr::distinct(genesymbol,.keep_all = T) %>%
  dplyr::select(mirnaid:genesymbol) %>%
  dplyr::ungroup(.)

Cp11_6h_miR_targets.df <- merge(miR_Cp11_6h$sig_df[,c(7,2,5)], Cp11_6h_miR_targets.df,
                                by.x = "miRname", by.y = "mirnaid", all = T)
Cp11_6h_miR_targets.filter <- merge(Cp11_6h.list$sig_df[,c(8,2,5)], Cp11_6h_miR_targets.df,
                                    by.x = "geneID", by.y = "genesymbol", suffixes = c("",".miR"))
# Cp11_6h_miR_targets.filter <- Cp11_6h_miR_targets.filter %>%
#   dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
#                   log2FoldChange < 0 & log2FoldChange.miR > 0)

# C.parapsilosis MOI 5:1 (6h)
Cp51_6h_miR_targets.file <- require_file(paste(mir_data_dir, "Cp51_6h_miRWalk_targets.csv", sep = "/"))
Cp51_6h_miR_targets.df <- Cp51_6h_miR_targets.file %>%
  dplyr::group_by(mirnaid) %>%
  dplyr::distinct(genesymbol,.keep_all = T) %>%
  dplyr::select(mirnaid:genesymbol) %>%
  dplyr::ungroup(.)

Cp51_6h_miR_targets.df <- merge(miR_Cp51_6h$sig_df[,c(7,2,5)], Cp51_6h_miR_targets.df,
                                by.x = "miRname", by.y = "mirnaid", all = T)
Cp51_6h_miR_targets.filter <- merge(Cp51_6h.list$sig_df[,c(8,2,5)], Cp51_6h_miR_targets.df,
                                    by.x = "geneID", by.y = "genesymbol", suffixes = c("",".miR"))
# Cp51_6h_miR_targets.filter <- Cp51_6h_miR_targets.filter %>%
#   dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
#                 log2FoldChange < 0 & log2FoldChange.miR > 0)

##############################################
# 6.) GSEA with (validated) miRNA targets    #
##############################################

# C.albicans MOI 1:1 (6h)
Ca11_6h_mir_GSEA.genelist <- subset(Ca11_6h.list$df, geneID %in% Ca11_6h_miR_targets.df$genesymbol) %>%
  na.omit(.) %>%
  dplyr::pull(log2FoldChange, name = entrezID)
#KEGG
Ca11_6h_mir_GSEA.kegg <- kegga(names(Ca11_6h_mir_GSEA.genelist))
Ca11_6h_mir_GSEA.kegg <- topKEGG(Ca11_6h_mir_GSEA.kegg)
#GO
Ca11_6h_mir_GSEA.go <- goana(names(Ca11_6h_mir_GSEA.genelist))
Ca11_6h_mir_GSEA.go <- Ca11_6h_mir_GSEA.go  %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column('ID') %>%
  dplyr::filter(P.DE < 0.05 & 15 < N & N < 500) %>%
  dplyr::group_split(Ont) 

names(Ca11_6h_mir_GSEA.go) <- c('BP','CC','MF')

Ca11_6h_miR_simm <- list()
sapply(c('BP','CC','MF'), function(x){
  Ca11_6h_miR_simm[[x]] <- make_simMatrix(Ca11_6h_mir_GSEA.go[[x]], as.character(x), treshold = .9)
})
