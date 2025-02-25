############################################## 
# 1.) Load analysis data                     #
##############################################
####################################################
# 1.) Set up working directory and directory tree  #
##################################################### Set downstream path
miR.folder <- "miRNA"
# Create the sub folders for: results, data, and pictures
if (!dir.exists(file.path(data_dir,miR.folder))) {
  dir.create(file.path(data_dir,miR.folder)) # create the data folder
}

# Create the dated results folder
# get the current date
date <- format(Sys.Date(), "%Y-%m-%d")
if (!dir.exists(file.path(results_dir, miR.folder, date))) {
  dir.create(file.path(results_dir, miR.folder, date)) # create the dated results folder
  dir.create(file.path(results_dir, miR.folder, date, tables_dir)) # create the tables folder
  dir.create(file.path(results_dir, miR.folder, date, plots_dir)) # create the plots folder
}

#######################################
# 2.) Load functions for the analyses #
#######################################
source("scripts/functions.R")
source("scripts/packages.R")

########################## 
# 3.) Load analysis data #
##########################
# Copy readcounts from their source folder to the data folder
if (exists("miR.readcounts") == F) {
  miR.reads <- lapply(list.files(file.path(data_dir, miR.folder), 
                                  pattern = ".xlsx$",
                                  full.names = T), 
                       read.xlsx, colNames = TRUE)
  
  miR.counts <- lapply(miR.reads, function(x) {
    x  %>% 
      select(c(1:2)) %>%
      setNames(c("miRname","Counts")) %>%
      dplyr::group_by(miRname) %>%
      dplyr::summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
  })
  miR.counts <- merge.rec(miR.counts, by = "miRname", suffixes = c("",""))
  
  miR.names <- list.files(file.path(data_dir, miR.folder),
                           pattern = ".xlsx$", full.names = F) %>% 
    gsub(".xlsx", "", .) %>% 
    gsub("miRNAs_expressed_", "", .)
  
  colnames(miR.counts)[-1] <- miR.names
  write.table(miR.counts, file.path(data_dir,miR.folder, "miR_readcounts.csv"),
              sep =",", na = "NA", dec = ".", row.names = F, col.names = T)
  
  rm(ls = ls()[grep("miR.", ls())])
}

##############################################
# 2.) Make differential analyses with DESeq2 #
##############################################
# Tissue specific mean expression of miRNAs
miR.tissue_expression <- read_parquet(file.path(data_dir,miR.folder,"miRNA_tissue_data.parquet"))

# Load the count data
miR.tissue_expression <- miR.tissue_expression %>%
  dplyr::filter(Biotype == "tissue" & Metric == "mean") %>% 
  dplyr::filter(Tissue %in% c("esophagus", "lymph_node", 
                              "thyroid", "tongue")) %>%
  dplyr::select(-c(1,2,4)) %>%
  tidyr::pivot_longer(cols = !Tissue, names_to = "miRname", values_to = "Expression.Mean") %>%
  dplyr::mutate(Tissue = factor(Tissue, 
                                levels = c("esophagus", "lymph_node", 
                                           "thyroid", "tongue")),
                miRname = factor(miRname)) %>% 
  tidyr::pivot_wider(names_from = Tissue, names_glue = "{Tissue} (mean expr.)",
                     values_from = Expression.Mean)


# Load the count data
miR.readcounts <- read.csv(file = file.path(data_dir, miR.folder, "miR_readcounts.csv"),
                       sep = ",", header = T, na.strings = NA, row.names = 1)
miR.readcounts <- as.matrix(miR.readcounts)
# Split the data at the different times of infection
miR.readcounts <- list("1h" = miR.readcounts[,1:6], "6h" = miR.readcounts[,7:18])

# Run DESeq2 on the 1 hour infection model data
miR.DESeq.1h <- calc_DiffExp(matrix = miR.readcounts[["1h"]], 
                             coldata = HSC2.coldata[["1h"]],
                             design = "condition")
# Run DESeq2 on the 6 hour infection model data
miR.DESeq.6h <- calc_DiffExp(matrix = miR.readcounts[["6h"]],
                         coldata = HSC2.coldata[["6h"]],
                         design = "condition")

(miR.pca <- list(
  # PCA analysis of the 1 hour infection model
  "1h" = make_pca(miR.DESeq.1h$dds_norm, "samples", HSC2.coldata[["1h"]]$condition),
  # PCA analysis of the 6 hour infection model
  "6h" = make_pca(miR.DESeq.6h$dds_norm, "samples", HSC2.coldata[["6h"]]$condition)))
# Save the PCA plots
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_1h_pca.png"),
       plot = miR.pca$`1h`, width = 8, height = 8, units = 'in')
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_6h_pca.png"),
       plot = miR.pca$`6h`, width = 8, height = 8, units = 'in')

# Extract significant results for the 1 hour infection model
miR.Ca11_1h.res <- miR_results(miR.DESeq.1h$dds,
                           contrast = list("condition_CA11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1,
                           p_treshold = 0.05,
                           tissue = miR.tissue_expression)
# Check the number of DEGs:
table(miR.Ca11_1h.res$sig_df$significance) # 0 down-regulated, 3 up-regulated

# save to excel
openxlsx::write.xlsx(miR.Ca11_1h.res[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Ca11_1h_miR_expr.xlsx"))

# Extract significant results for the 6 hour infection model
## C.albicans MOI 1:1 (6h)
miR.Ca11_6h.res <- miR_results(miR.DESeq.6h$dds,
                            contrast = list("condition_CA11_vs_ctrl"),
                            name = c(''),
                            lfc_treshold = 1,
                            p_treshold = 0.05,
                            tissue = miR.tissue_expression)
# Check the number of DEGs:
table(miR.Ca11_6h.res$sig_df$significance) # 15 down-regulated, 8 up-regulated

# Save to excel
openxlsx::write.xlsx(miR.Ca11_6h.res[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Ca11_6h_miR_expr.xlsx"))

## C.parapsilosis MOI 1:1 (6h)
miR.Cp11_6h.res <- miR_results(miR.DESeq.6h$dds,
                           contrast = list("condition_CP11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1,
                           p_treshold = 0.05,
                           tissue = miR.tissue_expression)
# Check the number of DEGs:
table(miR.Cp11_6h.res$sig_df$significance) # 4 down-regulated, 2 up-regulated

# Save to excel
openxlsx::write.xlsx(miR.Cp11_6h.res[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Cp11_6h_miR_expr.xlsx"))

## C.parapsilosis MOI 5:1 (6h)
miR.Cp51_6h.res <- miR_results(miR.DESeq.6h$dds,
                           contrast = list("condition_CP51_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1,
                           p_treshold = 0.05,
                           tissue = miR.tissue_expression)
# Check the number of DEGs:
table(miR.Cp51_6h.res$sig_df$significance) # 1 up-regulated

# Save to excel
openxlsx::write.xlsx(miR.Cp51_6h.res[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Cp51_6h_miR_expr.xlsx"))


##############################################
# 3.) Visualise DESeq2 results               #
##############################################
(miR.Ca11_1h.res$plot <- miRNA_plot(miR.Ca11_1h.res$sig_df))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Ca11_1h_plot.png"),
       plot = miR.Ca11_1h.res$plot, width = 10, height = 8, units = 'in')

(miR.Ca11_6h.res$plot <- miRNA_plot(miR.Ca11_6h.res$sig_df))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Ca11_6h_plot.png"),
              plot = miR.Ca11_6h.res$plot, width = 10, height = 14, units = 'in')

(miR.Cp11_6h.res$plot <- miRNA_plot(miR.Cp11_6h.res$sig_df))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp11_6h_plot.png"),
       plot = miR.Cp11_6h.res$plot, width = 10, height = 10, units = 'in')

(miR.Cp51_6h.res$plot <- miRNA_plot(miR.Cp51_6h.res$sig_df))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp51_6h_plot.png"),
       plot = miR.Cp51_6h.res$plot, width = 10, height = 4, units = 'in')

##############################################
# 3.) Check miRNA target genes               #
##############################################
# C.albicans (1h)
miR.Ca11_1h.targets <- list(
  predicted = require_file(file.path(data_dir,miR.folder,"miRWalk_targets_CA11_1h_padj.csv"))
  )
miR.Ca11_1h.targets$validated <- miR.Ca11_1h.targets$predicted %>% 
  dplyr::filter(validated != "")

miR.Ca11_1h.targets <- lapply(miR.Ca11_1h.targets, function(x){
  x %>% 
    dplyr::mutate(binding_site = paste0("(",start,":",end,")", sep = ""),
                  validated = ifelse(validated == "", NA, validated),
                  position = factor(position, levels = c("3UTR","CDS","5UTR"))) %>% 
    dplyr::rename(TarBaseAccession = validated) %>% 
    dplyr::select(mirnaid, genesymbol, binding_site, position, TarBaseAccession) %>%
    group_by(mirnaid, genesymbol) %>%
    dplyr::summarise(
      binding_site = paste0(na.omit(binding_site), collapse = ", "),
      TarBaseAccession = paste0(na.omit(TarBaseAccession), collapse = ", "),
      position = paste0(na.omit(position), collapse = ", ")
    )
})

miR.Ca11_1h.targets_filtered <- lapply(miR.Ca11_1h.targets, function(x){
  x %>% 
    dplyr::inner_join(miR.Ca11_1h.res$sig_df[,c(7,1,2,5,6,8,9,10,11,12)], by = c("mirnaid" = "miRname")) %>% 
    dplyr::inner_join(Ca11_1h.res$sig_df[,c(1,2,5,6,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR, padj.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue, padj,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

openxlsx::write.xlsx(miR.Ca11_1h.targets_filtered[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Ca11_1h_miR_targets.xlsx"))

# C. albicans (6h)
miR.Ca11_6h.targets <- list(
  predicted = require_file(file.path(data_dir,miR.folder, "miRWalk_targets_CA11_6h_padj.csv"))
)
miR.Ca11_6h.targets$validated <- miR.Ca11_6h.targets$predicted %>% 
  dplyr::filter(validated != "")

miR.Ca11_6h.targets <- lapply(miR.Ca11_6h.targets, function(x){
  x %>% 
    dplyr::mutate(binding_site = paste0("(",start,":",end,")", sep = ""),
                  validated = ifelse(validated == "", NA, validated),
                  position = factor(position, levels = c("3UTR","CDS","5UTR"))) %>% 
    dplyr::rename(TarBaseAccession = validated) %>% 
    dplyr::select(mirnaid, genesymbol, binding_site, position, TarBaseAccession) %>%
    group_by(mirnaid, genesymbol) %>%
    dplyr::summarise(
      binding_site = paste0(na.omit(binding_site), collapse = ", "),
      TarBaseAccession = paste0(na.omit(TarBaseAccession), collapse = ", "),
      position = paste0(na.omit(position), collapse = ", ")
    )
})

miR.Ca11_6h.targets_filtered <- lapply(miR.Ca11_6h.targets, function(x){
  x %>% 
    dplyr::inner_join(miR.Ca11_6h.res$sig_df[,c(7,1,2,5,6,8,9,10,11,12)], by = c("mirnaid" = "miRname")) %>% 
    dplyr::inner_join(Ca11_6h.res$sig_df[,c(1,2,5,6,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR, padj.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue, padj,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

openxlsx::write.xlsx(miR.Ca11_6h.targets_filtered[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Ca11_6h_miR_targets.xlsx"))

# C.parapsilosis MOI 1:1 (6h)
miR.Cp11_6h.targets <- list(
  predicted = require_file(file.path(data_dir, miR.folder, "miRWalk_targets_CP11_6h_padj.csv"))
)

miR.Cp11_6h.targets$validated <- miR.Cp11_6h.targets$predicted %>% 
  dplyr::filter(validated != "")

miR.Cp11_6h.targets <- lapply(miR.Cp11_6h.targets, function(x){
  x %>% 
    dplyr::mutate(binding_site = paste0("(",start,":",end,")", sep = ""),
                  validated = ifelse(validated == "", NA, validated),
                  position = factor(position, levels = c("3UTR","CDS","5UTR"))) %>% 
    dplyr::rename(TarBaseAccession = validated) %>% 
    dplyr::select(mirnaid, genesymbol, binding_site, position, TarBaseAccession) %>%
    group_by(mirnaid, genesymbol) %>%
    dplyr::summarise(
      binding_site = paste0(na.omit(binding_site), collapse = ", "),
      TarBaseAccession = paste0(na.omit(TarBaseAccession), collapse = ", "),
      position = paste0(na.omit(position), collapse = ", ")
    )
})

miR.Cp11_6h.targets_filtered <- lapply(miR.Cp11_6h.targets, function(x){
  x %>% 
    dplyr::inner_join(miR.Cp11_6h.res$sig_df[,c(7,1,2,5,6,8,9,10,11,12)], by = c("mirnaid" = "miRname")) %>% 
    dplyr::inner_join(Cp11_6h.res$sig_df[,c(1,2,5,6,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR, padj.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue, padj,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

openxlsx::write.xlsx(miR.Cp11_6h.targets_filtered[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Cp11_6h_miR_targets.xlsx"))


# C.parapsilosis MOI 5:1 (6h)
miR.Cp51_6h.targets <- list(
  predicted = require_file(file.path(data_dir, miR.folder, "miRWalk_targets_CP51_6h_padj.csv"))
)

miR.Cp51_6h.targets$validated <- miR.Cp51_6h.targets$predicted %>% 
  dplyr::filter(validated != "")

miR.Cp51_6h.targets <- lapply(miR.Cp51_6h.targets, function(x){
  x %>% 
    dplyr::mutate(binding_site = paste0("(",start,":",end,")", sep = ""),
                  validated = ifelse(validated == "", NA, validated),
                  position = factor(position, levels = c("3UTR","CDS","5UTR"))) %>% 
    dplyr::rename(TarBaseAccession = validated) %>% 
    dplyr::select(mirnaid, genesymbol, binding_site, position, TarBaseAccession) %>%
    group_by(mirnaid, genesymbol) %>%
    dplyr::summarise(
      binding_site = paste0(na.omit(binding_site), collapse = ", "),
      TarBaseAccession = paste0(na.omit(TarBaseAccession), collapse = ", "),
      position = paste0(na.omit(position), collapse = ", ")
    )
})

miR.Cp51_6h.targets_filtered <- lapply(miR.Cp51_6h.targets, function(x){
  x %>% 
    dplyr::inner_join(miR.Cp51_6h.res$sig_df[,c(7,1,2,5,6,8,9,10,11,12)], by = c("mirnaid" = "miRname")) %>% 
    dplyr::inner_join(Cp51_6h.res$sig_df[,c(1,2,5,6,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR, padj.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue, padj,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

openxlsx::write.xlsx(miR.Cp51_6h.targets_filtered[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Cp51_6h_miR_targets.xlsx"))


##############################################
# 3.) Visualize miRNA targets                #
##############################################
(miR.Ca11_1h.plot <- get_CategoryExpressionPlot(results = Ca11_1h.res$df,
                                                sig_log2FC = 1.5, 
                                                sig_pval = 0.05, 
                                                targets = miR.Ca11_1h.targets_filtered,
                                                labels = "predicted"))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Ca11_1h_predicted_targets_plot.png"),
       plot = miR.Ca11_1h.plot, width = 10, height = 8, units = 'in')

(miR.Ca11_6h.plot <- get_CategoryExpressionPlot(results = Ca11_6h.res$df,
                                                sig_log2FC = 1.5, 
                                                sig_pval = 0.05, 
                                                targets = miR.Ca11_6h.targets_filtered,
                                                labels = "validated"))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Ca11_6h_validated_targets_plot.png"),
       plot = miR.Ca11_6h.plot, width = 10, height = 8, units = 'in')


(miR.Cp11_6h.plot <- get_CategoryExpressionPlot(results = Cp11_6h.res$df,
                                                sig_log2FC = 1.5, 
                                                sig_pval = 0.05, 
                                                targets = miR.Cp11_6h.targets_filtered,
                                                label = "predicted"))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp11_6h_predicted_targets_plot.png"),
       plot = miR.Cp11_6h.plot, width = 10, height = 8, units = 'in')

(miR.Cp51_6h.plot <- get_CategoryExpressionPlot(results = Cp51_6h.res$df,
                                                sig_log2FC = 1.5, 
                                                sig_pval = 0.05, 
                                                targets = miR.Cp51_6h.targets_filtered,
                                                label = "predicted"))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp51_6h_Predicted,targets_plot.png"),
       plot = miR.Cp51_6h.plot, width = 10, height = 8, units = 'in')

################################################################################
# 6.) Overrepresentation analysis (ORA) of miRNA target genes                  #
################################################################################
miR.Ca11_6h.genelist <- get_genelist(.df = Ca11_6h.res$df,
                                     .filter = Ca11_6h.res$df[["geneID"]] %in% 
                                      unique(miR.Ca11_6h.targets_filtered$predicted$genesymbol),
                                     .value = "stat",
                                     .name = "entrezID")
# Ca11 - 6h
miR.Ca11_6h.KEGG <- run_ora(.interest = miR.Ca11_6h.genelist$interest,
                            .background = miR.Ca11_6h.genelist$background,
                            .pathways = kegg_pathways)
miR.Ca11_6h.KEGG <- c(miR.Ca11_6h.KEGG,
                      extract_ora_results(.ora = miR.Ca11_6h.KEGG$ora,
                                          .db =  kegg_pathways))

miR.Ca11_6h.GO <- run_ora(.interest = miR.Ca11_6h.genelist$interest,
                          .background = miR.Ca11_6h.genelist$background,
                          .pathways = go_terms)
miR.Ca11_6h.GO <- c(miR.Ca11_6h.GO,
                    extract_ora_results(.ora = miR.Ca11_6h.GO$ora,
                                        .db = go_terms))

miR.Ca11_6h.GO <- c(miR.Ca11_6h.GO,
                    miR.Ca11_6h.GO$sig_df %>%
                      dplyr::group_split(Database) %>% 
                      setNames(c('BP','CC','MF')))

miR.Ca11_6h.simm <- list()
for (i in c('BP','CC','MF')) {
  miR.Ca11_6h.simm[[i]] <- make_simMatrix(miR.Ca11_6h.GO[[i]], 
                                          as.character(i), treshold = .9)
}

for (i in c('BP','CC','MF')) {
  # write similarity matrix
  write.xlsx(miR.Ca11_6h.simm[[i]]$simM,
             file.path(results_dir, miR.folder, date, tables_dir, 
                       paste0("miR_targets_Ca11_6h_",i,"_similarity_matrix.xlsx")))
  # write reduced similarity matrix
  write.xlsx(miR.Ca11_6h.simm[[i]]$reduced,
             file.path(results_dir, miR.folder, date, tables_dir, 
                       paste0("miR_targets_Ca11_6h_",i,"_reduced_clusters.xlsx")))
}

miR.Ca11_6h.GO_network <- list()
for (i in c('BP','CC','MF')) {
  miR.Ca11_6h.GO_network[[i]] <- get_cluster_representative(
    .cluster = inner_join(miR.Ca11_6h.GO[[i]], miR.Ca11_6h.simm[[i]]$reduced, 
                          by = c("ID" = "go")) %>% 
      dplyr::rename(Name = term, core_enrichment = geneID, geneRatio = GeneRatio),
    .degs = miR.Ca11_6h.targets_filtered$predicted %>% 
      dplyr::rename(geneID = genesymbol))
} 
#Save data frames to excel sheets
for (i in c('BP','CC','MF')) {
  # write similarity matrix
  write.xlsx(miR.Ca11_6h.GO_network[[i]], 
             file.path(results_dir, miR.folder, date, tables_dir, 
                       paste0("miR_targets_Ca11_6h_",i,"cluster_representatives.xlsx")))
}

# add GO semsim plot to each category
(miR.Ca11_6h.GO_network$BP_simplot <- make_GO_simplot(miR.Ca11_6h.GO_network$BP))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Ca11_6h_BP_simplot.png"),
       plot = miR.Ca11_6h.GO_network$BP_simplot, width = 10, height = 8, units = 'in')

(miR.Ca11_6h.GO_network$CC_simplot <- make_GO_simplot(miR.Ca11_6h.GO_network$CC))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Ca11_6h_CC_simplot.png"),
       plot = miR.Ca11_6h.GO_network$CC_simplot, width = 10, height = 8, units = 'in')

(miR.Ca11_6h.GO_network$MF_simplot <- make_GO_simplot(miR.Ca11_6h.GO_network$MF))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Ca11_6h_MF_simplot.png"),
       plot = miR.Ca11_6h.GO_network$MF_simplot, width = 10, height = 8, units = 'in')

# # C.parapsilosis MOI 1:1 (6h)
# miR.Cp11_6h.genelist <- get_genelist(.df = Cp11_6h.res$df,
#                                     .filter = Cp11_6h.res$df[["geneID"]] %in% 
#                                       miR.Cp11_6h.targets_filtered$predicted$genesymbol,
#                                     .value = "stat",
#                                     .name = "entrezID")
# 
# # No enriched term under specified p.adjust cutoff
# miR.Cp11_6h.KEGG <- run_ora(.interest = miR.Cp11_6h.genelist$interest,
#                            .background = miR.Cp11_6h.genelist$background,
#                            .pathways = kegg_pathways)
# miR.Cp11_6h.KEGG <- c(miR.Cp11_6h.KEGG$ora,
#                        extract_ora_results(.ora = miR.Cp11_6h.KEGG$ora,
#                                            .db =  kegg_pathways))
# 
# # No enriched term under specified p.adjust cutoff
# miR.Cp11_6h.GO <- run_ora(.interest = miR.Cp11_6h.genelist$interest,
#                          .background = miR.Cp11_6h.genelist$background,
#                          .pathways = go_terms)
# miR.Cp11_6h.GO <- c(miR.Cp11_6h.GO$ora,
#                     extract_ora_results(.ora = miR.Cp11_6h.GO$ora,
#                                         .db = go_terms))

# C.parapsilosis MOI 5:1 (6h)
# miR.Cp51_6h.genelist <- get_genelist(.df = Cp51_6h.res$df,
#                                     .filter = Cp51_6h.res$df[["geneID"]] %in% 
#                                       miR.Cp51_6h.targets_filtered$predicted$genesymbol,
#                                     .value = "stat",
#                                     .name = "entrezID")
# 
# miR.Cp51_6h.KEGG <- run_ora(.interest = miR.Cp51_6h.genelist$interest,
#                             .background = miR.Cp51_6h.genelist$background,
#                             .pathways = kegg_pathways)
# miR.Cp51_6h.KEGG <- c(miR.Cp51_6h.KEGG$ora,
#                       extract_ora_results(.ora = miR.Cp51_6h.KEGG$ora,
#                                           .db =  kegg_pathways))
# 
# miR.Cp51_6h.GO <- run_ora(.interest = miR.Cp51_6h.genelist$interest,
#                           .background = miR.Cp51_6h.genelist$background,
#                           .pathways = go_terms)
# miR.Cp51_6h.GO <- c(miR.Cp51_6h.GO$ora,
#                     extract_ora_results(.ora = miR.Cp51_6h.GO$ora,
#                                         .db = go_terms))
# 
# miR.Cp51_6h.GO <- c(miR.Cp51_6h.GO,
#                     miR.Cp51_6h.GO$sig_df %>%
#                       dplyr::group_split(Database) %>% 
#                       setNames(c('BP','CC','MF')))
# 
# miR.Cp51_6h.simm <- list()
# for (i in c('BP','CC','MF')) {
#   miR.Cp51_6h.simm[[i]] <- make_simMatrix(miR.Cp51_6h.GO[[i]], 
#                                           as.character(i), treshold = .9)
# }
# 
# miR.Cp51_6h.GO_network <- list()
# for (i in c('BP','CC','MF')) {
#   miR.Cp51_6h.GO_network[[i]] <- get_cluster_representative(
#     .cluster = inner_join(miR.Cp51_6h.GO[[i]], miR.Cp51_6h.simm[[i]]$reduced, 
#                           by = c("ID" = "go")) %>% 
#       dplyr::rename(Name = term, core_enrichment = geneID, geneRatio = GeneRatio),
#     .degs = miR.Cp51_6h.targets_filtered$predicted %>% 
#       dplyr::rename(Symbol = genesymbol))
# } 
# 
# # add GO semsim plot to each category
# (miR.Cp51_6h.GO_network$BP_simplot <- make_GO_simplot(miR.Cp51_6h.GO_network$BP))
# ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp51_6h_BP_simplot.png"),
#        plot = miR.Cp51_6h.GO_network$BP_simplot, width = 10, height = 8, units = 'in')
# 
# (miR.Cp51_6h.GO_network$CC_simplot <- make_GO_simplot(miR.Cp51_6h.GO_network$CC))
# ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp51_6h_CC_simplot.png"),
#        plot = miR.Cp51_6h.GO_network$CC_simplot, width = 10, height = 8, units = 'in')
# 
# (miR.Cp51_6h.GO_network$MF_simplot <- make_GO_simplot(miR.Cp51_6h.GO_network$MF))
# ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp51_6h_MF_simplot.png"),
#        plot = miR.Cp51_6h.GO_network$MF_simplot, width = 10, height = 8, units = 'in')


################################################################################
# 7.) Network analysis miRNA target genes                                      #
################################################################################
# C.albicans MOI 1:1 (6h)
predicted_target_ID <- Ca11_6h.cluster$df %>%
  dplyr::select(ID, Name, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::filter(core_enrichment %in% miR.Ca11_6h.targets_filtered$predicted$genesymbol) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(Name = first(Name), n = n()) %>%
  dplyr::arrange(desc(n)) %>% 
  dplyr::mutate(rank = row_number())

validated_target_ID <- Ca11_6h.cluster$df %>%
  dplyr::select(ID, Name, core_enrichment) %>%
  tidyr::separate_rows(core_enrichment, sep = "/") %>%
  dplyr::filter(core_enrichment %in% miR.Ca11_6h.targets_filtered$validated$genesymbol) %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(Name = first(Name), n = n()) %>%
  dplyr::arrange(desc(n)) %>% 
  dplyr::mutate(rank = row_number())

interest_ID_df <- dplyr::inner_join(predicted_target_ID, validated_target_ID, 
                                 by = c("ID", "Name"), suffix = c(".predicted",".validated")) %>%
  dplyr::mutate(rank = rank.predicted + rank.validated) %>%
  dplyr::arrange(rank)

write.xlsx(interest_ID_df,
           file.path(results_dir, miR.folder, date, tables_dir, "Ca11_6h_miR_targets_ID.xlsx"))

interest_ID <- dplyr::pull(interest_ID_df, Name, name = ID)[c(1:5)]
interest_genes <- unique(miR.Ca11_6h.targets_filtered$validated$genesymbol)

library(scales)
col_fun <- col_numeric(palette = c("blue", "white", "red"), domain = c(-1*max(Ca11_6h.combinedGSEA$NES),
                                                                      0, max(Ca11_6h.combinedGSEA$NES)))
colors <- col_fun(Ca11_6h.combinedGSEA[which(Ca11_6h.combinedGSEA$ID %in% names(interest_ID)),]$NES)


cluster_palette <- setNames(colors, 
                            c("TNFA_SIGNALING_VIA_NFKB","NEGATIVE_REGULATION_OF_CELL_POPULATION_PROLIFERATION",
                              "RESPONSE_TO_BACTERIUM","CELLULAR_RESPONSE_TO_LIPID",
                              "REGULATION_OF_PROTEIN_KINASE_ACTIVITY"))

# PANC1 - chronic acidosis
Ca11_6h.circplot <- getCircplotData(.cluster = Ca11_6h.cluster$df,
                                  .deg = Ca11_6h.res$sig_df,
                                  .interest_cluster = interest_ID, 
                                  .interest_cluster_genes = interest_genes, 
                                  .palette = cluster_palette)

plotCircplot(.path = file.path(results_dir, miR.folder, date, plots_dir, "Ca11_6h_miR_targets_circosplot.png"),
             .data = Ca11_6h.circplot$data.mat,
             .color = Ca11_6h.circplot$grid.col,
             .links = Ca11_6h.circplot$border.mat,
             .labels = c(interest_ID, interest_genes))
# Write data matrix
Ca11_6h.circplot_df <- merge(
  as.data.frame(Ca11_6h.circplot$data.mat) %>% dplyr::mutate(genesymbol = row.names(.)),
  as.data.frame(Ca11_6h.circplot$data.mat) %>% dplyr::mutate(genesymbol = row.names(.)),
  by = "genesymbol", suffixes = c("",".weight")) %>% 
  dplyr::mutate(across(TNFA_SIGNALING_VIA_NFKB:REGULATION_OF_PROTEIN_KINASE_ACTIVITY, 
                       ~ifelse(.x == 0, F, T))) %>% 
  dplyr::left_join(Ca11_6h.res$sig_df[,c("geneID","log2FoldChange","padj")],
                   by = c("genesymbol" = "geneID")) %>% 
  dplyr::left_join(miR.Ca11_6h.targets_filtered$validated[,c("genesymbol","miRname","log2FoldChange.miR", "padj.miR")],
                   by = c("genesymbol")) %>% 
  dplyr::relocate(genesymbol, log2FoldChange, padj, miRname, log2FoldChange.miR, padj.miR,
                  .before = 1)

write.xlsx(Ca11_6h.circplot_df,
           file.path(results_dir, miR.folder, date, tables_dir, 
                     "Ca11_6h_miR_targets_score_circosplot.xlsx"))

miR.Ca11_6h.barplot <- miR.Ca11_6h.targets_filtered$validated %>% 
  dplyr::filter(genesymbol %in% c("TRIB1","LDLR","BTG1","SOD2","MXD1", "PER2", "ABL2","SERTAD1")) %>%
  dplyr::select(miRname, genesymbol, log2FoldChange, log2FoldChange.miR) %>%
  tidyr::pivot_longer(cols = c(log2FoldChange, log2FoldChange.miR), 
               names_to = "Type",
               values_to = "log2FC") %>%
  dplyr::mutate(Type = ifelse(Type == "log2FoldChange.miR", "miRNA", "mRNA")) %>% 
  dplyr::mutate(genesymbol = ifelse(Type == "mRNA", genesymbol, miRname)) %>% 
  dplyr::mutate(genesymbol = factor(genesymbol, 
                                    levels = c("hsa-let-7a-3p","hsa-let-7b-3p","hsa-miR-26a-1-3p",
                                               "hsa-miR-27b-5p","hsa-miR-301a-5p","hsa-miR-30c-1-3p",
                                               "TRIB1","LDLR","BTG1","SOD2","MXD1",
                                               "PER2", "ABL2","SERTAD1")))

(miRNA_plots <- plot_miRNA_targets(miR.Ca11_6h.barplot))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Ca11_6h_targets_barplot.png"),
       plot = miRNA_plots, width = 6, height = 10, units = 'in')
