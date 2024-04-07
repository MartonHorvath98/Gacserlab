############################################## 
# 1.) Load analysis data                     #
##############################################
####################################################
# 1.) Set up working directory and directory tree  #
####################################################
# Set downstream path
miR.folder <- "miRNA"
if (!dir.exists(miR.folder)) {
  dir.create(miR.folder) # create the main results folder
}

# Create the sub folders for: results, data, and pictures
#data folder
data_dir <- "data"
if (!dir.exists(file.path(miR.folder,data_dir))) {
  dir.create(file.path(miR.folder,data_dir)) # create the data folder
} 

# Create the dated results folder
#plots directory
plots_dir <- "plots"
#results directory
results_dir <- "results"
# get the current date
date <- format(Sys.Date(), "%Y-%m-%d")
if (!dir.exists(file.path(miR.folder,date))) {
  dir.create(file.path(miR.folder, date)) # create the dated results folder
  dir.create(file.path(miR.folder, date, results_dir)) # create the results folder
  dir.create(file.path(miR.folder, date, plots_dir)) # create the plots folder
}

#######################################
# 2.) Load functions for the analyses #
#######################################
source("functions.R")
source("packages.R")

########################## 
# 3.) Load analysis data #
##########################
if (exists("miR.readcounts") == F) {
  miR.path <- c("F:/Labor/Projektek/HSC2_Reni/miRNA/mirdeep2/counts")
    #choose.dir(getwd(), "Select the directory containing the count files")
  miR.files <- list.files(miR.path, pattern = ".xlsx$", full.names = T)
  file.copy(from = miR.files, to = file.path(miR.folder, data_dir))
  
  miR.reads <- lapply(miR.files, read.xlsx, colNames = TRUE)
  
  miR.counts <- lapply(miR.reads, function(x) {
    x  %>% 
      select(c(1:2)) %>%
      setNames(c("miRname","Counts")) %>%
      dplyr::group_by(miRname) %>%
      dplyr::summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))
  })
  miR.counts <- merge.rec(miR.counts, by = "miRname", suffixes = c("",""))
  
  miR.names <- list.files(file.path(miR.folder, data_dir),
                           pattern = ".xlsx$", full.names = F) %>% 
    gsub(".xlsx", "", .) %>% 
    gsub("miRNAs_expressed_", "", .)
  
  colnames(miR.counts)[-1] <- miR.names
  write.table(miR.counts, file.path(miR.folder, data_dir,"miR_readcounts.csv"),
              sep =",", na = "NA", dec = ".", row.names = F, col.names = T)
  
}

##############################################
# 2.) Make differential analyses with DESeq2 #
##############################################
# Tissue specific mean expression of miRNAs
miR.tissue_expression <- list(
  "esophagus" = require_file(file.path(miR.folder, data_dir, "expression", "esophagus.csv")),
  "lymph_node" = require_file(file.path(miR.folder, data_dir, "expression", "lymph_node.csv")),
  "salivary_gland" = require_file(file.path(miR.folder, data_dir, "expression", "salivary_gland.csv")),
  "submandibular" = require_file(file.path(miR.folder, data_dir, "expression", "submandibular.csv")),
  "thyroid" = require_file(file.path(miR.folder, data_dir, "expression", "thyroid.csv")),
  "tongue" = require_file(file.path(miR.folder, data_dir, "expression", "tongue.csv"))
)

for(i in names(miR.tissue_expression)){
  miR.tissue_expression[[i]] <- miR.tissue_expression[[i]] %>%
    dplyr::select(!Expression.STD) %>% 
    dplyr::rename(miRname = "sncRNA") %>% 
    dplyr::mutate(tissue = i)
}

miR.tissue_expression <- do.call(rbind, miR.tissue_expression) %>% 
  dplyr::mutate(Expression.Mean = round(Expression.Mean, 3),
                miRname = as.factor(miRname),
                tissue = stringr::str_replace_all(tissue, "_", " "),
                tissue = paste0(toupper(substr(tissue, 1, 1)), substr(tissue, 2, nchar(tissue))),
                tissue = as.factor(tissue)) %>% 
  tidyr::pivot_wider(names_from = tissue, names_glue = "{tissue} (mean expr.)",
                     values_from = Expression.Mean)

# Load the count data
miR.readcounts <- read.csv(file = file.path(miR.folder, data_dir,"miR_readcounts.csv"),
                       sep = ",", header = T, na.strings = NA, row.names = 1)
miR.readcounts <- as.matrix(miR.readcounts)
# Split the data at the different times of infection
miR.readcounts <- list("1h" = miR.readcounts[,1:6], "6h" = miR.readcounts[,7:18])

# Run DESeq2 on the 1 hour infection model data
miR.DESeq.1h <- calc_DiffExp(matrix = miR.readcounts[["1h"]], 
                         coldata = coldata[["1h"]],
                         design = "condition")
# Run DESeq2 on the 6 hour infection model data
miR.DESeq.6h <- calc_DiffExp(matrix = miR.readcounts[["6h"]],
                         coldata = coldata[["6h"]],
                         design = "condition")

(miR.pca <- list(
  # PCA analysis of the 1 hour infection model
  "1h" = make_pca(miR.DESeq.1h$dds_norm, "samples", coldata[["1h"]]$condition),
  # PCA analysis of the 6 hour infection model
  "6h" = make_pca(miR.DESeq.6h$dds_norm, "samples", coldata[["6h"]]$condition)))
# Save the PCA plots
ggsave(file.path(miR.folder, date, plots_dir,"miR_1h_pca.png"),
       plot = miR.pca.1h, width = 8, height = 8, units = 'in')
ggsave(file.path(miR.folder, date, plots_dir,"miR_6h_pca.png"),
       plot = miR.pca.6h, width = 8, height = 8, units = 'in')

# Extract significant results for the 1 hour infection model
miR.Ca11_1h.res <- miR_results(miR.DESeq.1h$dds,
                           contrast = list("condition_CA11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1,
                           fdr_treshold = 0.1,
                           tissue = miR.tissue_expression)
# save to excel
sapply(names(miR.Ca11_1h.res), function(x){
  openxlsx::write.xlsx(miR.Ca11_1h.res[[x]],
                       file.path(miR.folder, date, results_dir,
                                 paste0("Ca11_1h_miR_",x,".xlsx")), 
                       rowNames = T)
})
# Extract significant results for the 6 hour infection model
## C.albicans MOI 1:1 (6h)
miR.Ca11_6h.res <- miR_results(miR.DESeq.6h$dds,
                            contrast = list("condition_CA11_vs_ctrl"),
                            name = c(''),
                            lfc_treshold = 1,
                            fdr_treshold = 0.1,
                            tissue = miR.tissue_expression)
sapply(names(miR.Ca11_6h.res), function(x){
  openxlsx::write.xlsx(miR.Ca11_6h.res[[x]], 
                       file.path(miR.folder, date, results_dir,
                                 paste0("Ca11_6h_miR_",x,".xlsx")),
                                 rowNames = T)
})
## C.parapsilosis MOI 1:1 (6h)
miR.Cp11_6h.res <- miR_results(miR.DESeq.6h$dds,
                           contrast = list("condition_CP11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1,
                           fdr_treshold = 0.1,
                           tissue = miR.tissue_expression)
sapply(names(miR.Cp11_6h.res), function(x){
  openxlsx::write.xlsx(miR.Cp11_6h.res[[x]], 
                       file.path(miR.folder, date, results_dir,
                                 paste0("Cp11_6h_miR_",x,".xlsx")),
                       rowNames = T)
})
## C.parapsilosis MOI 5:1 (6h)
miR.Cp51_6h.res <- miR_results(miR.DESeq.6h$dds,
                           contrast = list("condition_CP51_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1,
                           fdr_treshold = 0.1,
                           tissue = miR.tissue_expression)
sapply(names(miR.Cp51_6h.res), function(x){
  openxlsx::write.xlsx(miR.Cp51_6h.res[[x]], 
                       file.path(miR.folder, date, results_dir,
                                 paste0("Cp51_6h_miR_",x,".xlsx")),
                       rowNames = T)
})

##############################################
# 3.) Visualise DESeq2 results               #
##############################################

(miR.Ca11_1h.res$plot <- miRNA_plot(miR.Ca11_1h.res$sig_df))
ggsave(file.path(miR.folder, date, plots_dir,"miR_Ca11_1h_plot.png"),
       plot = miR.Ca11_1h.res$plot, width = 10, height = 8, units = 'in')

(miR.Ca11_6h.res$plot <- miRNA_plot(miR.Ca11_6h.res$sig_df))
ggsave(file.path(miR.folder, date, plots_dir,"miR_Ca11_6h_plot.png"),
              plot = miR.Ca11_6h.res$plot, width = 10, height = 8, units = 'in')

(miR.Cp11_6h.res$plot <- miRNA_plot(miR.Cp11_6h.res$sig_df))
ggsave(file.path(miR.folder, date, plots_dir,"miR_Cp11_6h_plot.png"),
       plot = miR.Cp11_6h.res$plot, width = 10, height = 8, units = 'in')

(miR.Cp51_6h.res$plot <- miRNA_plot(miR.Cp51_6h.res$sig_df))
ggsave(file.path(miR.folder, date, plots_dir,"miR_Cp51_6h_plot.png"),
       plot = miR.Cp51_6h.res$plot, width = 10, height = 3, units = 'in')

##############################################
# 3.) Check miRNA target genes               #
##############################################

# C.albicans (1h)
miR.Ca11_1h.targets <- list(
  predicted = require_file(file.path(miR.folder, data_dir, 
                                     "targets", "miRWalk_targets_CA11_1h.csv"))
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
    dplyr::inner_join(miR.Ca11_1h.res$sig_df[,c(7,1,2,5,9,10,11,12,13,14)], by = c("mirnaid" = "miRname")) %>% 
    dplyr::inner_join(Ca11_1h.res$sig_df[,c(1,2,5,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

sapply(names(miR.Ca11_1h.targets_filtered), function(x){
   openxlsx::write.xlsx(miR.Ca11_1h.targets_filtered[[x]], 
                        file.path(miR.folder, date, results_dir,
                                  paste0("Ca11_1h_miR_targets_",x,".xlsx")),
                        rowNames = T)
})

# C. albicans (6h)
miR.Ca11_6h.targets <- list(
  predicted = require_file(file.path(miR.folder, data_dir, 
                                     "targets", "miRWalk_targets_CA11_6h.csv"))
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
    dplyr::inner_join(miR.Ca11_6h.res$sig_df[,c(7,1,2,5,9,10,11,12,13,14)], by = c("mirnaid" = "miRname")) %>% 
    dplyr::inner_join(Ca11_6h.res$sig_df[,c(1,2,5,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

sapply(names(miR.Ca11_6h.targets_filtered), function(x){
  openxlsx::write.xlsx(miR.Ca11_6h.targets_filtered[[x]], 
                       file.path(miR.folder, date, results_dir,
                                 paste0("Ca11_6h_miR_targets_",x,".xlsx")),
                       rowNames = F)
})

# C.parapsilosis MOI 1:1 (6h)
miR.Cp11_6h.targets <- list(
  predicted = require_file(file.path(miR.folder, data_dir, 
                                     "targets", "miRWalk_targets_CP11_6h.csv"))
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
    dplyr::inner_join(miR.Cp11_6h.res$sig_df[,c(7,1,2,5,9,10,11,12,13,14)], by = c("mirnaid" = "miRname")) %>% 
    dplyr::inner_join(Cp11_6h.res$sig_df[,c(1,2,5,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

sapply(names(miR.Cp11_6h.targets_filtered), function(x){
  openxlsx::write.xlsx(miR.Cp11_6h.targets_filtered[[x]], 
                       file.path(miR.folder, date, results_dir,
                                 paste0("Cp11_6h_miR_targets_",x,".xlsx")),
                       rowNames = F)
})

# C.parapsilosis MOI 5:1 (6h)
miR.Cp51_6h.targets <- list(
  predicted = require_file(file.path(miR.folder, data_dir, 
                                     "targets", "miRWalk_targets_CP51_6h.csv"))
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
    dplyr::inner_join(miR.Cp51_6h.res$sig_df[,c(7,1,2,5,9,10,11,12,13,14)], by = c("mirnaid" = "miRname")) %>% 
    dplyr::inner_join(Cp51_6h.res$sig_df[,c(1,2,5,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})


sapply(names(miR.Cp51_6h.targets_filtered), function(x){
  openxlsx::write.xlsx(miR.Cp51_6h.targets_filtered[[x]], 
                       file.path(miR.folder, date, results_dir,
                                 paste0("Cp51_6h_miR_targets_",x,".xlsx")),
                       rowNames = F)
})

##############################################
# 3.) Visualize miRNA targets                #
##############################################
(miR.Ca11_1h.plot <- get_CategoryExpressionPlot(results = Ca11_1h.res$df,
                                                sig_log2FC = 1.5, 
                                                sig_pval = 0.05, 
                                                targets = miR.Ca11_1h.targets_filtered,
                                                labels = "predicted"))
ggsave(file.path(miR.folder, date, plots_dir,"miR_Ca11_1h_targets_plot.png"),
       plot = miR.Ca11_1h.plot, width = 10, height = 8, units = 'in')

(miR.Ca11_6h.plot <- get_CategoryExpressionPlot(results = Ca11_6h.res$df,
                                                sig_log2FC = 1.5, 
                                                sig_pval = 0.05, 
                                                targets = miR.Ca11_6h.targets_filtered,
                                                labels = "validated"))
ggsave(file.path(miR.folder, date, plots_dir,"miR_Ca11_6h_targets_plot.png"),
       plot = miR.Ca11_6h.plot, width = 10, height = 8, units = 'in')


(miR.Cp11_6h.plot <- get_CategoryExpressionPlot(results = Cp11_6h.res$df,
                                                sig_log2FC = 1.5, 
                                                sig_pval = 0.05, 
                                                targets = miR.Cp11_6h.targets_filtered,
                                                label = "predicted"))
ggsave(file.path(miR.folder, date, plots_dir,"miR_Cp11_6h_targets_plot.png"),
       plot = miR.Cp11_6h.plot, width = 10, height = 8, units = 'in')

(miR.Cp51_6h.plot <- get_CategoryExpressionPlot(results = Cp51_6h.res$df,
                                                sig_log2FC = 1.5, 
                                                sig_pval = 0.05, 
                                                targets = miR.Cp51_6h.targets_filtered,
                                                label = "predicted"))
ggplot2::ggsave(file.path(miR.folder, date, plots_dir,"miR_Cp51_6h_targets_plot.png"),
       plot = miR.Cp51_6h.plot, width = 10, height = 8, units = 'in')

##############################################
# 6.) GSEA with (validated) miRNA targets    #
##############################################

# C.albicans MOI 1:1 (6h)
miR.Ca11_6h.genelist <- subset(Ca11_6h.res$df, geneID %in% miR.Ca11_6h.targets_filtered$predicted$genesymbol) %>%
  na.omit(.) %>%
  dplyr::pull(log2FoldChange, name = entrezID)

miR.Ca11_6h.GSEA <- list()
miR.Ca11_6h.GSEA$KEGG <- gsea_results(miR.Ca11_6h.genelist, kegg_pathways)
miR.Ca11_6h.GSEA$GO <- gsea_results(miR.Ca11_6h.genelist, go_terms)
miR.Ca11_6h.GSEA$MSigDB <- gsea_results(miR.Ca11_6h.genelist, hallmark_sets)

miR.Ca11_6h.GO <- go_results(miR.Ca11_6h.genelist,
                             Ca11_6h.entrez$background,
                             type = "ALL")

miR.Ca11_6h.GO <- miR.Ca11_6h.GO %>%
  dplyr::group_split(ONTOLOGY) %>% 
  setNames(c('BP','MF'))

#names(miR.Ca11_6h.limma$go) <- c('BP','CC','MF')

miR.Ca11_6h.simm <- list()
for (i in c('BP','MF')) {
  miR.Ca11_6h.simm[[i]] <- make_simMatrix(miR.Ca11_6h.GO[[i]], 
                                          as.character(i), treshold = .9)
}

miR.Ca11_6h.simm$MF$subset <- miR.Ca11_6h.simm$MF$reduced[ miR.Ca11_6h.simm$MF$reduced$termDispensability < 0.1, ]


make_GO_simplot(miR.Ca11_6h.simm$MF$reduced, miR.Ca11_6h.simm$MF$subset)
