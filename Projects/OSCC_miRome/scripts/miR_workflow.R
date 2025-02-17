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
source("functions.R")
source("packages.R")

########################## 
# 3.) Load analysis data #
##########################
# Copy readcounts from their source folder to the data folder
miR.path <- choose.dir(getwd(), "Select the directory containing the count files")
miR.files <- list.files(miR.path, pattern = ".xlsx$", full.names = T)
file.copy(from = miR.files, to = file.path(data_dir, miR.folder))

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
  
}
rm(ls = ls()[grep("miR.", ls())])

##############################################
# 2.) Make differential analyses with DESeq2 #
##############################################
# Tissue specific mean expression of miRNAs
install.packages("arrow")
library(arrow)

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
                           lfc_treshold = 1.5,
                           p_treshold = 0.05,
                           tissue = miR.tissue_expression)
# save to excel
openxlsx::write.xlsx(miR.Ca11_1h.res[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Ca11_1h_miR_expr.xlsx"))

# Extract significant results for the 6 hour infection model
## C.albicans MOI 1:1 (6h)
miR.Ca11_6h.res <- miR_results(miR.DESeq.6h$dds,
                            contrast = list("condition_CA11_vs_ctrl"),
                            name = c(''),
                            lfc_treshold = 1.5,
                            p_treshold = 0.05,
                            tissue = miR.tissue_expression)
openxlsx::write.xlsx(miR.Ca11_6h.res[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Ca11_6h_miR_expr.xlsx"))

## C.parapsilosis MOI 1:1 (6h)
miR.Cp11_6h.res <- miR_results(miR.DESeq.6h$dds,
                           contrast = list("condition_CP11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1.5,
                           p_treshold = 0.05,
                           tissue = miR.tissue_expression)
openxlsx::write.xlsx(miR.Cp11_6h.res[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Cp11_6h_miR_expr.xlsx"))

## C.parapsilosis MOI 5:1 (6h)
miR.Cp51_6h.res <- miR_results(miR.DESeq.6h$dds,
                           contrast = list("condition_CP51_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1.5,
                           p_treshold = 0.05,
                           tissue = miR.tissue_expression)
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
              plot = miR.Ca11_6h.res$plot, width = 10, height = 8, units = 'in')

(miR.Cp11_6h.res$plot <- miRNA_plot(miR.Cp11_6h.res$sig_df))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp11_6h_plot.png"),
       plot = miR.Cp11_6h.res$plot, width = 10, height = 8, units = 'in')

(miR.Cp51_6h.res$plot <- miRNA_plot(miR.Cp51_6h.res$sig_df))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp51_6h_plot.png"),
       plot = miR.Cp51_6h.res$plot, width = 10, height = 8, units = 'in')

##############################################
# 3.) Check miRNA target genes               #
##############################################
# C.albicans (1h)
miR.Ca11_1h.targets <- list(
  predicted = require_file(file.path(data_dir,miR.folder,"miRWalk_targets_CA11_1h.csv"))
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
    dplyr::inner_join(Ca11_1h.res$sig_df[,c(1,2,5,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

openxlsx::write.xlsx(miR.Ca11_1h.targets_filtered[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Ca11_1h_miR_targets.xlsx"))

# C. albicans (6h)
miR.Ca11_6h.targets <- list(
  predicted = require_file(file.path(data_dir,miR.folder, "miRWalk_targets_CA11_6h.csv"))
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
    dplyr::inner_join(Ca11_6h.res$sig_df[,c(1,2,5,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

openxlsx::write.xlsx(miR.Ca11_6h.targets_filtered[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Ca11_6h_miR_targets.xlsx"))

# C.parapsilosis MOI 1:1 (6h)
miR.Cp11_6h.targets <- list(
  predicted = require_file(file.path(data_dir, miR.folder, "miRWalk_targets_CP11_6h.csv"))
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
    dplyr::inner_join(Cp11_6h.res$sig_df[,c(1,2,5,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

openxlsx::write.xlsx(miR.Cp11_6h.targets_filtered[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Cp11_6h_miR_targets.xlsx"))


# C.parapsilosis MOI 5:1 (6h)
miR.Cp51_6h.targets <- list(
  predicted = require_file(file.path(data_dir, miR.folder, "miRWalk_targets_CP51_6h.csv"))
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
    dplyr::inner_join(Cp51_6h.res$sig_df[,c(1,2,5,8)], by = c("genesymbol" = "geneID"),
                      suffix = c(".miR","")) %>% 
    dplyr::rename(miRname = mirnaid) %>%
    dplyr::relocate(miRname, baseMean.miR, log2FoldChange.miR, pvalue.miR,
                  genesymbol, baseMean, log2FoldChange, pvalue,
                  binding_site, TarBaseAccession, .before = 1) %>% 
    dplyr::filter(log2FoldChange > 0 & log2FoldChange.miR < 0 |
                    log2FoldChange < 0 & log2FoldChange.miR > 0)
})

openxlsx::write.xlsx(miR.Cp11_6h.targets_filtered[c(1:2)],
                     file.path(results_dir, miR.folder, date, tables_dir, "Cp11_6h_miR_targets.xlsx"))


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
                                      miR.Ca11_6h.targets_filtered$predicted$genesymbol,
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

miR.Ca11_6h.GO_network <- list()
for (i in c('BP','CC','MF')) {
  miR.Ca11_6h.GO_network[[i]] <- get_cluster_representative(
    .cluster = inner_join(miR.Ca11_6h.GO[[i]], miR.Ca11_6h.simm[[i]]$reduced, 
                          by = c("ID" = "go")) %>% 
      dplyr::rename(Name = term, core_enrichment = geneID, geneRatio = GeneRatio),
    .degs = miR.Ca11_6h.targets_filtered$predicted %>% 
      dplyr::rename(Symbol = genesymbol))
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

# C.parapsilosis MOI 1:1 (6h)
miR.Cp11_6h.genelist <- get_genelist(.df = Cp11_6h.res$df,
                                    .filter = Cp11_6h.res$df[["geneID"]] %in% 
                                      miR.Cp11_6h.targets_filtered$predicted$genesymbol,
                                    .value = "stat",
                                    .name = "entrezID")

# No enriched term under specified p.adjust cutoff
miR.Cp11_6h.KEGG <- run_ora(.interest = miR.Cp11_6h.genelist$interest,
                           .background = miR.Cp11_6h.genelist$background,
                           .pathways = kegg_pathways)
miR.Cp11_6h.KEGG <- c(miR.Cp11_6h.KEGG$ora,
                       extract_ora_results(.ora = miR.Cp11_6h.KEGG$ora,
                                           .db =  kegg_pathways))

# No enriched term under specified p.adjust cutoff
miR.Cp11_6h.GO <- run_ora(.interest = miR.Cp11_6h.genelist$interest,
                         .background = miR.Cp11_6h.genelist$background,
                         .pathways = go_terms)
miR.Cp11_6h.GO <- c(miR.Cp11_6h.GO$ora,
                    extract_ora_results(.ora = miR.Cp11_6h.GO$ora,
                                        .db = go_terms))

# C.parapsilosis MOI 5:1 (6h)
miR.Cp51_6h.genelist <- get_genelist(.df = Cp51_6h.res$df,
                                    .filter = Cp51_6h.res$df[["geneID"]] %in% 
                                      miR.Cp51_6h.targets_filtered$predicted$genesymbol,
                                    .value = "stat",
                                    .name = "entrezID")

miR.Cp51_6h.KEGG <- run_ora(.interest = miR.Cp51_6h.genelist$interest,
                            .background = miR.Cp51_6h.genelist$background,
                            .pathways = kegg_pathways)
miR.Cp51_6h.KEGG <- c(miR.Cp51_6h.KEGG$ora,
                      extract_ora_results(.ora = miR.Cp51_6h.KEGG$ora,
                                          .db =  kegg_pathways))

miR.Cp51_6h.GO <- run_ora(.interest = miR.Cp51_6h.genelist$interest,
                          .background = miR.Cp51_6h.genelist$background,
                          .pathways = go_terms)
miR.Cp51_6h.GO <- c(miR.Cp51_6h.GO$ora,
                    extract_ora_results(.ora = miR.Cp51_6h.GO$ora,
                                        .db = go_terms))

miR.Cp51_6h.GO <- c(miR.Cp51_6h.GO,
                    miR.Cp51_6h.GO$sig_df %>%
                      dplyr::group_split(Database) %>% 
                      setNames(c('BP','CC','MF')))

miR.Cp51_6h.simm <- list()
for (i in c('BP','CC','MF')) {
  miR.Cp51_6h.simm[[i]] <- make_simMatrix(miR.Cp51_6h.GO[[i]], 
                                          as.character(i), treshold = .9)
}

miR.Cp51_6h.GO_network <- list()
for (i in c('BP','CC','MF')) {
  miR.Cp51_6h.GO_network[[i]] <- get_cluster_representative(
    .cluster = inner_join(miR.Cp51_6h.GO[[i]], miR.Cp51_6h.simm[[i]]$reduced, 
                          by = c("ID" = "go")) %>% 
      dplyr::rename(Name = term, core_enrichment = geneID, geneRatio = GeneRatio),
    .degs = miR.Cp51_6h.targets_filtered$predicted %>% 
      dplyr::rename(Symbol = genesymbol))
} 

# add GO semsim plot to each category
(miR.Cp51_6h.GO_network$BP_simplot <- make_GO_simplot(miR.Cp51_6h.GO_network$BP))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp51_6h_BP_simplot.png"),
       plot = miR.Cp51_6h.GO_network$BP_simplot, width = 10, height = 8, units = 'in')

(miR.Cp51_6h.GO_network$CC_simplot <- make_GO_simplot(miR.Cp51_6h.GO_network$CC))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp51_6h_CC_simplot.png"),
       plot = miR.Cp51_6h.GO_network$CC_simplot, width = 10, height = 8, units = 'in')

(miR.Cp51_6h.GO_network$MF_simplot <- make_GO_simplot(miR.Cp51_6h.GO_network$MF))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Cp51_6h_MF_simplot.png"),
       plot = miR.Cp51_6h.GO_network$MF_simplot, width = 10, height = 8, units = 'in')


################################################################################
# 7.) Network analysis miRNA target genes                                      #
################################################################################

# interest_ID <- Ca11_6h.combinedGSEA %>% 
#   dplyr::select(ID, Name, core_enrichment) %>%
#   tidyr::separate_rows(core_enrichment, sep = "/") %>% 
#   dplyr::filter(core_enrichment %in% miR.Ca11_6h.targets_filtered$validated$genesymbol) %>%
#   dplyr::group_by(ID) %>% 
#   dplyr::summarise(Name = first(Name), n = n()) %>% 
#   dplyr::arrange(desc(n))

interest_ID <- c("GO:0002521" = "LEUKOCYTE_DIFFERENTIATION",
                 "GO:0070848" = "RESPONSE_TO_GROWTH_FACTOR",
                 "M5890" = "TNFA_SIGNALING_VIA_NFKB",
                 "GO:0001816" = "CYTOKINE_PRODUCTION",
                 "GO:0006954" = "INFLAMMATORY_RESPONSE")
interest_genes <- unique(miR.Ca11_6h.targets_filtered$validated$genesymbol)

library(scales)
col_fun <- col_numeric(palette = c("blue", "white", "red"), domain = range(Ca11_6h.combinedGSEA$NES))
colors <- col_fun(Ca11_6h.combinedGSEA[which(Ca11_6h.combinedGSEA$ID %in% names(interest_ID)),]$NES)


cluster_palette <- setNames(colors, 
                            c("LEUKOCYTE_DIFFERENTIATION","RESPONSE_TO_GROWTH_FACTOR",
                              "TNFA_SIGNALING_VIA_NFKB","CYTOKINE_PRODUCTION",
                              "INFLAMMATORY_RESPONSE"))

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

miR.Ca11_6h.barplot <- miR.Ca11_6h.targets_filtered$validated %>% 
  dplyr::filter(genesymbol %in% c("EPHA2","MAP2K3","RORA","TRIB1","VEGFA","TUBB2A")) %>%
  dplyr::select(miRname, genesymbol, log2FoldChange, log2FoldChange.miR) %>%
  tidyr::pivot_longer(cols = c(log2FoldChange, log2FoldChange.miR), 
               names_to = "Type",
               values_to = "log2FC") %>%
  dplyr::mutate(Type = ifelse(Type == "log2FoldChange.miR", "miRNA", "mRNA")) %>% 
  dplyr::mutate(genesymbol = ifelse(Type == "mRNA", genesymbol, miRname)) %>% 
  dplyr::mutate(genesymbol = factor(genesymbol, 
                                    levels = c("hsa-let-7c-3p","hsa-miR-1277-5p","hsa-miR-128-1-5p","hsa-miR-195-3p",
                                               "EPHA2","MAP2K3","RORA","TRIB1","VEGFA","TUBB2A")))

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

(miRNA_plots <- plot_miRNA_targets(miR.Ca11_6h.barplot))
ggsave(file.path(results_dir, miR.folder, date, plots_dir,"miR_Ca11_6h_targets_barplot.png"),
       plot = miRNA_plots, width = 6, height = 10, units = 'in')
