####################################################
# 1.) Set up working directory and directory tree  #
####################################################
# Set downstream path
folder <- "mRNA"
if (!dir.exists(folder)) {
  dir.create(folder) # create the main results folder
}

# Create the sub folders for: results, data, and pictures
#data folder
data_dir <- "data"
if (!dir.exists(file.path(folder,data_dir))) {
  dir.create(file.path(folder,data_dir)) # create the data folder
} 

# Create the dated results folder
#plots directory
plots_dir <- "plots"
#results directory
results_dir <- "results"
# get the current date
date <- format(Sys.Date(), "%Y-%m-%d")
if (!dir.exists(file.path(folder,date))) {
  dir.create(file.path(folder, date)) # create the dated results folder
  dir.create(file.path(folder, date, results_dir)) # create the results folder
  dir.create(file.path(folder, date, plots_dir)) # create the plots folder
}

#######################################
# 2.) Load functions for the analyses #
#######################################
source("functions.R")
source("packages.R")

########################## 
# 3.) Load analysis data #
##########################
if (exists("readcounts") == F) {
  mrna.path <- choose.dir(getwd(), "Select the directory containing the count files")
  mrna.files <- list.files(mrna.path, pattern = "counts.txt$", full.names = T)
  file.copy(from = mrna.files, to = file.path(folder, data_dir))
  
  mrna.reads <- lapply(list.files(file.path(folder, data_dir), 
                                  pattern = "counts.txt$",
                                  full.names = T), 
                       # read csv: with header, tab-delimited, skip first row =>
                       # contains quantification parameters
                       read.csv, header = T, sep = "\t", skip = 1)
  
  mrna.counts <- lapply(mrna.reads, 
                        function(x) { x  %>%
                            select(c(1,7)) %>%
                            setNames(c("Geneid","Counts"))}
  )
  
  mrna.counts <- merge.rec(mrna.counts, by = "Geneid",  all = T, suffixes = c("",""))
  file.names <- list.files(file.path(folder, data_dir),
                           pattern = "counts.txt$", full.names = F)
  
  mrna.counts <- mrna.counts %>%
    # set the gene ID as row names
    tibble::column_to_rownames("Geneid") %>% 
    # clean up the column names
    setNames(str_remove(file.names, "_counts.txt"))
  
  write.table(mrna.counts, paste(data_dir,"readcounts.csv", sep = "/"), 
              sep =",", na = "NA", dec = ".", row.names = T, col.names = T)
}

##########################
# 4.) Ready count tables #
##########################
readcounts <- read.csv(file = file.path(folder, data_dir,"readcounts.csv"),
                       sep = ",", header = T, na.strings = NA, row.names = 1)
readcounts <- as.matrix(readcounts)

# split the readcounts into 1h and 6h
readcounts <- list("1h" = readcounts[,1:6], "6h" = readcounts[,7:18])

##############################
# 5.) Prepare metadata table #
##############################

# Create the metadata table
coldata <- data.frame("samples" = file.names) %>%
  dplyr::mutate(time = dplyr::case_when(
    stringr::str_detect(file.names, "H1") ~ "1h",
    stringr::str_detect(file.names, "H6") ~ "6h"),
    time = factor(time, 
                  levels = c("1h","6h"))) %>% 
  dplyr::group_split(time, .keep = F)

coldata <- lapply(coldata, function(x){
  x %>%
    dplyr::mutate(
      samples = as.factor(samples),
      # extract the treatment from the sample name
      treatment = dplyr::case_when(
        stringr::str_detect(samples, "_CA") ~ "CA",
        stringr::str_detect(samples, "_CP") ~ "CP",
        T ~ "ctrl"),
      treatment = relevel(as.factor(treatment), ref = "ctrl"),
      # extract the moiety from the sample name
      moi = dplyr::case_when(
        stringr::str_detect(samples, "_CA11") ~ "11",
        stringr::str_detect(samples, "_CP11") ~ "11",
        stringr::str_detect(samples, "_CP51") ~ "51",
        T ~ ""),
      moi = relevel(as.factor(moi), ref = ""),
      # create the condition column gluing the treatment and moi
      condition = dplyr::case_when(
        stringr::str_detect(treatment, "ctrl") ~ stringr::str_glue("{treatment}"),
        T ~ stringr::str_glue("{treatment}{moi}")
      ),
      condition = factor(condition),
      condition = relevel(condition, ref = "ctrl"))
})

names(coldata) <- c("1h","6h")

##############################################
# 6.) Make differential analyses with DESeq2 #
##############################################

# Run DESeq2 on the 1 hour infection model data
DESeq.1h <- calc_DiffExp(matrix = readcounts[["1h"]], 
                         coldata = coldata[["1h"]],
                         design = "condition")
# Run DESeq2 on the 6 hours infection model data
DESeq.6h <- calc_DiffExp(matrix = readcounts[["6h"]],
                         coldata = coldata[["6h"]],
                         design = "condition")

# PCA analysis of the 1 hour infection model data 
(pca.1h <- make_pca(DESeq.1h$dds_norm, "samples", coldata[["1h"]]$condition))
ggsave(file.path(folder, date, plots_dir, "1h_pca.png"), plot = pca.1h,
       width = 8, height = 8, units = 'in')

# PCA analysis of the 6 hours infection model data
(pca.6h <- make_pca(DESeq.6h$dds_norm, "samples", coldata[["6h"]]$condition))
ggsave(file.path(folder, date, plots_dir, "6h_pca.png"), plot = pca.6h,
       width = 8, height = 8, units = 'in')

# Extract significant results for the 1 hour infection model
Ca11_1h.res <- get_results(DESeq.1h$dds,
                           contrast = list("condition_CA11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1.5,
                           pval_treshold = 0.05)
# Save results to excel files
sapply(names(Ca11_1h.res), function(x){
  openxlsx::write.xlsx(Ca11_1h.res[[x]], 
                       file.path(folder, date, results_dir, 
                                 paste0("Ca11_1h_",x,".xlsx")),
                       rowNames = T)
})
# Extract significant results for the 6 hours infection model
## 1. Candida albicans - MOI 1:1
Ca11_6h.res <- get_results(DESeq.6h$dds,
                           contrast = list("condition_CA11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1.5,
                           pval_treshold = 0.05)
# Save results to excel files
sapply(names(Ca11_6h.res), function(x){
  openxlsx::write.xlsx(Ca11_6h.res[[x]], 
                       file.path(folder, date, results_dir, 
                                 paste0("Ca11_6h_",x,".xlsx")),
                       rowNames = T)
})
## 2. Candida parapsilosis - MOI 1:1
Cp11_6h.res <- get_results(DESeq.6h$dds,
                           contrast = list("condition_CP11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1.5,
                           pval_treshold = 0.05)

sapply(names(Cp11_6h.res), function(x){
  openxlsx::write.xlsx(Cp11_6h.res[[x]], 
                       file.path(folder, date, results_dir, 
                                 paste0("Cp11_6h_",x,".xlsx")),
                       rowNames = T)
})
## 3. Candida parapsilosis - MOI 5:1
Cp51_6h.res <- get_results(DESeq.6h$dds,
                           contrast = list("condition_CP51_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1.5,
                           pval_treshold = 0.05)

sapply(names(Cp51_6h.res), function(x){
  openxlsx::write.xlsx(Cp51_6h.res[[x]], 
                       file.path(folder, date, results_dir, 
                                 paste0("Cp51_6h_",x,".xlsx")),
                       rowNames = T)
})

##############################################
# 7.) Visualise DESeq2 results               #
##############################################
# 1. Vulcano plots
## 1.1 1 hour infection model
(Ca11_1h.res$plot <- make_vulcanplot(Ca11_1h.res$df, Ca11_1h.res$sig_df))

# save the plot
ggsave(file.path(folder, date, plots_dir, "Ca11_1h_vulcan.png"), 
       plot = Ca11_1h.res$plot, width = 10, height = 8, units = 'in')

## 1.2 6 hours infection model
# CA11
(Ca11_6h.res$plot <- make_vulcanplot(Ca11_6h.res$df, Ca11_6h.res$sig_df))
ggsave(file.path(folder, date, plots_dir, "Ca11_6h_vulcan.png"), 
       plot = Ca11_6h.res$plot, width = 10, height = 8, units = 'in')

# CP11
(Cp11_6h.res$plot <- make_vulcanplot(Cp11_6h.res$df, Cp11_6h.res$sig_df))
ggsave(file.path(folder, date, plots_dir, "Cp11_6h_vulcan.png"), 
       plot = Cp11_6h.res$plot, width = 10, height = 8, units = 'in')

# CP51
(Cp51_6h.res$plot <- make_vulcanplot(Cp51_6h.res$df, Cp51_6h.res$sig_df))
ggsave(file.path(folder, date, plots_dir, "Cp51_6h_vulcan.png"), 
       plot = Cp51_6h.res$plot, width = 10, height = 8, units = 'in')

# 2. Heatmaps
## 2.1 1 hour infection model
(heatmap.1h <- make_heatmap(DESeq.1h$dds_norm, "condition"))
ggsave(file.path(folder, date, plots_dir, "1h_heatmap.png"), 
       plot = heatmap.1h, width = 8, height = 10, units = 'in')

## 2.2 6 hours infection model
(heatmap.6h <- make_heatmap(DESeq.6h$dds_norm, "condition"))
ggsave(file.path(folder, date, plots_dir, "6h_heatmap.png"), 
       plot = heatmap.6h, width = 8, height = 10, units = 'in')

##############################################
# 8.) Overrepresentation analyses            #
##############################################

# 1 KEGG pathways - ORA analysis
## 1.1 Get entrez IDs: gene set of interest and background
Ca11_1h.entrez <- get_genelist(Ca11_1h.res) # CA11 1 hour
Ca11_6h.entrez <- get_genelist(Ca11_6h.res) # CA11 6 hours
Cp11_6h.entrez <- get_genelist(Cp11_6h.res) # CP11 6 hours
Cp51_6h.entrez <- get_genelist(Cp51_6h.res) # CP51 6 hours

## 1.2 Run the ORA analysis
Ca11_1h.KEGG <- list()
Ca11_1h.KEGG$df <- kegg_results(Ca11_1h.entrez$interest,
                             Ca11_1h.entrez$background)
openxlsx::write.xlsx(Ca11_1h.KEGG$df, 
                     file.path(folder, date, results_dir, "Ca11_1h_KEGG.xlsx"))

Ca11_6h.KEGG <- list()
Ca11_6h.KEGG$df <- kegg_results(Ca11_6h.entrez$interest,
                             Ca11_6h.entrez$background)
openxlsx::write.xlsx(Ca11_6h.KEGG$df,
                     file.path(folder, date, results_dir, "Ca11_6h_KEGG.xlsx"))

Cp11_6h.KEGG <- list()
Cp11_6h.KEGG$df <- kegg_results(Cp11_6h.entrez$interest,
                             Cp11_6h.entrez$background) # no results in CP11 H6
# openxlsx::write.xlsx(Cp11_6h_KEGG$df,
#                      file.path(folder, date, results_dir, "Cp11_6h_KEGG.xlsx"))

Cp51_6h.KEGG <- list()
Cp51_6h.KEGG$df <- kegg_results(Cp51_6h.entrez$interest,
                             Cp51_6h.entrez$background)
openxlsx::write.xlsx(Cp51_6h.KEGG$df,
                     file.path(folder, date, results_dir, "Cp51_6h_KEGG.xlsx"))

## 1.3 Graphical analysis
(Ca11_1h.KEGG$plot <- make_dotplot(Ca11_1h.KEGG$df, Ca11_1h.KEGG$df$Count, type = "KEGG"))
ggsave(file.path(folder, date, plots_dir,"Ca11_1h_KEGG.png"), 
       plot = Ca11_1h.KEGG$plot, width = 10, height = 8, units = 'in')

(Ca11_6h.KEGG$plot <- make_dotplot(Ca11_6h.KEGG$df, Ca11_6h.KEGG$df$Count, type = "KEGG"))
ggsave(file.path(folder, date, plots_dir,"Ca11_6h_KEGG.png"), 
       plot = Ca11_6h.KEGG$plot, width = 10, height = 8, units = 'in')

# (Cp11_6h.KEGG$plot <- make_dotplot(Cp11_6h.KEGG$df, Cp11_6h.KEGG$df$Count, type = "KEGG"))
# ggsave(file.path(folder, date, plots_dir,"Cp11_6h_KEGG.png"), 
#        plot = Cp11_6h.KEGG$plot, width = 10, height = 8, units = 'in')
# 0 hits

(Cp51_6h.KEGG$plot <- make_dotplot(Cp51_6h.KEGG$df, Cp51_6h.KEGG$df$Count, type = "KEGG"))
ggsave(file.path(folder, date, plots_dir,"Cp51_6h_KEGG.png"), 
       plot = Cp51_6h.KEGG$plot, width = 10, height = 8, units = 'in')

# 2 GO terms - ORA analysis
## 2.1 Run the ORA analysis
Ca11_1h.GO <- list()
# Ca11_1h.GO$df <- go_results(Ca11_1h.entrez$interest,
#                              Ca11_1h.entrez$background,
#                              type = 'ALL')
# openxlsx::write.xlsx(Ca11_1h_GO_ORA, 
#                      file.path(folder, date, results_dir, "Ca11_1h_GO.xlsx"))
# 0 hits

Ca11_6h.GO <- list()
Ca11_6h.GO$df <- go_results(Ca11_6h.entrez$interest,
                            Ca11_6h.entrez$background,
                            type = 'ALL')
openxlsx::write.xlsx(Ca11_6h.GO$df,
                     file.path(folder, date, results_dir, "Ca11_6h_GO.xlsx"))

Cp11_6h.GO <- list()
Cp11_6h.GO$df <- go_results(Cp11_6h.entrez$interest,
                            Cp11_6h.entrez$background,
                            type = 'ALL')
openxlsx::write.xlsx(Cp11_6h.GO$df,
                     file.path(folder, date, results_dir, "Cp11_6h_GO.xlsx"))

Cp51_6h.GO <- list()
Cp51_6h.GO$df <- go_results(Cp51_6h.entrez$interest,
                            Cp51_6h.entrez$background,
                            type = 'ALL')
openxlsx::write.xlsx(Cp51_6h.GO$df,
                     file.path(folder, date, results_dir, "Cp51_6h_GO.xlsx"))

## 2.2 Graphical analysis

# CA11 1 hour - 0 hits

(Ca11_6h.GO$plot <- make_dotplot(Ca11_6h.GO$df, Ca11_6h.GO$df$Count, type = "GO"))
ggsave(file.path(folder, date, plots_dir,"Ca11_6h_GO.png"), 
       plot = Ca11_6h.GO$plot, width = 10, height = 8, units = 'in')

(Cp11_6h.GO$plot <- make_dotplot(Cp11_6h.GO$df, Cp11_6h.GO$df$Count, type = "GO"))
ggsave(file.path(folder, date, plots_dir,"Cp11_6h_GO.png"), 
       plot = Cp11_6h.GO$plot, width = 10, height = 8, units = 'in')

(Cp51_6h.GO$plot <- make_dotplot(Cp51_6h.GO$df, Cp51_6h.GO$df$Count, type = "GO"))
ggsave(file.path(folder, date, plots_dir,"Cp51_6h_GO.png"), 
       plot = Cp51_6h.GO$plot, width = 10, height = 8, units = 'in')

##############################################
# 9.) Gene set enrichment analyses           #
##############################################
msigdbr_df <- msigdbr(species = "Homo sapiens")

# create reference databases
kegg_pathways <- msigdbr_df %>%
  dplyr::filter(
    gs_cat == "C2", # only canonical representations (compiled by experts)
    gs_subcat == "CP:KEGG" # KEGG pathways
  )
  
got_terms <- msigdbr_df %>%
  dplyr::filter(gs_cat == "C5" &  # ontology gene sets
                  gs_subcat != "HPO") # excluding human ontology phenotype gs

hallmark_sets <- msigdbr_df %>%
  dplyr::filter(gs_cat == "H")  # only hallmark gene sets
  
## Ca11 1h - no term enriched under specific pvalueCutoff
#####
# Ca11_1h.GSEA <- list()
# Ca11_1h.GSEA$KEGG <- gsea_results(Ca11_1h.entrez$interest, kegg_pathways)
# Ca11_1h.GSEA$GO <- gsea_results(Ca11_1h.entrez$interest, go_terms) 
# Ca11_1h.GSEA$MSigDB <- gsea_results(Ca11_1h.entrez$interest, hallmark_sets) 
# 
# Ca11_1h.GSEA.df <- list(
#   kegg = Ca11_1h.GSEA$KEGG$gsea_df, # GSEA against the KEGG pathway datasets
#   GO = Ca11_1h.GSEA$GO$gsea_df,     # GSEA against the GO term datasets
#   hallmark = Ca11_1h.GSEA$MSigDB$gsea_df # GSEA against the MSigDB halmark sets
# )
# 
# # save files to excel
# sapply(names(Ca11_1h.GSEA.df), function(x){
#   openxlsx::write.xlsx(Ca11_1h.GSEA.df[[x]], 
#                        file.path(folder, date, results_dir, 
#                                  paste0("Ca11_1h_GSEA_",x,".xlsx")),
#                        rowNames = T)
# })
#####
## Ca11 6h
Ca11_6h.GSEA <- list()
Ca11_6h.GSEA$KEGG <- gsea_results(Ca11_6h.entrez$interest, kegg_pathways)
Ca11_6h.GSEA$GO <- gsea_results(Ca11_6h.entrez$interest, go_terms)
Ca11_6h.GSEA$MSigDB <- gsea_results(Ca11_6h.entrez$interest, hallmark_sets)

Ca11_6h.GSEA.df <- list(
  kegg = Ca11_6h.GSEA$KEGG$gsea_df, # GSEA against the KEGG pathway datasets
  GO = Ca11_6h.GSEA$GO$gsea_df,     # GSEA against the GO term datasets
  hallmark = Ca11_6h.GSEA$MSigDB$gsea_df # GSEA against the MSigDB halmark sets
)

# save files to excel
sapply(names(Ca11_6h.GSEA.df), function(x){
  openxlsx::write.xlsx(Ca11_6h.GSEA.df[[x]], 
                       file.path(folder, date, results_dir, 
                                 paste0("Ca11_6h_GSEA_",x,".xlsx")),
                       rowNames = T)
})

## Cp11 6h - no terms were enriched under specific pvalueCutoff
#####
# Cp11_6h.GSEA <- list()
# Cp11_6h.GSEA$KEGG <- gsea_results(Cp11_6h.entrez$interest, kegg_pathways)
# Cp11_6h.GSEA$GO <- gsea_results(Cp11_6h.entrez$interest, go_terms)
# Cp11_6h.GSEA$MSigDB <- gsea_results(Cp11_6h.entrez$interest, hallmark_sets)
# 
# Cp11_6h.GSEA.df <- list(
#   kegg = Cp11_6h.GSEA$KEGG$gsea_df, # GSEA against the KEGG pathway datasets
#   GO = Cp11_6h.GSEA$GO$gsea_df,     # GSEA against the GO term datasets
#   hallmark = Cp11_6h.GSEA$MSigDB$gsea_df # GSEA against the MSigDB halmark sets
# )
# 
# # save files to excel
# sapply(names(Cp11_6h.GSEA.df), function(x){
#   openxlsx::write.xlsx(Cp11_6h.GSEA.df[[x]], 
#                        file.path(folder, date, results_dir, 
#                                  paste0("Cp11_6h_GSEA_",x,".xlsx")),
#                        rowNames = T)
# })
#####
## Cp51 6h - no terms were enriched under specific pvalueCutoff
#####
# Cp51_6h.GSEA <- list()
# Cp51_6h.GSEA$KEGG <- gsea_results(Cp51_6h.entrez$interest, kegg_pathways)
# Cp51_6h.GSEA$GO <- gsea_results(Cp51_6h.entrez$interest, go_terms)
# Cp51_6h.GSEA$MSigDB <- gsea_results(Cp51_6h.entrez$interest, hallmark_sets)
# 
# Cp51_6h.GSEA.df <- list(
#   kegg = Cp51_6h.GSEA$KEGG$gsea_df, # GSEA against the KEGG pathway datasets
#   GO = Cp51_6h.GSEA$GO$gsea_df,     # GSEA against the GO term datasets
#   hallmark = Cp51_6h.GSEA$MSigDB$gsea_df # GSEA against the MSigDB halmark sets
# )
# 
# # save files to excel
# sapply(names(Cp51_6h.GSEA.df), function(x){
#   openxlsx::write.xlsx(Cp51_6h.GSEA.df[[x]], 
#                        file.path(folder, date, results_dir, 
#                                  paste0("Cp51_6h_GSEA_",x,".xlsx")),
#                        rowNames = T)
# })

##############################################
# 10.) Overlap with cancer hallmark genesets #
##############################################
if (exists("pancancer.genes") == F) {
  pancancer.path <- choose.files(getwd(), "Select the directory containing the count files",
                                 multi = F)
  file.copy(from = pancancer.path, to = file.path(folder, data_dir),
            overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
}

pancancer.file <- list.files(file.path(folder, data_dir), pattern = "cancer",
                             full.names = F)
  
pancancer.df <- require_file(file.path(folder, data_dir, pancancer.file),
                             header = T, na.strings = NA)
pancancer.df <- pancancer.df %>%
  dplyr::mutate(across(where(is.character) & !c(Genes), ~ifelse(.x == '+',T,F)))

pancancer.category <- colnames(pancancer.df)[-1]

pancancer.genes <- list()
for(i in pancancer.category){
  pancancer.genes[[i]] <- pancancer.df %>%
    .[which(.[,i] == T),] %>%
    dplyr::select(Genes) %>%
    dplyr::rename("geneID" = Genes) %>%
    dplyr::mutate(PANC = TRUE)
}


pancancer.res <- lapply(pancancer.genes, function(x){
  return(merge.rec(list(Ca11_6h.res$sig_df, Cp51_6h.res$sig_df, x), 
                   by = "geneID", all = T, suffixes=c(".CA11",".CP51")))
})

pancancer.upset <- lapply(pancancer.res, function(x) {
  x %>% 
    # select the geneID and the significance columns
    dplyr::select(c(1,8,16,18)) %>% 
    # add columns for regions
    dplyr::mutate(CA11 = ifelse(is.na(significance.CA11), FALSE, TRUE)) %>%
    dplyr::mutate(CP51 = ifelse(is.na(significance.CP51), FALSE, TRUE)) %>%
    dplyr::mutate(PANC = ifelse(is.na(PANC), FALSE, TRUE)) %>%
    relocate(where(is.logical), .before = where(is.character)) %>%
    # add columns for trends in expression changes
    dplyr::mutate(Trend = case_when(
      significance.CA11 == 'Signif. down-regulated' | significance.CP51 == 'Signif. down-regulated' ~ 'Signif. down-regulated',
      significance.CA11 == 'Signif. up-regulated' | significance.CP51 == 'Signif. up-regulated' ~ 'Signif. up-regulated',
      significance.CA11 == 'Signif. down-regulated' & significance.CP51 == 'Signif. up-regulated' |
        significance.CA11 == 'Signif. up-regulated' & significance.CP51 == 'Signif. down-regulated' ~ 'Changed regulation',
      TRUE ~ NA
    ),
    Trend = as.factor(Trend)) %>%
    # remove significance columns
    dplyr::select(-c(5,6))}
)

pancancer.arranged <- list()
for(i in names(pancancer.upset)){
  pancancer.arranged[[i]] <- arrange_venn(pancancer.upset[[i]],
                                               c("CA11", "CP51", "PANC"),
                                          extract_regions = T) %>% 
    dplyr::arrange(., region)
}


##############################################
# 10.) Venn diagrams with cancer hallmark    #
##############################################

# Create the venn diagrams' base data frames

pancancer.vennbase <- list()
for(i in names(pancancer.arranged)){
  pancancer.vennbase[[i]] <- make_upsetvennbase(arranged = pancancer.arranged[[i]],
                                                df = pancancer.upset[[i]])
}

# Create the venn diagrams
pancancer.labels <- c(
  "Angiogenesis" = "Angiogenic genes",
  "Cancer.Metabolism" = "Metabolic genes",
  "ECM.Layers" = "ECM layers",
  "ECM.Remodeling" = "ECM remodeling",
  "EMT" = "Epithelial-mesenchymal transition",
  "Hypoxia" = "Hypoxia-related genes",             
  "Metastasis" = "Metastatic genes",
  "Transcription.Factor" = "Tumor transcription factors (TFs)",
  "Tumor.Growth" = "Tumor growth",
  "Tumor.Invasion" = "Tumor invasion"
)

pancancer.venn <- list()
for(i in names(pancancer.vennbase)){
  pancancer.venn[[i]] <- make_upsetvenn(data = pancancer.vennbase[[i]],
                                        sets = c("CA11", "CP51", "PANC"),
                                        labels = c(
                                          "C.albicans (MOI 1:1)",
                                          "C.parapsilosis (MOI 5:1)",
                                          pancancer.labels[i]))
}

# Save the venn diagrams
sapply(names(pancancer.venn), function(x){
  ggsave(file.path(folder, date, plots_dir, paste0("HSC2_vs_",x,".png")),
         plot = pancancer.venn[[x]], width = 16, height = 10, units = 'in')
})


ggsave(paste(plots_dir,"/HSC2_vs_Angiogenesis.png"), plot = HSC2_Angiogenesis.plot,
       width = 16, height = 10, units = 'in')

# Extract the geneIDs from the venn diagrams' subsets
HSC2_vs_pancancer <- list()
for(i in names(pancancer.vennbase)){
  HSC2_vs_pancancer[[i]] <- pancancer.vennbase[[i]] %>%
    dplyr::filter(!region %in% c("CA11","CP51","PANC","CA11-CP51")) %>%
    dplyr::pull(geneID, name = region)
}

##############################################
# 11.) Compare C.albi & C.para cancer sets   #
##############################################

# [1] "Angiogenesis"         "Cancer.Metabolism"    "ECM.Layers"           "ECM.Remodeling"      
# [5] "EMT"                  "Hypoxia"              "Metastasis"           "Transcription.Factor"
# [9] "Tumor.Growth"         "Tumor.Invasion" 

# Calculate gene ratio and activation z-score
pancancer.spider <- lapply(pancancer.upset, make_spiderbase)

# Create data frame for the spider plot
pancancer.circ <- lapply(names(pancancer.spider), function(x){
  pancancer.spider[[x]] %>% 
    dplyr::mutate(category = as.character(x),
                  category = as.factor(str_replace_all(category,"[.]"," "))) %>%
    dplyr::select(c(category, condition, geneRation, zscore))
  })
pancancer.circ <- do.call(rbind.fill, pancancer.circ)

(pancancer.circ_plot <- make_circPlot(pancancer.circ))
ggsave(file.path(folder, date, plots_dir, "pancancer_circplot.png"), 
       plot = pancancer.circ_plot, width = 17, height = 10, units = 'in')



