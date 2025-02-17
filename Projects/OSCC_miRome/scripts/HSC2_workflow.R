####################################################
# 1.) Set up working directory and directory tree  #
####################################################
# Create the sub folders for: results, data, and pictures
#data folder
# Set downstream path
pw <- getwd()
folder <- "mRNA"
# Create the dated results folder
data_dir <- "data"
if (!dir.exists(file.path(data_dir, folder))) {
  dir.create(file.path(data_dir, folder)) # create the main results folder
}
results_dir <- "results"
if (!dir.exists(file.path(results_dir, folder))) {
  dir.create(file.path(results_dir, folder)) # create the main results folder
}
# Copy readcounts from their source folder to the data folder
mrna.path <- choose.dir(getwd(), "Select the directory containing the count files")
mrna.files <- list.files(mrna.path, pattern = "counts.txt$", full.names = T)
file.copy(from = mrna.files, to = file.path(data_dir, folder))

# get the current date
date <- format(Sys.Date(), "%Y-%m-%d")
#plots directory
plots_dir <- "plots"
#tables directory
tables_dir <- "tables"
if (!dir.exists(file.path(results_dir, folder, date))) {
  dir.create(file.path(results_dir, folder, date), recursive = T) # create the dated results folder
  dir.create(file.path(results_dir, folder, date, tables_dir)) # create the tables folder
  dir.create(file.path(results_dir, folder, date, plots_dir)) # create the plots folder
}

#######################################
# 2.) Load functions for the analyses #
#######################################
source("scripts/functions.R")
source("scripts/packages.R")

########################## 
# 3.) Load analysis data #
##########################
if (exists("readcounts") == F) {
  mrna.reads <- lapply(list.files(file.path(data_dir, folder), 
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
  file.names <- list.files(file.path(data_dir, folder),
                           pattern = "counts.txt$", full.names = F)
  
  mrna.counts <- mrna.counts %>%
    # set the gene ID as row names
    tibble::column_to_rownames("Geneid") %>% 
    # clean up the column names
    setNames(gsub(".counts.txt", "", file.names))
  
  write.table(mrna.counts, paste(data_dir,folder,"readcounts.csv", sep = "/"), 
              sep =",", na = "NA", dec = ".", row.names = T, col.names = T)
}
# clear up workplace
rm(list = ls()[grepl("mrna.", ls())])

##########################
# 4.) Ready count tables #
##########################
HSC2.readcounts <- read.csv(file = file.path(pw, data_dir,folder,"readcounts.csv"),
                       sep = ",", header = T, na.strings = NA, row.names = 1) %>% 
  setNames(gsub(".counts.txt", "", colnames(.)))
HSC2.readcounts <- as.matrix(HSC2.readcounts) 

# split the readcounts into 1h and 6h
HSC2.readcounts <- list("1h" = HSC2.readcounts[,1:6], 
                        "6h" = HSC2.readcounts[,7:18])

##############################
# 5.) Prepare metadata table #
##############################

# Create the metadata table
HSC2.coldata <- data.frame("samples" = gsub(".counts.txt", "", file.names)) %>%
  dplyr::mutate(time = dplyr::case_when(
    stringr::str_detect(file.names, "H1") ~ "1h",
    stringr::str_detect(file.names, "H6") ~ "6h"),
    time = factor(time, 
                  levels = c("1h","6h"))) %>% 
  dplyr::group_split(time, .keep = F)

HSC2.coldata <- lapply(HSC2.coldata, function(x){
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

names(HSC2.coldata) <- c("1h","6h")

##############################################
# 6.) Make differential analyses with DESeq2 #
##############################################

# Run DESeq2 on the 1 hour infection model data
DESeq.1h <- calc_DiffExp(matrix = HSC2.readcounts[["1h"]], 
                         coldata = HSC2.coldata[["1h"]],
                         design = "condition")
# Run DESeq2 on the 6 hours infection model data
DESeq.6h <- calc_DiffExp(matrix = HSC2.readcounts[["6h"]],
                         coldata = HSC2.coldata[["6h"]],
                         design = "condition")

# PCA analysis of the 1 hour infection model data 
(pca.1h <- make_pca(DESeq.1h$dds_norm, "samples", HSC2.coldata[["1h"]]$condition))
ggsave(file.path(results_dir, folder, date, plots_dir, "1h_pca.png"), plot = pca.1h,
       width = 8, height = 8, units = 'in')

# PCA analysis of the 6 hours infection model data
(pca.6h <- make_pca(DESeq.6h$dds_norm, "samples", HSC2.coldata[["6h"]]$condition))
ggsave(file.path(results_dir, folder, date, plots_dir, "6h_pca.png"), plot = pca.6h,
       width = 8, height = 8, units = 'in')

# Extract significant results for the 1 hour infection model
Ca11_1h.res <- get_results(DESeq.1h$dds,
                           contrast = list("condition_CA11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1.5,
                           pval_treshold = 0.05)
# Check the number of DEGs:
table(Ca11_1h.res$sig_df$significance) # 9 down-regulated, 4 up-regulated

# Save results to excel files
sapply(names(Ca11_1h.res), function(x){
  openxlsx::write.xlsx(Ca11_1h.res[[x]], 
                       file.path(results_dir, folder, date, tables_dir, 
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
# Check the number of DEGs:
table(Ca11_6h.res$sig_df$significance) # 770 down-regulated, 480 up-regulated

# Save results to excel files
sapply(names(Ca11_6h.res), function(x){
  openxlsx::write.xlsx(Ca11_6h.res[[x]], 
                       file.path(results_dir, folder, date, tables_dir, 
                                 paste0("Ca11_6h_",x,".xlsx")),
                       rowNames = T)
})
## 2. Candida parapsilosis - MOI 1:1
Cp11_6h.res <- get_results(DESeq.6h$dds,
                           contrast = list("condition_CP11_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1.5,
                           pval_treshold = 0.05)
# Check the number of DEGs:
table(Cp11_6h.res$sig_df$significance) # 12 down-regulated, 26 up-regulated

# Save results to excel files
sapply(names(Cp11_6h.res), function(x){
  openxlsx::write.xlsx(Cp11_6h.res[[x]], 
                       file.path(results_dir, folder, date, tables_dir, 
                                 paste0("Cp11_6h_",x,".xlsx")),
                       rowNames = T)
})

## 3. Candida parapsilosis - MOI 5:1
Cp51_6h.res <- get_results(DESeq.6h$dds,
                           contrast = list("condition_CP51_vs_ctrl"),
                           name = c(''),
                           lfc_treshold = 1.5,
                           pval_treshold = 0.05)
# Check the number of DEGs:
table(Cp51_6h.res$sig_df$significance) # 25 down-regulated, 66 up-regulated

sapply(names(Cp51_6h.res), function(x){
  openxlsx::write.xlsx(Cp51_6h.res[[x]], 
                       file.path(results_dir, folder, date, tables_dir, 
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
ggsave(file.path(results_dir, folder, date, plots_dir, "Ca11_1h_vulcan.png"), 
       plot = Ca11_1h.res$plot, width = 10, height = 8, units = 'in')

## 1.2 6 hours infection model
# CA11
(Ca11_6h.res$plot <- make_vulcanplot(Ca11_6h.res$df, Ca11_6h.res$sig_df))
# save the plot
ggsave(file.path(results_dir, folder, date, plots_dir, "Ca11_6h_vulcan.png"), 
       plot = Ca11_6h.res$plot, width = 10, height = 8, units = 'in')

# CP11
(Cp11_6h.res$plot <- make_vulcanplot(Cp11_6h.res$df, Cp11_6h.res$sig_df))
# save the plot
ggsave(file.path(results_dir, folder, date, plots_dir, "Cp11_6h_vulcan.png"), 
       plot = Cp11_6h.res$plot, width = 10, height = 8, units = 'in')

# CP51
(Cp51_6h.res$plot <- make_vulcanplot(Cp51_6h.res$df, Cp51_6h.res$sig_df))
# save the plot
ggsave(file.path(results_dir, folder, date, plots_dir, "Cp51_6h_vulcan.png"), 
       plot = Cp51_6h.res$plot, width = 10, height = 8, units = 'in')

# 2. Heatmaps
## 2.1 1 hour infection model
(heatmap.1h <- make_heatmap(DESeq.1h$dds, "condition", "condition_CA11_vs_ctrl"))
ggsave(file.path(results_dir, folder, date, plots_dir, "1h_heatmap.png"), 
       plot = heatmap.1h, width = 8, height = 10, units = 'in')

## 2.2 6 hours infection model
(heatmap.6h <- make_heatmap(DESeq.6h$dds, "condition", "condition_CA11_vs_ctrl"))
ggsave(file.path(results_dir, folder, date, plots_dir, "6h_heatmap.png"), 
       plot = heatmap.6h, width = 8, height = 10, units = 'in')

##############################################
# 8.) Overrepresentation analyses            #
##############################################
library(msigdbr)
msigdbr_df <- msigdbr(species = "Homo sapiens")

# create reference databases
kegg_pathways <- msigdbr_df %>%
  dplyr::filter(
    gs_cat == "C2", # only canonical representations (compiled by experts)
    gs_subcat == "CP:KEGG" # KEGG pathways
  )

go_terms <- msigdbr_df %>%
  dplyr::filter(gs_cat == "C5" &  # ontology gene sets
                  gs_subcat != "HPO") # excluding human ontology phenotype gs

hallmark_sets <- msigdbr_df %>%
  dplyr::filter(gs_cat == "H") %>%  # only hallmark gene sets 
  dplyr::mutate(gs_subcat = "MSIGDB",
                gs_exact_source = gs_id)

# KEGG pathways
## 1 Get entrez IDs: gene set of interest and background
Ca11_1h.entrez <- get_genelist(.df = Ca11_1h.res$df,
                               .filter = Ca11_1h.res$df[["significance"]] %in% c('Signif. up-regulated',
                                                                                'Signif. down-regulated'),
                               .value = "stat",
                               .name = "entrezID") # CA11 1 hour
Ca11_6h.entrez <- get_genelist(.df = Ca11_6h.res$df,
                               .filter = Ca11_6h.res$df[["significance"]] %in% c('Signif. up-regulated',
                                                                                'Signif. down-regulated'),
                               .value = "stat",
                               .name = "entrezID") # CA11 6 hours
Cp11_6h.entrez <- get_genelist(.df = Cp11_6h.res$df,
                               .filter = Cp11_6h.res$df[["significance"]] %in% c('Signif. up-regulated',
                                                                                'Signif. down-regulated'),
                               .value = "stat",
                               .name = "entrezID") # CP11 6 hours
Cp51_6h.entrez <- get_genelist(.df = Cp51_6h.res$df,
                               .filter = Cp51_6h.res$df[["significance"]] %in% c('Signif. up-regulated',
                                                                                'Signif. down-regulated'),
                               .value = "stat",
                               .name = "entrezID") # CP51 6 hours

## 2 Run the ORA analysis
Ca11_1h.ORA.KEGG <- list()
Ca11_1h.ORA.KEGG <- run_ora(.interest = Ca11_1h.entrez$interest,
                        .background = Ca11_1h.entrez$background,
                        .pathways = kegg_pathways)
Ca11_1h.ORA.KEGG <- c(Ca11_1h.ORA.KEGG,
                      extract_ora_results(.ora = Ca11_1h.ORA.KEGG$ora,
                                          .db = kegg_pathways))
# Save the results
openxlsx::write.xlsx(Ca11_1h.ORA.KEGG[c(2:3)], 
                     file.path(results_dir, folder, date, tables_dir, "Ca11_1h_ORA.KEGG.xlsx"))

Ca11_6h.ORA.KEGG <- list()
Ca11_6h.ORA.KEGG <- run_ora(.interest = Ca11_6h.entrez$interest,
                        .background = Ca11_6h.entrez$background,
                        .pathways = kegg_pathways)
Ca11_6h.ORA.KEGG <- c(Ca11_6h.ORA.KEGG, extract_ora_results(.ora = Ca11_6h.ORA.KEGG$ora, .db = kegg_pathways))
# Save the results
openxlsx::write.xlsx(Ca11_6h.ORA.KEGG[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Ca11_6h_ORA.KEGG.xlsx"))

Cp11_6h.ORA.KEGG <- list()
Cp11_6h.ORA.KEGG <- run_ora(.interest = Cp11_6h.entrez$interest,
                         .background = Cp11_6h.entrez$background,
                         .pathways = kegg_pathways)
Cp11_6h.ORA.KEGG <- c(Cp11_6h.ORA.KEGG, 
                      extract_ora_results(.ora = Cp11_6h.ORA.KEGG$ora, .db = kegg_pathways))
# Save the results
openxlsx::write.xlsx(Cp11_6h.ORA.KEGG[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Cp11_6h_ORA.KEGG.xlsx"))

Cp51_6h.ORA.KEGG <- list()
Cp51_6h.ORA.KEGG <- run_ora(.interest = Cp51_6h.entrez$interest,
                            .background = Cp51_6h.entrez$background,
                            .pathways = kegg_pathways)
Cp51_6h.ORA.KEGG <- c(Cp51_6h.ORA.KEGG, 
                      extract_ora_results(.ora = Cp51_6h.ORA.KEGG$ora, .db = kegg_pathways))
# Save the results
openxlsx::write.xlsx(Cp51_6h.ORA.KEGG[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Cp51_6h_ORA.KEGG.xlsx"))

## 3 Graphical analysis
(Ca11_1h.ORA.KEGG$plot <- make_dotplot(Ca11_1h.ORA.KEGG$sig_df, 
                                       Ca11_1h.ORA.KEGG$sig_df$Count, type = "KEGG"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Ca11_1h_ORA.KEGG.png"), 
       plot = Ca11_1h.ORA.KEGG$plot, width = 10, height = 8, units = 'in')

(Ca11_6h.ORA.KEGG$plot <- make_dotplot(Ca11_6h.ORA.KEGG$sig_df, 
                                       Ca11_6h.ORA.KEGG$sig_df$Count, type = "KEGG"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Ca11_6h_ORA.KEGG.png"), 
       plot = Ca11_6h.ORA.KEGG$plot, width = 10, height = 8, units = 'in')

(Cp11_6h.ORA.KEGG$plot <- make_dotplot(Cp11_6h.ORA.KEGG$sig_df, 
                                   Cp11_6h.ORA.KEGG$sig_df$Count, type = "KEGG"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Cp11_6h_ORA.KEGG.png"), 
      plot = Cp11_6h.ORA.KEGG$plot, width = 10, height = 8, units = 'in')

# no significant hits
#(Cp51_6h.ORA.KEGG$plot <- make_dotplot(Cp51_6h.ORA.KEGG$sig_df, 
#                                       Cp51_6h.ORA.KEGG$sig_df$Count, type = "KEGG"))
#ggsave(file.path(results_dir, folder, date, plots_dir,"Cp51_6h_KEGG.png"), 
#       plot = Cp51_6h.KEGG$plot, width = 10, height = 8, units = 'in')


##############################################
# 9.) Gene set enrichment analyses           #
##############################################
# ----------------------------- KEGG pathways -------------------------------- #
# Run the GSEA analysis
Ca11_1h.KEGG <- list()
Ca11_1h.KEGG <- run_gsea(
  .geneset = Ca11_1h.entrez$background[!is.na(names(Ca11_1h.entrez$background))], 
  .terms = kegg_pathways)
Ca11_1h.KEGG <- c(Ca11_1h.KEGG, 
                  extract_gsea_results(.gsea = Ca11_1h.KEGG$gsea, .db = kegg_pathways))
# Save the results
openxlsx::write.xlsx(Ca11_1h.KEGG[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Ca11_1h_GSEA.KEGG.xlsx"))

Ca11_6h.KEGG <- list()
Ca11_6h.KEGG <- run_gsea(
  .geneset = Ca11_6h.entrez$background[!is.na(names(Ca11_6h.entrez$background))], 
  .terms = kegg_pathways)

Ca11_6h.KEGG <- c(Ca11_6h.KEGG, 
                  extract_gsea_results(.gsea = Ca11_6h.KEGG$gsea, .db = kegg_pathways))

# Save the results
openxlsx::write.xlsx(Ca11_6h.KEGG[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Ca11_6h_GSEA.KEGG.xlsx"))

Cp11_6h.KEGG <- list()
Cp11_6h.KEGG <- run_gsea(
  .geneset = Cp11_6h.entrez$background[!is.na(names(Cp11_6h.entrez$background))], 
  .terms = kegg_pathways)

Cp11_6h.KEGG <- c(Cp11_6h.KEGG, 
                  extract_gsea_results(.gsea = Cp11_6h.KEGG$gsea, .db = kegg_pathways))

# Save the results
openxlsx::write.xlsx(Cp11_6h.KEGG[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Cp11_6h_GSEA.KEGG.xlsx"))

Cp51_6h.KEGG <- list()
Cp51_6h.KEGG <- run_gsea(
  .geneset = Cp51_6h.entrez$background[!is.na(names(Cp51_6h.entrez$background))], 
  .terms = kegg_pathways)

Cp51_6h.KEGG <- c(Cp51_6h.KEGG, 
                  extract_gsea_results(.gsea = Cp51_6h.KEGG$gsea, .db = kegg_pathways))

# Save the results
openxlsx::write.xlsx(Cp51_6h.KEGG[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Cp51_6h_GSEA.KEGG.xlsx"))

# Graphical analysis
## CA11 - 1h - no significant hits

## CA11 - 6h
(Ca11_6h.KEGG$plot <- make_gseaplot(data =Ca11_6h.KEGG$sig_df, type = "Hallmark"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Ca11_6h_GSEA.KEGG.png"), 
       plot = Ca11_6h.KEGG$plot, width = 10, height = 8, units = 'in')

## CP11 - 6h
(Cp11_6h.KEGG$plot <- make_gseaplot(data = Cp11_6h.KEGG$sig_df, type = "Hallmark"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Cp11_6h_GSEA.KEGG.png"), 
       plot = Cp11_6h.KEGG$plot, width = 10, height = 8, units = 'in')

## CP51 - 6h
(Cp51_6h.KEGG$plot <- make_gseaplot(data = Cp51_6h.KEGG$sig_df, type = "Hallmark"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Cp51_6h_GSEA.KEGG.png"), 
       plot = Cp51_6h.KEGG$plot, width = 10, height = 8, units = 'in')

# ----------------------------- GO term -------------------------------- #
# Run the GSEA analysis
Ca11_1h.GO <- list()
Ca11_1h.GO <- run_gsea(.geneset = Ca11_1h.entrez$background[!is.na(names(Ca11_1h.entrez$background))], 
                       .terms = go_terms)
Ca11_1h.GO <- c(Ca11_1h.GO, extract_gsea_results(.gsea = Ca11_1h.GO$gsea, .db = go_terms))
# Save the results
openxlsx::write.xlsx(Ca11_1h.GO[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Ca11_1h_GO.xlsx"))

Ca11_6h.GO <- list()
Ca11_6h.GO <- run_gsea(.geneset = Ca11_6h.entrez$background[!is.na(names(Ca11_6h.entrez$background))], 
                       .terms = go_terms)
Ca11_6h.GO <- c(Ca11_6h.GO, extract_gsea_results(.gsea = Ca11_6h.GO$gsea, .db = go_terms))
# Save the results
openxlsx::write.xlsx(Ca11_6h.GO[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Ca11_6h_GO.xlsx"))

Cp11_6h.GO <- list()
Cp11_6h.GO <- run_gsea(.geneset = Cp11_6h.entrez$background[!is.na(names(Cp11_6h.entrez$background))], 
                          .terms = go_terms)
Cp11_6h.GO <- c(Cp11_6h.GO, extract_gsea_results(.gsea = Cp11_6h.GO$gsea, .db = go_terms))
# Save the results
openxlsx::write.xlsx(Cp11_6h.GO[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Cp11_6h_GO.xlsx"))

Cp51_6h.GO <- list()
Cp51_6h.GO <- run_gsea(.geneset = Cp51_6h.entrez$background[!is.na(names(Cp51_6h.entrez$background))], 
                       .terms = go_terms)
Cp51_6h.GO <- c(Cp51_6h.GO, extract_gsea_results(.gsea = Cp51_6h.GO$gsea, .db = go_terms))
# Save the results
openxlsx::write.xlsx(Cp51_6h.GO[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Cp51_6h_GO.xlsx"))

# 1.2 Graphical analysis
## CA11 - 1h 
(Ca11_1h.GO$plot <- make_gseaplot(data = Ca11_1h.GO$sig_df, type = "GO"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Ca11_1h_GO.png"), 
       plot = Ca11_1h.GO$plot, width = 10, height = 8, units = 'in')
## CA11 - 6h
(Ca11_6h.GO$plot <- make_gseaplot(data = Ca11_6h.GO$sig_df, type = "GO"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Ca11_6h_GO.png"), 
       plot = Ca11_6h.GO$plot, width = 10, height = 8, units = 'in')
# CP11 - 6h
(Cp11_6h.GO$plot <- make_gseaplot(Cp11_6h.GO$sig_df, type = "GO"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Cp11_6h_GO.png"), 
       plot = Cp11_6h.GO$plot, width = 10, height = 8, units = 'in')
# CP51 - 6h
(Cp51_6h.GO$plot <- make_gseaplot(Cp51_6h.GO$sig_df, type = "GO"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Cp51_6h_GO.png"), 
       plot = Cp51_6h.GO$plot, width = 10, height = 8, units = 'in')

# ----------------------------- MSigDB Hallmark ------------------------------ #
#  Run the GSEA analysis
## Ca11 - 1h 
Ca11_1h.Hallmark <- list()
Ca11_1h.Hallmark <- run_gsea(
  .geneset = Ca11_1h.entrez$background[!is.na(names(Ca11_1h.entrez$background))], 
  .terms = hallmark_sets)
Ca11_1h.Hallmark <- c(Ca11_1h.Hallmark, 
                      extract_gsea_results(.gsea = Ca11_1h.Hallmark$gsea, 
                                           .db = hallmark_sets))
# Save the results
openxlsx::write.xlsx(Ca11_1h.Hallmark[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Ca11_1h_Hallmark.xlsx"))
## Ca11 - 6h
Ca11_6h.Hallmark <- list()
Ca11_6h.Hallmark <- run_gsea(
  .geneset = Ca11_6h.entrez$background[!is.na(names(Ca11_6h.entrez$background))], 
  .terms = hallmark_sets)
Ca11_6h.Hallmark <- c(Ca11_6h.Hallmark, 
                      extract_gsea_results(.gsea = Ca11_6h.Hallmark$gsea, 
                                           .db = hallmark_sets))
# Save the results
openxlsx::write.xlsx(Ca11_6h.Hallmark[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Ca11_6h_Hallmark.xlsx"))

## CP11 - 6h
Cp11_6h.Hallmark <- list()
Cp11_6h.Hallmark <- run_gsea(
  .geneset = Cp11_6h.entrez$background[!is.na(names(Cp11_6h.entrez$background))], 
  .terms = hallmark_sets)
Cp11_6h.Hallmark <- c(Cp11_6h.Hallmark, 
                      extract_gsea_results(.gsea = Cp11_6h.Hallmark$gsea, 
                                           .db = hallmark_sets))
# Save the results
openxlsx::write.xlsx(Cp11_6h.Hallmark[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Cp11_6h_Hallmark.xlsx"))

## CP51 - 6h
Cp51_6h.Hallmark <- list()
Cp51_6h.Hallmark <- run_gsea(
  .geneset = Cp51_6h.entrez$background[!is.na(names(Cp51_6h.entrez$background))], 
  .terms = hallmark_sets)
Cp51_6h.Hallmark <- c(Cp51_6h.Hallmark,
                      extract_gsea_results(.gsea = Cp51_6h.Hallmark$gsea,
                                           .db = hallmark_sets))
# Save the results
openxlsx::write.xlsx(Cp51_6h.Hallmark[c(2:3)],
                     file.path(results_dir, folder, date, tables_dir, "Cp51_6h_Hallmark.xlsx"))

## 2.2 Graphical analysis
## CA11 - 1h
(Ca11_1h.Hallmark$plot <- make_gseaplot(data = Ca11_1h.Hallmark$sig_df, type = "Hallmark"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Ca11_1h_Hallmark.png"), 
       plot = Ca11_1h.Hallmark$plot, width = 10, height = 8, units = 'in')
## CA11 - 6h
(Ca11_6h.Hallmark$plot <- make_gseaplot(data = Ca11_6h.Hallmark$sig_df, type = "Hallmark"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Ca11_6h_Hallmark.png"), 
       plot = Ca11_6h.Hallmark$plot, width = 10, height = 8, units = 'in')
## CP11 - 6h
(Cp11_6h.Hallmark$plot <- make_gseaplot(data = Cp11_6h.Hallmark$sig_df, type = "Hallmark"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Cp11_6h_Hallmark.png"), 
       plot = Cp11_6h.Hallmark$plot, width = 10, height = 8, units = 'in')
## CP51 - 6h
(Cp51_6h.Hallmark$plot <- make_gseaplot(data = Cp51_6h.Hallmark$sig_df, type = "Hallmark"))
ggsave(file.path(results_dir, folder, date, plots_dir,"Cp51_6h_Hallmark.png"), 
       plot = Cp51_6h.Hallmark$plot, width = 10, height = 8, units = 'in')

##########################################################
# 11.) Network analysis of enriched genesets (GSEA)      #
##########################################################
# Calculate Cohen's similarity between all the pathways, terms and hallmark#
if (!file.exists("./data/term_similarity_matrix.RData")){
  # extract a list of gene sets from every reference used:
  total.genesets <- c(split(go_terms$gene_symbol, # GO terms
                            go_terms$gs_exact_source),
                      split(kegg_pathways$gene_symbol, # KEGG and Reactome pathways
                            kegg_pathways$gs_exact_source),
                      split(hallmark_sets$gene_symbol, # Hallmark gene sets
                            hallmark_sets$gs_exact_source))
  # calculate the similarity matrix
  similarity.matrix <- cohen_kappa(total.genesets)
  # save the similarity matrix as RData
  save(similarity.matrix, file = "./data/term_similarity_matrix.RData")
  # save file
  write.xlsx(similarity.matrix, file = "./data/term_similarity_matrix.xlsx")
} else {
  load("./data/term_similarity_matrix.RData")
}

# combine GO and pathway enrichment results
### CA11 - 1h
Ca11_1h.combinedGSEA <- rbind(Ca11_1h.GO$sig_df, # GO terms
                              Ca11_1h.KEGG$sig_df, # KEGG pathways
                              Ca11_1h.Hallmark$sig_df) %>% # Hallmark gene sets
  dplyr::rowwise(.) %>% 
  # Calculate background ratio
  dplyr::mutate(
    geneRatio = length(unlist(strsplit(core_enrichment,"\\/")))/setSize
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(c("ID", "Name", "setSize", "geneRatio", "NES", "p.adjust", "core_enrichment"),
                  .before = everything())
# Cluster enriched terms on gene set similarity
Ca11_1h.cluster <- get_cluster(Ca11_1h.combinedGSEA, similarity.matrix, .threshold = 0.25)
# Select cluster representative terms based on the involved genes' importance
Ca11_1h.cluster$df <- get_cluster_representative(.cluster = Ca11_1h.cluster$df,
                                                 .degs = Ca11_1h.res$df)
V(Ca11_1h.cluster$graph)$Representative <- Ca11_1h.cluster$df$Representative
V(Ca11_1h.cluster$graph)$Description <- Ca11_1h.cluster$df$Name
# Save results
write.xlsx(Ca11_1h.cluster$df, 
           file.path(results_dir, folder, date, tables_dir, "Ca11_1h_total_GSEA_clusters.xlsx"))

### CA11 - 6h
Ca11_6h.combinedGSEA <- rbind(filter(Ca11_6h.GO$sig_df, p.adjust < 0.001), # GO terms
                              filter(Ca11_6h.KEGG$sig_df, p.adjust < 0.001),# KEGG pathways
                              filter(Ca11_6h.Hallmark$sig_df, p.adjust < 0.001)) %>% # Hallmark gene sets
  # Hallmark gene sets
  dplyr::rowwise(.) %>% 
  # Calculate background ratio
  dplyr::mutate(
    geneRatio = length(unlist(strsplit(core_enrichment,"\\/")))/setSize
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(c("ID", "Name", "setSize", "geneRatio", "NES", "p.adjust", "core_enrichment"),
                  .before = everything())

# Cluster enriched terms on gene set similarity
Ca11_6h.cluster <- get_cluster(Ca11_6h.combinedGSEA, similarity.matrix, .threshold = 0.25)
# Select cluster representative terms based on the involved genes' importance
Ca11_6h.cluster$df <- get_cluster_representative(.cluster = Ca11_6h.cluster$df,
                                                 .degs = Ca11_6h.res$df)
V(Ca11_6h.cluster$graph)$Representative <- Ca11_6h.cluster$df$Representative
V(Ca11_6h.cluster$graph)$Description <- Ca11_6h.cluster$df$Name
# Save results
write.xlsx(Ca11_6h.cluster$df, 
           file.path(results_dir, folder, date, tables_dir, "Ca11_6h_total_GSEA_clusters.xlsx"))

### CP11 - 6h
Cp11_6h.combinedGSEA <- rbind(Cp11_6h.GO$sig_df, # GO terms
                              Cp11_6h.KEGG$sig_df, # KEGG pathways
                              Cp11_6h.Hallmark$sig_df) %>% # Hallmark gene sets
  dplyr::rowwise(.) %>% 
  # Calculate background ratio
  dplyr::mutate(
    geneRatio = length(unlist(strsplit(core_enrichment,"\\/")))/setSize
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(c("ID", "Name", "setSize", "geneRatio", "NES", "p.adjust", "core_enrichment"),
                  .before = everything())

# Cluster enriched terms on gene set similarity
Cp11_6h.cluster <- get_cluster(Cp11_6h.combinedGSEA, similarity.matrix, .threshold = 0.25)
# Select cluster representative terms based on the involved genes' importance
Cp11_6h.cluster$df <- get_cluster_representative(.cluster = Cp11_6h.cluster$df,
                                                 .degs = Cp11_6h.res$df)
V(Cp11_6h.cluster$graph)$Representative <- Cp11_6h.cluster$df$Representative
V(Cp11_6h.cluster$graph)$Description <- Cp11_6h.cluster$df$Name
# Save results
write.xlsx(Cp11_6h.cluster$df, 
           file.path(results_dir, folder, date, tables_dir, "Cp11_6h_total_GSEA_clusters.xlsx"))

### CP51 - 6h
Cp51_6h.combinedGSEA <- rbind(Cp51_6h.GO$sig_df, # GO terms
                              Cp51_6h.KEGG$sig_df, # KEGG pathways
                              Cp51_6h.Hallmark$sig_df) %>% # Hallmark gene sets
  dplyr::rowwise(.) %>% 
  # Calculate background ratio
  dplyr::mutate(
    geneRatio = length(unlist(strsplit(core_enrichment,"\\/")))/setSize
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::relocate(c("ID", "Name", "setSize", "geneRatio", "NES", "p.adjust", "core_enrichment"),
                  .before = everything())

# Cluster enriched terms on gene set similarity
Cp51_6h.cluster <- get_cluster(Cp51_6h.combinedGSEA, similarity.matrix, .threshold = 0.25)
# Select cluster representative terms based on the involved genes' importance
Cp51_6h.cluster$df <- get_cluster_representative(.cluster = Cp51_6h.cluster$df,
                                                 .degs = Cp51_6h.res$df)
V(Cp51_6h.cluster$graph)$Representative <- Cp51_6h.cluster$df$Representative
V(Cp51_6h.cluster$graph)$Description <- Cp51_6h.cluster$df$Name
# Save results
write.xlsx(Cp51_6h.cluster$df, 
           file.path(results_dir, folder, date, tables_dir, "Cp51_6h_total_GSEA_clusters.xlsx"))

# Network visualization of the enriched clusters 
# CA11 - 1h
Ca11_1h.cluster$sub_graph <- filter_graph(Ca11_1h.cluster$graph, 5)
set.seed(42)
Ca11_1h.cluster$layout <- layout_with_fr(Ca11_1h.cluster$sub_graph)

(Ca11_1h.cluster$plot <- plot_network(.net = Ca11_1h.cluster$sub_graph,
                                      .layout = Ca11_1h.cluster$layout,
                                      .labels =  V(Ca11_1h.cluster$sub_graph)$Name,
                                      .df = Ca11_1h.cluster$df))
# Save plot
ggsave(file.path(results_dir, folder, date, plots_dir, "Ca11_1h_GSEA_network.png"),
       plot = Ca11_1h.cluster$plot, bg = "white",
       width = 20, height = 14, units = "in")

# CA11 - 6h
Ca11_6h.cluster$sub_graph <- filter_graph(Ca11_6h.cluster$graph, 5)
set.seed(42)
Ca11_6h.cluster$layout <- layout_with_fr(Ca11_6h.cluster$sub_graph)

(Ca11_6h.cluster$plot <- plot_network(.net = Ca11_6h.cluster$sub_graph,
                                      .layout = Ca11_6h.cluster$layout,
                                      .labels =  V(Ca11_6h.cluster$sub_graph)$Name,
                                      .df = Ca11_6h.cluster$df))

# Save plot
ggsave(file.path(results_dir, folder, date, plots_dir, "Ca11_6h_GSEA_network.png"),
       plot = Ca11_6h.cluster$plot, bg = "white",
       width = 20, height = 14, units = "in")

# CP11 - 6h
Cp11_6h.cluster$sub_graph <- filter_graph(Cp11_6h.cluster$graph, 5)
set.seed(42)
Cp11_6h.cluster$layout <- layout_with_fr(Cp11_6h.cluster$sub_graph)

(Cp11_6h.cluster$plot <- plot_network(.net = Cp11_6h.cluster$sub_graph,
                                      .layout = Cp11_6h.cluster$layout,
                                      .labels =  V(Cp11_6h.cluster$sub_graph)$Name,
                                      .df = Cp11_6h.cluster$df))

# Save plot
ggsave(file.path(results_dir, folder, date, plots_dir, "Cp11_6h_GSEA_network.png"),
       plot = Cp11_6h.cluster$plot, bg = "white",
       width = 20, height = 14, units = "in")

# CP51 - 6h
Cp51_6h.cluster$sub_graph <- filter_graph(Cp51_6h.cluster$graph, 5)
set.seed(42)
Cp51_6h.cluster$layout <- layout_with_fr(Cp51_6h.cluster$sub_graph)

(Cp51_6h.cluster$plot <- plot_network(.net = Cp51_6h.cluster$sub_graph,
                                      .layout = Cp51_6h.cluster$layout,
                                      .labels =  V(Cp51_6h.cluster$sub_graph)$Name,
                                      .df = Cp51_6h.cluster$df))

# Save plot
ggsave(file.path(results_dir, folder, date, plots_dir, "Cp51_6h_GSEA_network.png"),
       plot = Cp51_6h.cluster$plot, bg = "white",
       width = 20, height = 14, units = "in")

##############################################
# 11.) Overlap with cancer hallmark genesets #
##############################################
if (exists("pancancer.genes") == F) {
  pancancer.path <- choose.files(getwd(), "Select the directory containing the count files",
                                 multi = F)
  file.copy(from = pancancer.path, to = file.path(folder, data_dir),
            overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
}

pancancer.file <- list.files(file.path(data_dir, folder), pattern = "cancer",
                             full.names = F)

pancancer.df <- require_file(file.path(data_dir, folder, pancancer.file),
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
  ggsave(file.path(results_dir, folder, date, plots_dir, paste0("HSC2_vs_",x,".png")),
         plot = pancancer.venn[[x]], width = 16, height = 10, units = 'in')
})

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
ggsave(file.path(results_dir, folder, date, plots_dir, "pancancer_circplot.png"), 
       plot = pancancer.circ_plot, width = 17, height = 10, units = 'in')




