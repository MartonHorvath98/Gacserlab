####################################################
# 1.) Set up working directory and directory tree  #
####################################################
.libPaths("C:/Program Files/R/R-4.2.1/library")
setwd("D:/Evi/OSCC_markers/RStudio")

  #data folder
  data_dir <- "data"
  if (!dir.exists("data")) {
    dir.create("data")
  }
  
  #plots directory
  plots_dir <- "plots"
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
  }
  
  #results directory
  results_dir <- "results" 
  if (!dir.exists(results_dir)) {
    dir.create(results_dir)
  }
####################################################
# 2.) Load packages and functions for the analyses #
####################################################
source("packages.R")
source("functions.R")

########################## 
# 3.) Load analysis data #
##########################
# copy count data from Project #1 - ECE1 KO project: 
# C. albicans BWP17 and SC5314 strains
file.copy(from="D:/Evi/ECE1_KO_project/live_CA/RStudio/ECE_project_readcounts.csv", to=data_dir, 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)
# add file path to local variable
ECE_file <- file.path(data_dir,
                      "ECE_project_readcounts.csv")

# copy count data from Project #2 - Mate project:
# C. albicans SC5314 strain
file.copy(from="D:/Mate/Mate_in_vitro/RStudio/Mate_readcounts.csv", to=data_dir, 
          overwrite = TRUE, recursive = FALSE, 
          copy.mode = TRUE)
# add file path to local variable
MATE_file <- file.path(data_dir,
                      "Mate_readcounts.csv")

##########################
# 4.) Ready count tables #
##########################
library(dplyr)

# load count data
ECE.readcounts <- require_file(ECE_file, header = T, na.strings = NA)

ECE.readcounts <- ECE.readcounts %>%
  round() %>%
  data.frame() %>%
  dplyr::select(!contains("ECE1")) %>% # deselect ECE1 KO samples
  dplyr::filter(rowSums(.) >= 50) # filter genes with low counts

# load count data
Mate.readcounts <- require_file(MATE_file, header = T, 
                                row.names = 1, na.strings = NA)

Mate.readcounts <- Mate.readcounts %>%
  round() %>%
  data.frame() %>%
  dplyr::filter(rowSums(.) >= 50) # filter genes with low counts

##############################
# 5.) Prepare metadata table #
##############################
library(stringr)
library(dplyr)

ECE.names <- colnames(ECE.readcounts) # get sample names
ECE.coldata <- data.frame("samples" = ECE.names) %>%
  dplyr::mutate(samples = as.factor(samples),
         treatment = dplyr::case_when(
           stringr::str_detect(ECE.names, "_Bwp17_") ~ "BWP17",
           stringr::str_detect(ECE.names, "_SC5314_") ~ "SC5314",
           T ~ "ctrl"
         ), # assign treatment to samples
         treatment = factor(treatment, 
                            levels = c("ctrl","BWP17","SC5314")))

Mate.names <- colnames(Mate.readcounts) # get sample names
Mate.coldata <- data.frame("samples" = Mate.names) %>%
  dplyr::mutate(samples = as.factor(samples),
         treatment = dplyr::case_when(
           stringr::str_detect(Mate.names, "_Ca_") ~ "SC5314",
           T ~ "ctrl"
         ), # assign treatment to samples
         treatment = factor(treatment, 
                            levels = c("ctrl","SC5314")))

##############################################
# 6.) Make differential analyses with DESeq2 #
##############################################
library(DESeq2)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
# run differential expression analysis
ECE.DE <- calc_DiffExp(matrix = ECE.readcounts, 
                       coldata = ECE.coldata,
                       design = "treatment") # design formula
# run differential expression analysis
Mate.DE <- calc_DiffExp(matrix = Mate.readcounts, 
                       coldata = Mate.coldata,
                       design = "treatment") # design formula

(ECE.pca <- make_pca(ECE.DE$dds_norm, ECE.coldata$treatment)) # make PCA plot
ggsave(paste(plots_dir,"/ECE_pca.png"), plot = ECE.pca,
       width = 6, height = 6, units = 'in') # save PCA plot

(Mate.pca <- make_pca(Mate.DE$dds_norm, Mate.coldata$treatment)) # make PCA plot
ggsave(paste(plots_dir,"/Mate_pca.png"), plot = Mate.pca,
       width = 6, height = 6, units = 'in') # save PCA plot

library(fdrtool)
library(ashr)
library(org.Hs.eg.db)
# calculate differential expression
BWP17.ECE.results <- get_results(ECE.DE$dds, # DESeq2 object
                           contrast = list("treatment_BWP17_vs_ctrl"), 
                           # contrast of interest
                           name = c('')) 
# select significant results
BWP17.ECE.df <- sig_results(BWP17.ECE.results$resLFC,
                            sig_log2FC = 1.5, # effect size threshold
                            sig_pval = 0.05) # p-value threshold
# make volcano plot
BWP17.ECE.df$plot <- make_vulcanplot(BWP17.ECE.df$df, BWP17.ECE.df$sig_df) 
ggsave(paste(plots_dir,"/vulcanplot_BWP17.ECE.png"), plot = BWP17.ECE.df$plot,
       width = 10, height = 8, units = 'in') # save volcano plot

sapply(names(BWP17.ECE.df), function(x){
  openxlsx::write.xlsx(BWP17.ECE.df[[x]], paste0(results_dir, "/BWP17.ECE.",x,".xlsx"), rowNames = T)
}) # save results to excel

# calculate differential expression
SC5314.ECE.results <- get_results(ECE.DE$dds, # DESeq2 object
                                 contrast = list("treatment_SC5314_vs_ctrl"), 
                                 # contrast of interest 
                                 name = c(''))
# select significant results
SC5314.ECE.df <- sig_results(SC5314.ECE.results$resLFC,
                            sig_log2FC = 1.5, # effect size threshold
                            sig_pval = 0.05) # p-value threshold
# make volcano plot
SC5314.ECE.df$plot <- make_vulcanplot(SC5314.ECE.df$df, SC5314.ECE.df$sig_df)
ggsave(paste(plots_dir,"/vulcanplot_BWP17.ECE.png"), plot = BWP17.ECE.df$plot,
       width = 10, height = 8, units = 'in') # save volcano plot

sapply(names(SC5314.ECE.df), function(x){
  openxlsx::write.xlsx(SC5314.ECE.df[[x]], paste0(results_dir, "/SC5314.ECE.",x,".xlsx"), rowNames = T)
}) # save results to excel

# calculate differential expression
SC5314.Mate.results <- get_results(Mate.DE$dds, # DESeq2 object
                                 contrast = list("treatment_SC5314_vs_ctrl"), 
                                 # contrast of interest
                                 name = c(''))
# select significant results
SC5314.Mate.df <- sig_results(SC5314.Mate.results$resLFC,
                            sig_log2FC = 1.5, # effect size threshold
                            sig_pval = 0.05) # p-value threshold
# make volcano plot
SC5314.Mate.df$plot <- make_vulcanplot(SC5314.Mate.df$df, SC5314.Mate.df$sig_df)
ggsave(paste(plots_dir,"/vulcanplot_BWP17.ECE.png"), plot = BWP17.ECE.df$plot,
       width = 10, height = 8, units = 'in') # save volcano plot

sapply(names(SC5314.Mate.df), function(x){
       openxlsx::write.xlsx(SC5314.Mate.df[[x]], paste0(results_dir, "/SC5314.Mate.",x,".xlsx"), rowNames = T)
  }) # save results to excel


####################################
# 7.) Over representation analyses #
####################################
library(clusterProfiler)
library(msigdbr) 
library(org.Hs.eg.db)

msigdbr_df <- msigdbr(species = "Homo sapiens") # save local database

kegg_pathways <- msigdbr_df %>% 
  dplyr::filter(
    gs_cat == "C2", # only canonical representations (compiled by experts)
    gs_subcat == "CP:KEGG" # KEGG pathways
  )

# pull a list of unique gene names from the total background and set of interest 
# (DEGs) and convert them to ENTREZ format
BWP17.ECE.entrez <- get_entrez(BWP17.ECE.df) 
SC5314.ECE.entrez <- get_entrez(SC5314.ECE.df)
SC5314.Mate.entrez <- get_entrez(SC5314.Mate.df)

# Get over-represented KEGG pathways using the enricher() function from the 
# clusterProfiler package
# parameters: pValueCutoff = 0.1, pAdjustMethod = "Benjamini-Hochberg"
BWP17.ECE.KEGGs <- kegg_results(genes_of_interest = BWP17.ECE.entrez$interest, 
                                background_set = BWP17.ECE.entrez$background,
                                kegg_df = kegg_pathways) 
SC5314.ECE.KEGGs <- kegg_results(genes_of_interest = SC5314.ECE.entrez$interest, 
                                background_set = SC5314.ECE.entrez$background,
                                kegg_df = kegg_pathways)
SC5314.Mate.KEGGs <- kegg_results(genes_of_interest = SC5314.Mate.entrez$interest, 
                                 background_set = SC5314.Mate.entrez$background,
                                 kegg_df = kegg_pathways)

library(ggplot2)
library(forcats)
library(numbers)
BWP17.ECE.KEGGs$plot <- make_dotplot(BWP17.ECE.KEGGs$df) # make KEGG dot plot
ggsave(paste(plots_dir,"/KEGG_BWP17.ECE.png"), plot = BWP17.ECE.KEGGs$plot,
       width = 10, height = 6, units = 'in') # save plot

SC5314.ECE.KEGGs$plot <- make_dotplot(SC5314.ECE.KEGGs$df) # make KEGG dot plot
ggsave(paste(plots_dir,"/KEGG_SC5314.ECE.png"), plot = SC5314.ECE.KEGGs$plot,
       width = 10, height = 6, units = 'in') # save plot

SC5314.Mate.KEGGs$plot <- make_dotplot(SC5314.Mate.KEGGs$df) # make KEGG dot plot
ggsave(paste(plots_dir,"/KEGG_SC5314.Mate.png"), plot = SC5314.Mate.KEGGs$plot,
       width = 10, height = 6, units = 'in') # save plot


####################################
# 8.) Gene set enrichment analyses #
####################################
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)


hallmark_sets <- msigdbr_df %>%
  dplyr::filter(gs_cat == "H") # only hallmark gene sets

# Using the GSEA() function from the clusterProfiler package, get gene set 
# enrichment analysis results for the hallmark gene sets
# parameters: minGSSize = 25, maxGSSize = 500, pvalueCutoff = 0.05, 
#             pAdjustMethod = "BH"
BWP17.ECE.gsea <- gsea_results(BWP17.ECE.df$df, hallmark_sets)
BWP17.ECE.gsea_plot <- make_gseaplot(BWP17.ECE.gsea$gsea_df) # make GSEA plot
ggsave(paste(plots_dir,"/GSEA_BWP17.ECE.png"), plot = BWP17.ECE.gsea_plot,
       width = 12, height = 6, units = 'in') # save plot

SC5314.ECE.gsea <- gsea_results(SC5314.ECE.df$df, hallmark_sets)
SC5314.ECE.gsea_plot <- make_gseaplot(SC5314.ECE.gsea$gsea_df) # make GSEA plot
ggsave(paste(plots_dir,"/GSEA_SC5314.ECE.png"), plot = SC5314.ECE.gsea_plot,
       width = 12, height = 6, units = 'in') # save plot

SC5314.Mate.gsea <- gsea_results(SC5314.Mate.df$df, hallmark_sets)
SC5314.Mate.gsea_plot <- make_gseaplot(SC5314.Mate.gsea$gsea_df) # make GSEA plot
ggsave(paste(plots_dir,"/GSEA_SC5314.Mate.png"), plot = SC5314.Mate.gsea_plot,
       width = 12, height = 6, units = 'in') # save plot

#################################
# 9.) Gene sets in interception #
#################################
library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)
library(tibble)
library(VennDiagram)
library(GOplot)
library(ComplexUpset)
library(org.Hs.eg.db)
library(conicfit)


GO_BP <- msigdbr_df %>%
  dplyr::filter(gs_cat == "C5" & gs_subcat == "GO:BP") # biological processes (GO)
GO_MF <- msigdbr_df %>%
  dplyr::filter(gs_cat == "C5" & gs_subcat == "GO:MF") # molecular functions (GO)

# Create Venn diagram of shared differentially expressed genes 
Calbi.venn <- make_venn(BWP17.ECE.df$sig_df, SC5314.ECE.df$sig_df, SC5314.Mate.df$sig_df,
                        c("BWP17", "SC5314.ECE", "SC5314.Mate"), "venn_diagram") 

# total background gene set
HSC.gene <- unique(c(as.character(BWP17.ECE.df$df$entrezID),
                     as.character(SC5314.ECE.df$df$entrezID),
                     as.character(SC5314.Mate.df$df$entrezID)))

# get the membership of the genes in the regions
Calbi.upset <- make_upsetbase(Calbi.venn$table, c("BWP17", "SC5314.ECE", "SC5314.Mate"))
# get the centrum of the regions
Calbi.arranged <- arrange_venn(Calbi.upset, c("BWP17", "SC5314.ECE", "SC5314.Mate"),
                               extract_regions = T) 

set.seed(123)
# randomize the position of the genes in the regions
xy = rbind(
  calculateCircle(x = 0.00, y = 0.2886, r = 0.1, # 3-way intersection
                  noiseFun = function(x) (x + rnorm(1,0,0.1)), # add noise
                  steps = 579,randomDist = T, randomFun = rnorm), # 579 genes
  calculateEllipse(x = -0.7, y = -0.1266, a = 0.2, b = .2, angle = -240, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 21, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0.7, y = -0.1266, a = 0.2, b = .2, angle = -120, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 55, randomDist = T, randomFun = rnorm),
  calculateEllipse(x = 0, y = 1.0992, a = 0.1, b = .2, angle = 0, 
                   noiseFun = function(x) (x + rnorm(1,0,0.1)),
                   steps = 166, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 0.00, y = -1.2124, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 124, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 1.30, y = 1.0392, r = .3,  
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 175,randomDist = T, randomFun = rnorm),
  calculateCircle(x = -1.30, y = 1.0392, r = .3, 
                  noiseFun = function(x) (x + rnorm(1,0,0.2)),
                  steps = 69,randomDist = T, randomFun = rnorm)
)

Calbi.upset <- Calbi.upset %>%
  dplyr::mutate(region = dplyr::case_when(
    BWP17 & SC5314.ECE & SC5314.Mate ~ "BWP17-SC5314.ECE-SC5314.Mate",
    BWP17 & !SC5314.ECE & SC5314.Mate ~ "BWP17-SC5314.Mate",
    !BWP17 & SC5314.ECE & SC5314.Mate ~ "SC5314.ECE-SC5314.Mate",
    BWP17 & SC5314.ECE & !SC5314.Mate ~ "BWP17-SC5314.ECE",
    !BWP17 & !SC5314.ECE & SC5314.Mate ~ "SC5314.Mate",
    !BWP17 & SC5314.ECE & !SC5314.Mate ~ "SC5314.ECE",
    BWP17 & !SC5314.ECE & !SC5314.Mate ~ "BWP17"
  ), # relabel the regions
  region = as.factor(region),
  x = xy[,1],
  y = xy[,2],
  ) # add the x and y coordinates of the genes
Calbi.upset <- Calbi.upset %>% 
  # add their log2FoldChange values to the data frame
  dplyr::mutate(log2FoldChange = rowMeans(Calbi.upset[,5:7], na.rm = T)) %>% 
  # calculate mean log2FoldChange values for the intersection genes
  dplyr::arrange(desc(abs(log2FoldChange))*-1) 
  # sort the genes by their log2FoldChange values


(venn_plot <- ( # make the plot
  ggplot(Calbi.upset) 
  + coord_fixed()
  + theme_void() # remove the background
  + geom_point(aes(x = x, y = y, colour = log2FoldChange), size = 1) # add the genes
  + scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey50", guide = "colourbar", aesthetics = "colour"
  ) # add a color gradient to the points
  + geom_venn_region(Calbi.upset, sets = c("BWP17", "SC5314.ECE", "SC5314.Mate"), 
                     alpha = 0.3, show.legend = F) # add the regions
  + geom_venn_circle(Calbi.upset, sets = c("BWP17", "SC5314.ECE", "SC5314.Mate"),
                     size = .5) # add the circles
  + scale_fill_venn_mix(Calbi.upset, sets = c("BWP17", "SC5314.ECE", "SC5314.Mate"),
                        guide='none', highlight=c("BWP17-SC5314.ECE-SC5314.Mate"),
                        inactive_color='NA') # highlight the 3-way intersection
  + geom_venn_label_set(Calbi.upset, outwards_adjust = 2,
                        sets = c("BWP17", "SC5314.ECE", "SC5314.Mate"),
                        fill = alpha("black", .35), 
                        aes(label = c("SC5314.old", "SC5314.new", "BWP17"))) 
  # add the set labels
  + geom_venn_label_region(
    Calbi.upset, sets = c("BWP17", "SC5314.ECE", "SC5314.Mate"),
    aes(label=size), outwards_adjust=1.25, position=position_nudge(y=0.2)) 
  # add the number of genes in each region
  + theme(
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12)
  )))

ggsave(paste(plots_dir,"/venn_C.albicans_strain.png"), plot = venn_plot,
       width = 12, height = 6, units = 'in') # save plot


####################################
# 10.) Create C. albicans hallmark #
####################################
library(rrvgo)
library(GOSim)
library(GOSemSim)

Calbi.genelist <- Calbi.upset %>%
  # select the intersection genes
  dplyr::filter(BWP17 & SC5314.ECE & SC5314.Mate) %>% 
  # convert gene names from ENTREZ format to SYMBOL and ENSEMBL format
  dplyr::mutate(entrezID = mapIds(org.Hs.eg.db, geneID, keytype = "SYMBOL", 
                                 column = "ENTREZID", multiVals = "first"),
               ensembl = mapIds(org.Hs.eg.db, geneID, keytype = "SYMBOL", 
                                column = "ENSEMBL", multiVals = "first")) %>% 
  dplyr::select(geneID, entrezID, ensembl, log2FoldChange, Trend) %>%
  # add p-values to the data frame
  dplyr::left_join(., merge.rec(list(BWP17.ECE.df$sig_df[,c("geneID", "pvalue")],
                                     SC5314.ECE.df$sig_df[,c("geneID", "pvalue")],
                                     SC5314.Mate.df$sig_df[,c("geneID", "pvalue")]),
                                by = "geneID", all = F),
                   by = "geneID") %>%
  dplyr::rowwise(.) %>%
  # calculate mean p-values for the intersection genes
  dplyr::mutate(pvalue = mean(c_across(starts_with("pvalue")), na.rm = TRUE)) %>% 
  dplyr::ungroup(.) %>%
  # add minimum p-value, when the p-value is 0
  dplyr::mutate(pvalue = ifelse(pvalue == 0, pvalue + 1e-288, pvalue),
                # assign significance to the genes
                significance = case_when(
                  Trend == "UP" ~ "up-regulated",
                  Trend == "DOWN" ~ "down-regulated"
                ),
                significance = as.factor(significance)) %>% 
  # select the columns of interest
  dplyr::select(geneID, entrezID, ensembl, log2FoldChange, pvalue, significance) %>% 
  as.data.frame(.) # convert to data frame
 
# Get KEGG over representation for the C. albicans hallmark genes
Calbi.KEGG <- kegg_results(Calbi.genelist$entrezID, HSC.gene, kegg_pathways)

Calbi.KEGG_plot <- make_dotplot(Calbi.KEGG$df) # make KEGG dot plot
ggsave(paste(plots_dir,"/C.albicans_shared_kegg.png"), plot = Calbi.KEGG_plot,
       width = 12, height = 6, units = 'in') # save plot

# Get gene set enrichment analysis for the C. albicans hallmark genes
Calbi.hallmark <- gsea_results(Calbi.genelist, hallmark_sets)

Calbi.hallmark_plot <- make_gseaplot(Calbi.hallmark$gsea_df) # make GSEA plot
ggsave(paste(plots_dir,"/C.albicans_shared_gsea.png"), plot = Calbi.hallmark_plot,
       width = 12, height = 6, units = 'in') # save plot


# extract the expression of the C. albicans hallmark genes, name the vector values
# with the gene ENTREZ ID and sort them in descending order
Calbi.expr <- setNames(c(Calbi.genelist$log2FoldChange),
                       c(Calbi.genelist$entrezID)) %>%
  sort(., decreasing = T)

# get the gen set enrichment analysis for the hallmark genes - biological processes
Calbi.BP <- enrichGO(names(Calbi.expr),
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     ont = "BP")

# get the gen set enrichment analysis for the hallmark genes - molecular functions
Calbi.MF <- enrichGO(names(Calbi.expr),
                     OrgDb = "org.Hs.eg.db",
                     keyType = "ENTREZID",
                     ont = "MF")

Calbi.BP_simM <- make_simMatrix(Calbi.BP, "BP")
Calbi.BP_plot <- make_GO_simplot(Calbi.BP_simM$reduced, Calbi.BP_simM$subset)
ggsave(paste(plots_dir,"/Calbi.BP_plot.png"), plot = Calbi.BP_plot,
       width = 24, height = 14, units = 'in')

Calbi.MF_simM <- make_simMatrix(Calbi.MF, "MF")
Calbi.MF_plot <- make_GO_simplot(Calbi.MF_simM$reduced, Calbi.MF_simM$subset)
ggsave(paste(plots_dir,"/Calbi.MF_plot.png"), plot = Calbi.MF_plot,
       width = 16, height = 10, units = 'in')

Calbi.EMT_genes <- Calbi.hallmark$gsea_df["HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "core_enrichment"]
Calbi.EMT_genes <- unlist(stringr::str_split(Calbi.EMT_genes, "/"))

#################################
# 10.) TGF-b vs Ca EMT hallmark #
#################################
library(GEOquery)
library(limma)

GSE23952 <- getGEO("GSE23952")
GSE23952.expr <- exprs(GSE23952[[1]])

GSE23952.meta <- pData(GSE23952[[1]])
GSE23952.meta <- GSE23952.meta %>%
  rownames_to_column("samples") %>%
  dplyr::mutate(treatment = dplyr::case_when(
                  stringr::str_detect(`treatment group:ch1`, "Untreated") ~ "ctrl",
                  stringr::str_detect(`treatment group:ch1`, "TGF") ~ "TGF.b"
                ),
                treatment = factor(treatment, 
                                   levels = c("ctrl","TGF.b"))) %>%
  select(c("samples","treatment", "title"))

GSE23952.features <- fData(GSE23952[[1]]) %>%
  dplyr::select(c("ID","Gene Symbol","ENTREZ_GENE_ID","Gene Title")) %>%
  setNames(c("ID","geneID","entrezID","Gene Title"))

GSE23952.design <- model.matrix(~0+ GSE23952.meta$treatment)
colnames(GSE23952.design) <- c("ctrl","TGFb")

GSE23952.DE <- diffexp_GEO(GSE23952.expr, GSE23952.design)
GSE23952.DE$sig_df <- GSE23952.DE$df %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::inner_join(., GSE23952.features, by = "ID") %>%
  dplyr::filter(P.Value < 0.05 & abs(logFC) >= 1) %>%
  dplyr::arrange(desc(abs(logFC))) %>%
  dplyr::distinct(geneID, .keep_all = T)

shared_set <- merge(Calbi.genelist, GSE23952.DE$sig_df, 
                    by = "geneID", suffix = c("_Ca", "_TGFb")) %>% 
  dplyr::arrange(desc(logFC_Ca)) %>%
  dplyr::rename(log2FoldChange = logFC_Ca)

#
Calbi.hallmark_expression <- get_CategoryExpressionPlot(Calbi.genelist, shared_set$geneID)
ggsave(paste(plots_dir,"/C.albicans_EMT_set_expr.png"), plot = Calbi.hallmark_expression,
       width = 10, height = 8, units = 'in')

shared_KEGG <- kegg_results(genes_of_interest = shared_set$entrezID_Ca, 
                            background_set = unique(c(BWP17.ECE.entrez$background,
                                                      SC5314.ECE.entrez$background,
                                                      SC5314.Mate.entrez$background)),
                            kegg_df = kegg_pathways)

shared_KEGG_plot <- make_dotplot(shared_KEGG$df)
ggsave(paste(plots_dir,"/C.albicans_EMT_kegg.png"), plot = shared_KEGG_plot,
       width = 8, height = 6, units = 'in')

Calbi.EMT_genes <- shared_set %>%
  dplyr::select(c(3,1,2,5,12)) %>%
  setNames(c("entrezID","geneID","logFC_Ca","logFC_TGFb","Description")) %>%
  dplyr::arrange(desc(logFC_Ca)) %>%
  dplyr::mutate(rank_Ca = seq(1, n())) %>%
  dplyr::arrange(desc(logFC_TGFb)) %>%
  dplyr::mutate(rank_TGFb = seq(1, n())) %>%
  dplyr::arrange(desc(logFC_Ca))
  
#######################################
# 11.) WGCNA analyses of shared genes #
#######################################
library(DESeq2)
library(dplyr)
library(tibble)
library(WGCNA)
library(ggplot2)
library(RCy3)
library(RColorBrewer)
library(gProfileR)


tmp <- list(ECE = require_file(ECE_file, header = T, na.strings = NA),
            Mate = require_file(MATE_file, header = T, row.names = 1, na.strings = NA))
tmp <- lapply(tmp, function(x){
  x %>%
    rownames_to_column("geneID")
})

Ca.count_norm <- merge(tmp$ECE, tmp$Mate, by = "geneID") %>%
  column_to_rownames(., "geneID")
Ca.count_norm <- Ca.count_norm %>%
  round() %>%
  data.frame() %>%
  dplyr::select(!contains("ECE1")) %>%
  dplyr::filter(rowSums(.) >= 50)

Ca.metadata <- data.frame("samples" = colnames(Ca.count_norm)) %>%
  dplyr::mutate(samples = as.factor(samples),
                treatment = dplyr::case_when(
                  stringr::str_detect(samples, "_Bwp17_") ~ "C.albi",
                  stringr::str_detect(samples, "_SC5314_") ~ "C.albi",
                  stringr::str_detect(samples, "_Ca_") ~ "C.albi",
                  T ~ "ctrl"
                ),
                treatment = factor(treatment, 
                                   levels = c("ctrl","C.albi")))

Ca.dds <- DESeqDataSetFromMatrix(
  countData = Ca.count_norm,
  colData = Ca.metadata, 
  design = ~1 # Here we are not specifying a model
)
Ca.dds_norm <- vst(Ca.dds)
Ca.dds_norm <- assay(Ca.dds_norm) %>%
  t()
wgcna_soft <- pickSoftThreshold(Ca.dds_norm,
                         dataIsExpr = TRUE,
                         corFnc = cor,
                         networkType = "signed"
)
wgcna_soft_df <- data.frame(wgcna_soft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

ggplot(wgcna_soft_df, aes(x = Power, y = model_fit, label = Power)) +
  # Plot the points
  geom_point() +
  # We'll put the Power labels slightly above the data points
  geom_text(nudge_y = 0.1) +
  # We will plot what WGCNA recommends as an R^2 cutoff
  geom_hline(yintercept = 0.80, col = "red") +
  # Just in case our values are low, we want to make sure we can still see the 0.80 level
  ylim(c(min(wgcna_soft_df$model_fit), 1.05)) +
  # We can add more sensible labels for our axis
  xlab("Soft Threshold (power)") +
  ylab("Scale Free Topology Model Fit, signed R^2") +
  ggtitle("Scale independence") +
  # This adds some nicer aesthetics to our plot
  theme_classic()

Ca.bwnet <- blockwiseModules(Ca.dds_norm,
                          maxBlockSize = 5000, # What size chunks (how many genes) the calculations should be run in
                          TOMType = "signed", # topological overlap matrix
                          power = 16, # soft threshold for network construction
                          numericLabels = TRUE, # Let's use numbers instead of colors for module labels
                          randomSeed = 42) # there's some randomness associated with this calculation
                                           # so we should set a seed
readr::write_rds(Ca.bwnet,
                 file = file.path(results_dir, "SRP140558_wgcna_results.RDS")
)
Ca_eigengenes <- Ca.bwnet$MEs
Ca_eigenmat <- model.matrix(~ Ca.metadata$treatment)
fit <- limma::lmFit(t(Ca_eigengenes), design = Ca_eigenmat)
fit <- limma::eBayes(fit)


Ca_eigendf <- limma::topTable(fit, number = ncol(Ca_eigengenes)) %>%
  tibble::rownames_to_column("module")

module_df <- Ca_eigengenes %>%
  tibble::rownames_to_column("samples") %>%
  dplyr::inner_join(Ca.metadata, by = c("samples"))


#Analysis of the shared gene set
Ca_module_genes <- tibble::enframe(Ca.bwnet$colors, name = "ensembl", value = "module") %>%
  dplyr::mutate(module = paste0("ME", module),
                module = factor(module, levels = c("ME1","ME2","ME3","ME4","ME11"))) %>% 
  dplyr::inner_join(Calbi.genelist, by = c("ensembl"))

Ca_module_msig <- msigdbr_df %>%
  dplyr::filter(gs_cat == "C5" & gs_subcat == "GO:MF") %>%
  dplyr::inner_join(Ca_module_genes, by = c("gene_symbol" = "geneID"))

Ca_module_mf <- Ca_module_msig %>%
  dplyr::select(c("ensembl","entrezID","gene_symbol","log2FoldChange","pvalue",
                  "module","gs_name","gs_exact_source","gs_description")) %>%
  setNames(c("ensembl","entrezID","geneID","log2FC","pvalue","module","GO","ID","Description")) %>%
  dplyr::group_split(module, .keep = T)

Ca_module_simM <- lapply(Ca_module_mf, make_simMatrix, "MF")

Ca_module_mf_df <- list()
  
for (i in 1:length(Ca_module_mf)){
  Ca_module_mf_df[[length(Ca_module_mf_df)+1]] <- Ca_module_mf[[i]] %>% 
                                     dplyr::inner_join(Ca_module_simM[[i]][["reduced"]][,c("go","parent","term","parentTerm")],
                                                       by = c("ID" = "go"))
}

Ca_module_mf_df <- do.call(rbind, Ca_module_mf_df)
Ca_module_mf_df <- Ca_module_mf_df %>%
  dplyr::distinct(geneID, .keep_all = T) %>%
  merge(., Ca_module_genes[,c("geneID", "module")], 
        by = "geneID", suffixes = c("",""), all = T)

Ca_module_mf_df <- Ca_module_mf_df[,c(1,2,3,4,5,13,7:12)]
openxlsx::write.xlsx(Ca_module_mf_df, paste(results_dir, "Ca_module_mf_df.xlsx", sep = "/"))
Ca_module_stats <- Ca_module_mf_df %>%
  dplyr::group_by(parentTerm, module) %>%
  dplyr::summarise(term = n())
openxlsx::write.xlsx(Ca_module_stats, paste(results_dir, "Ca_module_stats.xlsx", sep = "/"))


library(ggplot2)
library(ggrepel)
Ca_module_mf_plot <- ggplot(Ca_module_stats, 
                            aes(fill=parentTerm, y=term, x=module,
                                label = paste(term, parentTerm, sep = " - "))) + 
  geom_bar(position="stack", stat="identity", colour = "grey25") + 
  geom_label_repel(size = 3, position = position_stack(vjust = 0.5)) + 
  theme(legend.position = 'none')

ggsave(paste(plots_dir,"/Ca_module_mf_plot.png"), plot = Ca_module_mf_plot,
       width = 8, height = 14, units = 'in')


#crate a cytoscape picture of interacting genes
cytoscapePing()
#commandsHelp("help string protein query")
Ca_module_interaction <- paste('string protein query taxonID=9606 limit=150 cutoff=0.9 query="',paste(Ca_module_genes$geneID, collapse=","),'"',sep="")
commandsGET(Ca_module_interaction)

Ca_module_attributes <- getTableColumns(table="node")
loadTableData(Ca_module_genes,table.key.column = "display name",data.key.column = "geneID")  #default data.frame key is row.names
#getTableColumnNames()

#setting a new display style
style.name = "EMT"
defaults.list <- list(NODE_SHAPE="ellipse",
                      NODE_SIZE=60,
                      NODE_FILL_COLOR="#AAAAAA",
                      EDGE_TRANSPARENCY=120)
node.label.map <- mapVisualProperty('node label','display name','p') # p for passthrough; nothing else needed
createVisualStyle(style.name, defaults.list, list(node.label.map))
setVisualStyle(style.name=style.name)

#setting node size
size.values <- c(min(Ca_module_genes$log2FoldChange),0,max(Ca_module_genes$log2FoldChange))
setNodeSizeMapping(table.column = 'log2FoldChange', 
                   table.column.values = size.values, 
                   sizes = c(30, 60, 150), mapping.type = "c", style.name = style.name)
setNodeShapeBypass(node.names = EMT_suid, new.shapes = "TRIANGLE")


color.values <- levels(Ca_module_genes$module)
node.colors <- c(brewer.pal(length(color.values), "Set1"))
setNodeColorMapping(table.column = "module", color.values, node.colors, 
                    style.name="EMT", mapping.type = "d", default.color = "#dddddd")

#different node shape for EMT genes
EMT_suid <- Ca_module_attributes$SUID[which(Ca_module_attributes$`display name` %in% shared_set$geneID)]
setNodeBorderWidthBypass(EMT_suid, 20)
setNodeBorderColorBypass(EMT_suid,  "#880808")

#export image
Ca_module_png_file <- file.path(getwd(),plots_dir, "EMT_module_network2.png")
if(file.exists(Ca_module_png_file)){
  file.remove(Ca_module_png_file)
} 

 #ME1 GO enrichments
ME1_nodes <- selectNodes("ME1", by.col="module")
ME1_subnet <- createSubnetwork(nodes="selected")
renameNetwork("ME1_Subnetwork", network=as.numeric(ME1_subnet))

ME1_genes <- getTableColumns(table= "node",network = as.numeric(ME1_subnet))[,"geneID"]
ME1_go <- gProfileR::gprofiler(ME1_genes, organism = "hsapiens",
                    significant=T, ordered_query = F,
                    exclude_iea= T ,max_set_size = 500,
                    min_set_size = 25, correction_method = "fdr",
                    src_filter = c("KEGG","GO:BP"))

ME1_go_df <- ME1_go %>%
  dplyr::filter(overlap.size >= 8) %>%
  dplyr::mutate(phenotype = 1, 
                q.value = p.value) %>%
  dplyr::select(c("term.id","term.name","p.value", "q.value",
                  "phenotype", "intersection")) %>%
  setNames(c("Name","Description","pvalue","qvalue","phenotype","genes"))

ME1_filename <-file.path(getwd(), results_dir, paste("gprofiler_cluster_ME1.txt",sep="_"))

write.table(ME1_go_df,ME1_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

#ME1 GO enrichments
ME2_nodes <- selectNodes("ME2", by.col="module")
ME2_subnet <- createSubnetwork(nodes="selected")
renameNetwork("ME2_Subnetwork", network=as.numeric(ME2_subnet))

ME2_genes <- getTableColumns(table= "node",network = as.numeric(ME2_subnet))[,"geneID"]
ME2_go <- gProfileR::gprofiler(ME2_genes, organism = "hsapiens",
                               significant=T, ordered_query = F,
                               exclude_iea= T ,max_set_size = 500,
                               min_set_size = 25, correction_method = "fdr",
                               src_filter = c("KEGG","GO:BP"))

ME2_go_df <- ME2_go %>%
  dplyr::filter(overlap.size >= 8) %>%
  dplyr::mutate(phenotype = 1, 
                q.value = p.value) %>%
  dplyr::select(c("term.id","term.name","p.value", "q.value",
                  "phenotype", "intersection")) %>%
  setNames(c("Name","Description","pvalue","qvalue","phenotype","genes"))

ME2_filename <-file.path(getwd(), results_dir, paste("gprofiler_cluster_ME2.txt",sep="_"))

write.table(ME2_go_df,ME2_filename,col.name=TRUE,sep="\t",row.names=FALSE,quote=FALSE)

#######################################
# 12.) TCGA univariate COX-regression #
#######################################
library(TCGAbiolinks)
library(org.Hs.eg.db)
library(DESeq2)
library(dplyr)
library(tibble)
library(survival)
library(survminer)
library(ggplot2)
library(ggfortify)
library(MASS)

tcga_path <- file.path(paste0(data_dir, "/tcga_sample_sheet.tsv"))
tcga_file <- read.csv(tcga_path, header = T, sep = "\t")
tcga_samples <- as.character(tcga_file$Sample.ID)
 
query <- GDCquery(
   project = "TCGA-HNSC",
   data.category = "Transcriptome Profiling",
   data.type = "Gene Expression Quantification",
   barcode = tcga_samples
)
GDCdownload(query)
TCGA.HNSC <- GDCprepare(query)

TCGA.HNSC.matrix <- assay(TCGA.HNSC, "unstranded") 

TCGA.features <- as.data.frame(TCGA.HNSC@rowRanges@elementMetadata)
TCGA.metafile <- read.csv(paste0(data_dir,"/tcga_meta/gdac/HNSC.clin.merged.txt"), 
                          header = F, row.names = 1, sep = "\t")
TCGA.clinical <- read.csv(paste0(data_dir,"/tcga_meta/clinical.tsv"), header = T, sep = "\t")
TCGA.coldata <- as.data.frame(TCGA.HNSC@colData)

TCGA.clinical.df <- TCGA.clinical %>%
  dplyr::filter(seq(1,dim(TCGA.clinical)[1]) %% 2 != 1) %>%
  dplyr::mutate(vital_status = ifelse(vital_status == "Dead",T, F),
                ajcc_clinical_m = ifelse(ajcc_clinical_m == "'--", NA, ajcc_clinical_m),
                ajcc_clinical_n = ifelse(ajcc_clinical_n == "'--", NA, ajcc_clinical_n),
                ajcc_clinical_t = ifelse(ajcc_clinical_t == "'--", NA, ajcc_clinical_t),
                ajcc_clinical_stage = ifelse(ajcc_clinical_stage == "'--", NA, ajcc_clinical_stage),
                time = ifelse(days_to_last_follow_up == "'--", days_to_death, days_to_last_follow_up)) %>%
  dplyr::select(c("case_submitter_id","gender","vital_status","ajcc_clinical_m",
                  "ajcc_clinical_n","ajcc_clinical_t","ajcc_clinical_stage",
                  "time","primary_diagnosis")) %>%
  setNames(c("key","gender","status","clinical_m",
             "clinical_n","clinical_t","stage",
             "time","primary_diagnosis"))

# TCGA.clinical.df <- TCGA.clinical.df %>%
#   dplyr::filter(complete.cases(.),
#                 time > 0)
TCGA.tissue.df <- TCGA.coldata %>%
  dplyr::select(c("barcode","patient", "shortLetterCode","definition",
                  "tissue_or_organ_of_origin")) %>%
  setNames(c("barcode","key", "type","definition",
             "tissue_origin"))

TCGA.drug.df <- TCGA.metafile %>%
  t(.) %>%
  as.data.frame(.) %>%
  dplyr::select(c("patient.bcr_patient_barcode","patient.hpv_test_results.hpv_test_result.hpv_status") | 
                  starts_with("patient.drugs.drug.")) %>%
  dplyr::mutate(patient.bcr_patient_barcode = toupper(patient.bcr_patient_barcode)) %>%
  dplyr::rename_with(., ~stringr::str_replace(.x, 'patient.drugs.drug.','')) %>%
  dplyr::select(c("patient.bcr_patient_barcode","patient.hpv_test_results.hpv_test_result.hpv_status",
                  "drug_name","clinical_trail_drug_classification",
                  "number_cycles","prescribed_dose","prescribed_dose_units", 
                  "measure_of_response","regimen_indication","therapy_types.therapy_type",
                  "therapy_ongoing","year_of_form_completion")) %>%
  setNames(c("key","hpv_status","drug_name","drug_classification","nr_cycles","dose","dose_units",
             "response","indicator","therapy","therapy_ongoing","year_of_completion"))


TCGA.coldata.df <- TCGA.clinical.df %>%
  dplyr::inner_join(TCGA.tissue.df, by = "key") %>%
  dplyr::inner_join(TCGA.drug.df, by = "key") %>%
  dplyr::mutate(gender = factor(gender, levels = c('male','female'), labels = c('male','female')),
                clinical_m = as.factor(clinical_m),
                clinical_m = as.factor(clinical_m),
                clinical_n = as.factor(clinical_n),
                clinical_t = as.factor(clinical_t),
                stage = as.factor(stage),
                time = as.numeric(time),
                primary_diagnosis = as.factor(primary_diagnosis),                     
                type = as.factor(type),
                tissue_origin = as.factor(tissue_origin),
                hpv_status = as.factor(hpv_status),
                drug_name = as.factor(drug_name),
                drug_classification = as.factor(drug_classification),
                nr_cycles = as.factor(nr_cycles),
                dose = as.numeric(dose),
                response = as.factor(response),            
                indicator = as.factor(indicator), 
                therapy = as.factor(therapy),
                therapy_ongoing = as.factor(therapy_ongoing),
                thyear_of_completionerapy = as.factor(year_of_completion)) %>%
  dplyr::filter(time > 0 & hpv_status == "negative")

TCGA.coldata.filter <- TCGA.coldata.df$key %in% TCGA.coldata.df$key[which(duplicated(TCGA.coldata.df$key))]
TCGA.coldata.unpaired <- TCGA.coldata.df[!TCGA.coldata.filter,]
TCGA.coldata.paired <- TCGA.coldata.df[TCGA.coldata.filter,] %>%
  dplyr::filter(key != "TCGA-UF-A71A")

TCGA.HNSC.df <- TCGA.HNSC.matrix %>%
  as.data.frame(.) %>%
  dplyr::filter(rowSums(.) > 50) %>%
  tibble::rownames_to_column("geneID") %>% 
  dplyr::mutate(geneID = gsub("\\..*", "", geneID)) %>%
  dplyr::filter(geneID %in% mapIds(org.Hs.eg.db, shared_set$geneID, 
                                   keytype = "SYMBOL",column = "ENSEMBL",
                                   multiVals = "first")) %>%
  tibble::column_to_rownames("geneID") %>%
  t(.) %>%
  as.data.frame(.) 

TCGA.HNSC.unpaired <- TCGA.HNSC.df %>%
  tibble::rownames_to_column("barcode") %>% 
  dplyr::inner_join(TCGA.coldata.unpaired[,c("status","time","key","barcode")], by = "barcode") %>%
  dplyr::select(!c(barcode)) %>%
  dplyr::relocate(where(is.factor), .before = where(is.numeric)) %>%
  dplyr::relocate(where(is.logical), .before = everything()) %>%
  dplyr::relocate(time, .after = status) %>%
  dplyr::distinct(key,.keep_all = T) %>%
  tibble::column_to_rownames("key")

TCGA.HNSC.filter <-  lapply(TCGA.HNSC.unpaired[,-c(1:2)],function(x){
  tmp <- sort(x, decreasing = T)
  l <- round(length(x)*.33)

  high <- tmp[l]
  low <- tmp[length(tmp)-l]

  r <- ifelse(x > high, T, ifelse(x < low, T, F))
  return(r)
})
lapply(TCGA.HNSC.filter, summary)

covariates <- colnames(TCGA.HNSC.unpaired[,3:length(TCGA.HNSC.unpaired)])
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))

TCGA.HNSC.unpaired.base <- TCGA.HNSC.unpaired

TCGA.HNSC.unpaired[,3:length(TCGA.HNSC.unpaired)] <- sapply(names(TCGA.HNSC.unpaired[,3:length(TCGA.HNSC.unpaired)]), function(x){
  TCGA.HNSC.unpaired[,x] <- ifelse(TCGA.HNSC.filter[[x]] == T,ifelse(TCGA.HNSC.unpaired[,x] > median(TCGA.HNSC.unpaired[,x]), 1,-1),0)
})


TCGA.univ_models <- lapply(covariates, function(x){coxph(univ_formulas[[x]], 
                                                         data = TCGA.HNSC.unpaired,
                                                         subset = TCGA.HNSC.filter[[x]])})
names(TCGA.univ_models) <- covariates
# Extract data 
TCGA.univ_results <- lapply(TCGA.univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
TCGA.univ_results <- t(as.data.frame(TCGA.univ_results, check.names = FALSE))
TCGA.univ_results <- as.data.frame(TCGA.univ_results) %>%
  tibble::rownames_to_column("ensembl") %>%
  tidyr::separate(`HR (95% CI for HR)`, sep = " ", into = c("HR","CI (.95)")) %>%
  dplyr::mutate(geneID = mapIds(org.Hs.eg.db, ensembl, keytype = "ENSEMBL", column = "SYMBOL"),
                beta = as.numeric(beta),
                HR = as.numeric(HR),
                wald.test = as.numeric(wald.test),
                p.value = as.numeric(p.value)) %>%
  dplyr::inner_join(shared_set[,c("geneID", "log2FoldChange")], by = "geneID") %>%
  dplyr::filter((HR > 1 & log2FoldChange > 1) | (HR < 1 & log2FoldChange < 1))



#########################################
# 13.) TCGA multivariate COX-regression #
#########################################
library(BiocParallel)
library(DESeq2)
library(fdrtool)
library(reshape2)
library(ggplot2)
library(ggbiplot)
library(ggrepel)
library(genefilter)
library(pheatmap)
TCGA.paired.counts <- TCGA.HNSC.matrix[,which(colnames(TCGA.HNSC.matrix) %in% TCGA.coldata.paired$barcode)]
TCGA.paired.counts <- TCGA.paired.counts %>%
  as.data.frame(.) %>%
  dplyr::filter(rowSums(.) > 50) 
TCGA.paired.counts <- TCGA.paired.counts[,match(TCGA.coldata.paired$barcode, colnames(TCGA.paired.counts))]

TCGA.paired.dds <- calc_DiffExp(matrix = TCGA.paired.counts,
                                   coldata = TCGA.coldata.paired,
                                   design = "1")


#pca of pairwise samples
(TCGA.paired.pca <- prcomp(t(assay(TCGA.paired.dds$dds_norm)), center =T, scale. = TRUE) %>%
    ggbiplot::ggbiplot(., choices = 1:2, scale = 1,
             groups = TCGA.paired.dds$dds_norm@colData$type, var.axes = F, circle = T) +
    geom_point(size = 3, aes(color = TCGA.paired.dds$dds_norm@colData$type)) +
    geom_label_repel(aes(label = TCGA.paired.dds$dds_norm@colData$key, 
                         color = TCGA.paired.dds$dds_norm@colData$type), size = 5) +
    scale_color_discrete(labels=c('Primary tumor','Normal tissue')) +
    labs(color = "Sample type") + 
    stat_ellipse(geom = "polygon", aes(color = TCGA.paired.dds$dds_norm@colData$type,
                                       fill = TCGA.paired.dds$dds_norm@colData$type),
               type = "norm", level = 0.68, alpha = .1, linetype = 2, show.legend = F))

ggsave(paste(plots_dir,"/TCGA_paired_samples_pca.png"), plot = TCGA.paired.pca,
       width = 10, height = 10, units = 'in')

#heatmap of genes with highes variance
mat <- assay(TCGA.paired.dds$dds_norm)
mat <- mat - rowMeans(mat)
mat <- mat %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column("geneID") %>% 
  dplyr::mutate(geneID = gsub("\\..*", "", geneID)) %>%
  dplyr::mutate(geneID = mapIds(org.Hs.eg.db, geneID, 
                                keytype = "ENSEMBL",column = "SYMBOL",
                                multiVals = "first")) %>%
  dplyr::filter(geneID %in% shared_set$geneID) %>%
  tibble::column_to_rownames("geneID") %>%
  as.matrix(.)
col_annot <- TCGA.paired.dds$dds_norm@colData %>%
  as.data.frame(.) %>%
  dplyr::select(type)
row_annot <- Ca_module_genes %>%
  dplyr::filter(geneID %in% row.names(mat)) %>%
  dplyr::select(c(geneID,module)) %>%
  tibble::column_to_rownames("geneID")
(TCGA.paired.heatmap <- pheatmap(mat, show_colnames = F, 
                                clustering_distance_cols = "euclidean",
                                clustering_distance_rows = "correlation", 
                                clustering_method = "ward.D2", 
                                cutree_cols = 2,cutree_rows = 2,
                                annotation_col = col_annot,
                                annotation_row = as.data.frame(row_annot),
                                annotation_colors = list(
                                  type=c("NT" = "lightgreen",
                                         "TP" = "maroon"),
                                  module=c("ME1" = "#90e0ef",
                                           "ME2" = "blue",
                                           "ME3" = "#00b4d8",
                                           "ME4" = "#03045e")
                                )))
ggsave(paste(plots_dir,"/TCGA_paired_heatmap.png"), plot = TCGA.paired.heatmap,
       width = 16, height = 12, units = 'in')

groups <- matrix(data = c(13,26,19,20,1,0,0,1), nrow = 4, ncol = 2, byrow = T,
                 dimnames = list(c("me1","me2","me3","me4"),c("cl1","cl2")))
chisq.test(groups)

register(SnowParam())
TCGA.paired.NT_vs_TP <- calc_DiffExp(matrix = TCGA.paired.counts,
                                     coldata = TCGA.coldata.paired,
                                     design = "type")
TCGA.paired.NT_vs_TP.df <- get_results(TCGA.paired.NT_vs_TP$dds,
                                    contrast = list("type_TP_vs_NT"),
                                    name = c(''))

TCGA.paired.NT_vs_TP.df <- TCGA.paired.NT_vs_TP.df[["resLFC"]] %>%
  as.data.frame(.) %>%
  tibble::rownames_to_column("ensembl") %>% 
  dplyr::mutate(ensembl = gsub("\\..*", "", ensembl)) %>%
  dplyr::distinct(ensembl, .keep_all = T) %>%
  tibble::column_to_rownames("ensembl")

TCGA.paitents.results.df <- sig_results(results = TCGA.paired.NT_vs_TP.df,
                                        sig_log2FC = 1.5,
                                        sig_pval = 0.05)

TCGA.paitents.results.df <- lapply(TCGA.paitents.results.df, function(x){
  x %>%
    dplyr::inner_join(TCGA.univ_results[,c("geneID", "ensembl")], by = "geneID")
})
TCGA.patients.LFC <- 2^TCGA.paitents.results.df$df$log2FoldChange
names(TCGA.patients.LFC) <- TCGA.paitents.results.df$df$geneID


#multivariable cox regression
TCGA.HNSC.paired <- TCGA.HNSC.df %>%
  dplyr::select(contains(TCGA.univ_results$ensembl)) %>%
  tibble::rownames_to_column("barcode") %>% 
  dplyr::inner_join(TCGA.coldata.unpaired[,c("status","time","key","barcode")], by = "barcode") %>%
  dplyr::select(!c(barcode)) %>%
  dplyr::relocate(where(is.factor), .before = where(is.numeric)) %>%
  dplyr::relocate(where(is.logical), .before = everything()) %>%
  dplyr::relocate(time, .after = status) %>%
  dplyr::distinct(key,.keep_all = T) %>%
  tibble::column_to_rownames("key")

treshold <- colMeans(TCGA.HNSC.paired[,3:length(TCGA.HNSC.paired)]) * TCGA.patients.LFC
TCGA.univ_results <- TCGA.univ_results %>% tibble::column_to_rownames("ensembl")
rows <- row.names(TCGA.HNSC.paired)

TCGA.HNSC.paired[,3:length(TCGA.HNSC.paired)] <- sapply(names(TCGA.HNSC.paired[,3:length(TCGA.HNSC.paired)]), function(x){
  TCGA.HNSC.paired[,x] <- ifelse(TCGA.univ_results[x,"HR"] > 1 & TCGA.HNSC.paired[,x] > treshold[x], 1,
                                 ifelse(TCGA.univ_results[x,"HR"] < 1 & TCGA.HNSC.paired[,x] < treshold[x], 1, 0))
})
row.names(TCGA.HNSC.paired) <- rows
TCGA.HNSC.paired <- TCGA.HNSC.paired %>%
  dplyr::mutate(status = ifelse(time > 5*365, F, status)) %>%
  dplyr::mutate(time = ifelse(time > 5*365, 5*365, time))

# Stepwise regression model
TCGA.multiv <-  coxph(Surv(time, status) ~ ., data = TCGA.HNSC.paired)
TCGA.step.model <- stepAIC(TCGA.multiv, direction = "both", 
                           trace = F)

TCGA.step_results <- summary(TCGA.step.model)
TCGA.step_results <- as.data.frame(TCGA.step_results$coefficients, 
                                     check.names = F) %>%
  dplyr::filter(`Pr(>|z|)` < 0.05)

TCGA.step.genes <- mapIds(org.Hs.eg.db,row.names(TCGA.step_results), keytype = "ENSEMBL", column = "SYMBOL")
names(TCGA.step.genes) <- row.names(TCGA.step_results)

set.seed(123)
TCGA.step_plot <- survminer::ggforest(coxph(Surv(time, status) ~ ENSG00000006327 + ENSG00000113070 + 
                                              ENSG00000133816 + ENSG00000144655 + ENSG00000183691,
                                              data = TCGA.HNSC.paired))

mapIds(org.Hs.eg.db, c("ENSG00000006327","ENSG00000113070","ENSG00000133816",
                "ENSG00000144655","ENSG00000183691"), keytype = "ENSEMBL", column = "SYMBOL")

ggsave(paste(plots_dir,"/TCGA_multiv_forest.png"), plot = TCGA.step_plot,
       width = 10, height = 10, units = 'in')


surv.plots <- list()
for(i in 1:length(TCGA.step.genes)){
  #fitting KM regression
    gene <- TCGA.step.genes[i]
    ensembl <- names(TCGA.step.genes[i])
    fit <- survfit(as.formula(paste("Surv(time, status) ~", ensembl)),
                 data = TCGA.HNSC.paired)
    #making the plot
    #png(file.path(plots_dir,paste0(TCGA.genes[i],"_surv.png")),
    #    width = 12, height = 8, units = 'in', res = 300)  
    plot <- ggsurvplot(fit, 
                       data = TCGA.HNSC.paired,
                       pval = T, pval.method = T, conf.int = T, 
                       risk.table = "abs_pct",
                       #add.all = T,
                       legend.labs = c(
                         #paste("OS"),
                         paste(gene, " = low"),
                         paste(gene, " = high")),
                       ggtheme = theme_bw()) 
      
    #dev.off()
    surv.plots[[i]] <- plot 
}

TCGA.HNSC.hazard <- rowSums(TCGA.HNSC.paired[,-2L:-1L])/52
TCGA.HNSC.hazard.df <- cbind(TCGA.HNSC.paired[,1:2],TCGA.HNSC.hazard) %>%
  setNames(c("status","time","hazard")) %>%
  dplyr::mutate(hazard = ifelse(hazard > 0.28846, 1, -1))


TCGA.HNSC.hazard.fit <- survfit(Surv(time, status) ~ hazard, 
                                      data = TCGA.HNSC.hazard.df)


TCGA.HNSC.hazard.plot <- ggsurvplot(TCGA.HNSC.hazard.fit, 
                   data = TCGA.HNSC.hazard.df,  surv.median.line = "hv",
                   pval = T, pval.method = T, conf.int = T,
                   risk.table = "abs_pct",
                   legend.labs = c(
                     "HR-","HR+"),
                   ggtheme = theme_bw(),
                   break.time.by = c(365))

png(file.path(plots_dir,"hazard_surv_plot.png"),
    width = 12, height = 8, units = 'in', res = 300)  
TCGA.HNSC.hazard.plot
dev.off()

hazard.cox <- coxph(Surv(time, status) ~ hazard, 
        data = TCGA.HNSC.hazard.df)
TCGA.HNSC.hazard.stat <- summary(hazard.cox)

dt = sig_results(TCGA.paitents.results$type_TP_vs_NT, 1, 0.05)

marker_plot <- (get_CategoryExpressionPlot(dt$df, 1, 0.05, TCGA.step.genes))
ggsave(paste(plots_dir,"/marker_expr.png"), plot = marker_plot,
       width = 10, height = 8, units = 'in')
