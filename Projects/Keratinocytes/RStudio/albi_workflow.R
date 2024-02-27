setwd("D:/Adam/Fungi/RStudio")
source("../R_functions.R")

#data folder
data_dir <- "data"
if (!dir.exists("data")) {
  dir.create("data")
}

#plots directory
plots_dir <- "pictures"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

#results directory
results_dir <- "tables" 
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}


##########################################################
library(dplyr)
library(stringr)
library(openxlsx)
library(countToFPKM)
#copy count files to data dir
albi.paths <- list.files("../Counts/", pattern = "counts.txt$", full.names = T)
albi.reads <- lapply(albi.paths, read.csv, header = T, sep = "\t", skip = 1)
albi.counts <- lapply(albi.reads, function(x){
  x %>% select(c(1,7))
})
gene.lengths <- lapply(albi.reads, function(x){
  return(setNames(x$Length, x$Geneid))
})
gene.lengths <- gene.lengths[[1]]



#merge & write count table
albi.counts <- merge.rec(albi.counts, by = "Geneid", all = T)
albi.counts <- albi.counts %>%
  tibble::column_to_rownames("Geneid") %>%
  setNames(str_sub(colnames(.), 7L, -10L))

albi.fpkm <- fpkm(as.matrix(albi.counts), gene.lengths, rep(150,12)) 

albi.countTable <- merge(albi.counts %>% as.data.frame() %>% tibble::rownames_to_column("Geneid"),
                         albi.fpkm %>% as.data.frame() %>% tibble::rownames_to_column("Geneid"),
                         by = "Geneid", suffixes = c("_RC","_FPKM"))
write.xlsx(albi.countTable, paste(results_dir,"total_readcounts.xlsx",sep = "/"),
           colNames = T, rowNames = T)

library(stringr)
library(DESeq2)
library(edgeR)

strain <- as.factor(c(rep("SC5314",6), rep("WO1",6)))
condition <- as.factor(rep(c(rep("ctrl",3),rep("HK",3)),2))
rep <- as.factor(rep(paste("run", seq(1,3), sep = ""), 4))
experiment <- as.factor(str_glue("{strain}_{condition}",
                                 strain = strain,
                                 condition = condition))
sample <- as.factor(str_glue("{strain}_{condition}_{rep}",
                             strain = strain,
                             condition = condition,
                             rep = rep))
coldata <- data.frame(
  "strain" = strain,
  "condition" = condition,
  "rep" = rep,
  "experiment" = experiment,
  "sample" = sample
)

albi.deseq <- make_deseq(matrix = albi.counts, 
                         coldata =  coldata,
                         design = "strain + condition")
strain.deseq <- make_deseq(matrix = albi.counts, 
                         coldata =  coldata,
                         design = "experiment")

##########################################################
library(ggplot2)
library(ggbiplot)
library(ggrepel)
library(FactoMineR)
albi.pca <- make_pca(albi.deseq$rld, group = experiment)
ggsave(paste(plots_dir,"albicans_pca.png",sep="/"), plot = albi.pca,
       width = 8, height = 8, units = 'in')

library(fdrtool)
library(openxlsx)
# cell_effect.res <- get_results(albi.deseq$dds, 0, 0.05, contrast = list("condition_HK_vs_ctrl"), name = c(''))
strain.res <- get_results(strain.deseq$dds, 1.5, 0.05, contrast = c("experiment","SC5314_HK","WO1_HK"), name = c(''))
SC5314.res <- get_results(strain.deseq$dds, 1.5, 0.05, contrast = c("experiment","SC5314_HK","SC5314_ctrl"), name = c(''))
write.xlsx(SC5314.res$sig_df, paste(results_dir,"SC5314_degs.xlsx",sep = "/"))

WO1.res <- get_results(strain.deseq$dds, 1.5, 0.05, contrast = c("experiment","WO1_HK","WO1_ctrl"), name = c(''))
write.xlsx(WO1.res$sig_df, paste(results_dir,"WO1_degs.xlsx",sep = "/"))

##########################################################
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
mat <- rlog(as.matrix(albi.counts))
count.mat <- make_matrix(mat, unique(c(rownames(SC5314.res$sig_df),
                                       rownames(WO1.res$sig_df))))#,
                                       #rownames(strain.res$sig_df))))


#c("silhouette", "wss", "gap_stat")
fviz_nbclust(count.mat, kmeans, method = "wss")
annot <- HeatmapAnnotation(Strain = as.factor(coldata$strain),
                           Experiment = as.factor(coldata$condition),
                           show_annotation_name = F)

#albi.h1 <- 
#method1 (dist) - "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski"  
#method2 (hclust) - "ward.D", "ward.D2", "complete", "average", "mcquitty"
albi.h1 <- expression_heatmap(mat = count.mat, columns = 2, clusters = 4, 
                   coldata = coldata, annot = annot, 
                   method1 = "euclidean", method2 = "complete") 

# albi.h2 <- expression_heatmap(mat = count.mat, clusters = 5, 
#                    coldata = coldata, annot = annot, 
#                    method1 = "minkowski", method2 = "ward.D2")
library(dplyr)
library(plyr)
library(tidyr)
library(tibble)
expr.mat <- list(
  data.frame("geneID" = SC5314.res$sig_df$geneID,
             `CTRL vs HK (SC5314)` = SC5314.res$sig_df$log2FoldChange, 
             row.names = row.names(SC5314.res$sig_df)),
  data.frame("geneID" = WO1.res$sig_df$geneID,
             `CTRL vs HK (WO1)` = WO1.res$sig_df$log2FoldChange,
             row.names = row.names(WO1.res$sig_df)))#,
  #data.frame("geneID" = cell_effect.res$sig_df$geneID,
  #           `CTRL vs HK` = cell_effect.res$sig_df$log2FoldChange,
  #           row.names = row.names(cell_effect.res$sig_df))),
  # data.frame("geneID" = strain.res$sig_df$geneID,
  #             `SC5314 vs WO1` = strain.res$sig_df$log2FoldChange,
  #             row.names = row.names(strain.res$sig_df)))

expr.mat <- merge.rec(expr.mat, by="geneID", all=T, suffixes=c("", ""))
expr.mat <- expr.mat %>% column_to_rownames("geneID") %>%
  setNames(., c("CTRL vs HK (SC5314)", "CTRL vs HK (WO1)"))#, "SC5314 vs WO1")) #"CTRL vs HK", 


albi.expr1 <- logFC_heatmap(mat = as.matrix(expr.mat), dend = albi.h1$dendrogram)

svg(paste(plots_dir,"total_heatmap.svg",sep="/"), width = 14, height = 8, onefile = T)
draw(albi.h1$heatmap + albi.expr1, merge_legend = T)
dev.off()

png(paste(plots_dir,"total_heatmap.png",sep="/"), width = 14, height = 8, units = 'in', res = 300)
draw(albi.h1$heatmap + albi.expr1, merge_legend = T)
dev.off()

##########################################################
library(openxlsx)
write.xlsx(expr.mat, paste(results_dir,"genes.xlsx",sep = "/"), rowNames = T)
cluster1_genes <- labels(albi.h1$dendrogram[[1]][[1]][[1]])
cluster2_genes <- labels(albi.h1$dendrogram[[1]][[1]][[2]])
cluster3_genes <- labels(albi.h1$dendrogram[[1]][[2]])
cluster4_genes <- labels(albi.h1$dendrogram[[2]])

fpkm.mat <- make_matrix(albi.fpkm, row.names(count.mat))
fpkm.df <- as.data.frame(fpkm.mat) %>%
  rownames_to_column("Geneid") %>%
  pivot_longer(!Geneid,
               names_to = "Condition",
               values_to = "FPKM") %>%
  dplyr::mutate(Cluster = case_when(
    Geneid %in% cluster1_genes ~ "CL1",
    Geneid %in% cluster2_genes ~ "CL2",
    Geneid %in% cluster3_genes ~ "CL3",
    Geneid %in% cluster4_genes ~ "CL4"
  ),
  Cluster = as.factor(Cluster),
  Condition = as.factor(Condition)) %>%
  dplyr::relocate(where(is.numeric), .after = everything())

fpkm.df$Condition <- factor(fpkm.df$Condition, levels = c("WO_C_3","WO_C_1","WO_C_2","WO_HK_3","WO_HK_1","WO_HK_2",
                                                         "SC5314_C_3","SC5314_C_2","SC5314_C_1","SC5314_HK_2","SC5314_HK_1","SC5314_HK_3"))  
cluster_plot <- ggplot(fpkm.df, aes(x = Condition, y = FPKM, group = Geneid)) + 
  geom_line(color = "grey", linewidth = 0.7) + 
  facet_wrap(~ Cluster, ncol = 4, scales = "free", labeller = "label_both") +
  geom_line(aes(y=with(fpkm.df, ave(FPKM, Cluster, Condition, FUN=mean))), 
            linewidth = 1.5, color="steelblue") + 
  theme_minimal() + 
  scale_x_discrete(guide = guide_axis(angle = 60)) + 
  theme(panel.spacing = unit(3, "cm"),
        plot.margin = margin(1,1,1.5,1.2, "cm"),
        strip.background = element_rect(color = "black", fill = "steelblue"),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(color = "grey25", size = 14),
        axis.text = element_text(color = "grey25", size = 14))

ggsave(paste(plots_dir,"cluster_plot.png",sep="/"), plot = cluster_plot,
       width = 16, height = 8, units = 'in')

svg(paste(plots_dir,"cluster_plot.svg",sep="/"), width = 18, height = 10, onefile = T)
(cluster_plot)
dev.off()

##########################################################
library(ComplexUpset)
library(GOplot)
library(ggplot2)
library(ggrepel)
library(conicfit)
venn.base <- expr.mat[,1:2] %>%
  filter(abs(rowSums(., na.rm = T)) > 0) %>%
  mutate(SC5314 = ifelse(!is.na(.$`CTRL vs HK (SC5314)`), T, F),
          WO1 = ifelse(!is.na(.$`CTRL vs HK (WO1)`), T, F)) %>%
  dplyr::rename("SC5314_LFC" = `CTRL vs HK (SC5314)`,
                 "WO1_LFC" = `CTRL vs HK (WO1)`) %>%
  dplyr::relocate(where(is.numeric), .after = where(is.logical))

venn <- GOVenn(venn.base %>% rownames_to_column("Geneid") %>% select(c(1,4)) %>% drop_na(.),
               venn.base %>% rownames_to_column("Geneid") %>% select(c(1,5)) %>% drop_na(.), 
               label = c("SC5314","WO1"), circle.col = c("salmon", "turquoise"), lfc.col = c("darkred","green","blue"))

venn.arranged <- arrange_venn(venn.base, sets = c("SC5314","WO1"), extract_regions = T)
set.seed(23)
xy = rbind(
  calculateCircle(x = 1.3, y = 0, r = 0.2, 
                  noiseFun = function(x) (x + rnorm(1,0,0.1)),
                  steps = 42,randomDist = T, randomFun = rnorm),
  calculateCircle(x = -1.3, y = 0, r = .2, 
                  noiseFun = function(x) (x + rnorm(1,0,0.1)),
                  steps = 44, randomDist = T, randomFun = rnorm),
  calculateCircle(x = 0, y = 0, r = .05, 
                  noiseFun = function(x) (x + rnorm(1,0,0.1)),
                  steps = 21,randomDist = T, randomFun = rnorm)
)

venn.base <- venn.base %>%
  dplyr::mutate(region = dplyr::case_when(
    SC5314 & WO1 ~ "SC5314-WO1",
    SC5314 & !WO1 ~ "SC5314",
    !SC5314 & WO1 ~ "WO1"
  ),
  region = factor(region, levels=c("WO1","SC5314","SC5314-WO1"))) %>%
  dplyr::arrange(region) %>%
  dplyr::mutate(
  x = xy[,1],
  y = xy[,2],
  log2FoldChange = rowMeans(venn.base[,3:4], na.rm = T)) #%>%
  #dplyr::arrange(desc(abs(log2FoldChange))*-1)

(venn_plot <- (
  ggplot(venn.base)
  + coord_fixed()
  + theme_void()
  + geom_point(aes(x = x, y = y, colour = log2FoldChange), size = 4)
  + scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    na.value = "grey50", guide = "colourbar", aesthetics = "colour"
  )
  + geom_venn_region(venn.base, sets = c("SC5314","WO1"), 
                     alpha = 0.0, show.legend = F)
  + geom_venn_circle(venn.base, sets = c("SC5314","WO1"),
                     size = .5)
  + scale_fill_venn_mix(venn.base, sets = c("SC5314","WO1"),
                        guide='none', inactive_color='white')
  + geom_venn_label_set(venn.base, outwards_adjust = 1.5, nudge_y = 1,
                        sets = c("SC5314","WO1"), size = 12,
                        fill = alpha("white", .35), aes(label = c("SC5314","WO1")))
  + geom_venn_label_region(venn.base, sets = c("SC5314","WO1"),aes(label=size), 
                           outwards_adjust=1.25, position=position_nudge(y=0.3),
                           size = 10, fill = alpha("white", .35))
  + theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )))
ggsave(paste(plots_dir,"venn_plot.png",sep = "/"), plot = venn_plot,
       width = 12, height = 8, units = 'in')

##########################################################
library(org.Sc.sgd.db)
library(clusterProfiler)
library(rrvgo)
library(GOSim)
library(stringr)
library(GOSemSim)
library(AnnotationHub)
library(biomaRt)
library(GenomeInfoDb)
library(AnnotationDbi)
library(Biobase)
library(AnnotationForge)
library(GO.db)
library(ggplot2)
library(ggrepel)
library(tibble)
library(ggplot2)
library(dplyr)
library(forcats)
library(numbers)
library(plyr)
library(stringr)

##########################################################
WO1_GO <- list(
  BP = read.xlsx(paste(results_dir,"WO1vsHK_BP.xlsx", sep = "/"), colNames = T),
  MF = read.xlsx(paste(results_dir,"WO1vsHK_MF.xlsx", sep = "/"), colNames = T),
  CC = read.xlsx(paste(results_dir,"WO1vsHK_CC.xlsx", sep = "/"), colNames = T)
)

WO1_simM <- make_simMatrix(WO1_GO)
WO1_scores <- lapply(WO1_GO, function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(WO1_scores)){
  WO1_scores[[var]] <- WO1_scores[[var]][which(names(WO1_scores[[var]]) %in% row.names(WO1_simM[[var]]))]
} 
WO1_reduced <- list()
for(var in names(WO1_simM)){
  WO1_reduced[[var]] <- get_reducedTerms(simm = WO1_simM[[var]], scores = WO1_scores[[var]], 0.9)
} 

WO1_plotBP <- make_GO_simplot(WO1_reduced$BP$reduced, WO1_reduced$BP$subset)
ggsave(paste(plots_dir,"WO1vsHK_BP.png", sep = "/"), WO1_plotBP,
       width = 24, height = 12, units = 'in')
WO1_plotBP_plot <- make_FGplot(WO1_GO$BP, type = "GO")
ggsave(paste(plots_dir,"WO1vsHK_BP_dotplot.png", sep = "/"), WO1_plotBP_plot,
       width = 12, height = 10, units = 'in')

WO1_plotMF <- make_GO_simplot(WO1_reduced$MF$reduced, WO1_reduced$MF$subset)
ggsave(paste(plots_dir,"WO1vsHK_MF.png", sep = "/"), WO1_plotMF,
       width = 24, height = 12, units = 'in')
WO1_plotMF_plot <- make_FGplot(WO1_GO$MF)
ggsave(paste(plots_dir,"WO1vsHK_MF_dotplot.png", sep = "/"), WO1_plotMF_plot,
       width = 12, height = 10, units = 'in')

WO1_plotCC <- make_GO_simplot(WO1_reduced$CC$reduced, WO1_reduced$CC$subset)
ggsave(paste(plots_dir,"WO1vsHK_CC.png", sep = "/"), WO1_plotCC,
       width = 24, height = 12, units = 'in')
WO1_plotCC_plot <- make_FGplot(WO1_GO$CC)
ggsave(paste(plots_dir,"WO1vsHK_CC_dotplot.png", sep = "/"), WO1_plotCC_plot,
       width = 12, height = 10, units = 'in')

WO1_KEGG_df <- read.xlsx(paste(results_dir,"WO1vsHK_KEGG.xlsx", sep = "/"), colNames = T)
WO1_KEGG_plot <- make_FGplot(WO1_KEGG_df)
ggsave(paste(plots_dir,"WO1vsHK_KEGG.png", sep = "/"), WO1_KEGG_plot,
       width = 12, height = 10, units = 'in')

##########################################################
library(GOplot)

for (i in names (WO1_GO)){
  WO1_GO[[i]] <- WO1_GO[[i]] %>%
    dplyr::mutate(Category = i)
}
WO1_KEGG_df$Category = "KEGG"
WO1.df <- do.call(rbind.fill, append(WO1_GO, list(KEGG=WO1_KEGG_df)))
WO1.df$category <- factor(WO1.df$category, levels=c("BP", "MF", "CC", "KEGG")) 

WO1.df <- make_GObase(WO1.df, WO1.res$sig_df)
(WO1.plot <- compound_GOplot(WO1.df))
ggsave(paste(plots_dir,'WO1_GO_plot.png',sep = "/"),plot = WO1.plot,device = "png",width = 18,height = 10, units = 'in')

##########################################################
SC5314_GO <- list(
  BP = read.xlsx(paste(results_dir,"SC5314vsHK_BP.xlsx", sep = "/"), colNames = T),
  MF = read.xlsx(paste(results_dir,"SC5314vsHK_MF.xlsx", sep = "/"), colNames = T),
  CC = read.xlsx(paste(results_dir,"SC5314vsHK_CC.xlsx", sep = "/"), colNames = T)
)

SC5314_simM <- make_simMatrix(SC5314_GO)
SC5314_scores <- lapply(SC5314_GO, function(x) setNames(-log10(x$pvalue), x$ID))
for(var in names(SC5314_scores)){
  SC5314_scores[[var]] <- SC5314_scores[[var]][which(names(SC5314_scores[[var]]) %in% row.names(SC5314_simM[[var]]))]
} 
SC5314_reduced <- list()
for(var in names(SC5314_simM)){
  SC5314_reduced[[var]] <- get_reducedTerms(simm = SC5314_simM[[var]], scores = SC5314_scores[[var]], 0.9)
} 

(SC5314_plotBP <- make_GO_simplot(SC5314_reduced$BP$reduced, SC5314_reduced$BP$subset))
ggsave(paste(plots_dir,"SC5314vsHK_BP.png", sep = "/"), SC5314_plotBP,
       width = 24, height = 12, units = 'in')
(SC5314_plotBP_plot <- make_FGplot(SC5314_GO$BP))
ggsave(paste(plots_dir,"SC5314vsHK_BP_dotplot.png", sep = "/"), SC5314_plotBP_plot,
       width = 12, height = 10, units = 'in')

(SC5314_plotMF <- make_GO_simplot(SC5314_reduced$MF$reduced, SC5314_reduced$MF$subset))
ggsave(paste(plots_dir,"SC5314vsHK_MF.png", sep = "/"), SC5314_plotMF,
       width = 24, height = 12, units = 'in')
(SC5314_plotMF_plot <- make_FGplot(SC5314_GO$MF))
ggsave(paste(plots_dir,"SC5314vsHK_MF_dotplot.png", sep = "/"), SC5314_plotMF_plot,
       width = 12, height = 10, units = 'in')

(SC5314_plotCC <- make_GO_simplot(SC5314_reduced$CC$reduced, SC5314_reduced$CC$subset))
ggsave(paste(plots_dir,"SC5314vsHK_CC.png", sep = "/"), SC5314_plotCC,
       width = 24, height = 12, units = 'in')
(SC5314_plotCC_plot <- make_FGplot(SC5314_GO$CC))
ggsave(paste(plots_dir,"SC5314vsHK_CC_dotplot.png", sep = "/"), SC5314_plotCC_plot,
       width = 12, height = 10, units = 'in')

SC5314_KEGG_df <- read.xlsx(paste(results_dir,"SC5314vsHK_KEGG.xlsx", sep = "/"), colNames = T)
(SC5314_KEGG_plot <- make_FGplot(SC5314_KEGG_df))
ggsave(paste(plots_dir,"SC5314vsHK_KEGG.png", sep = "/"), SC5314_KEGG_plot,
       width = 12, height = 10, units = 'in')

##########################################################
library(GOplot)

for (i in names (SC5314_GO)){
  SC5314_GO[[i]] <- SC5314_GO[[i]] %>%
    dplyr::mutate(Category = i)
}
SC5314_KEGG_df$Category = "KEGG"
SC5314.df <- do.call(rbind.fill, append(SC5314_GO, list(KEGG=SC5314_KEGG_df)))
SC5314.df$Category <- factor(SC5314.df$Category, levels=c("BP", "MF", "CC", "KEGG")) 

SC5314.df <- make_GObase(SC5314.df, SC5314.res$sig_df)
(SC5314.plot <- compound_GOplot(SC5314.df))
ggsave(paste(plots_dir,'SC5314_GO_plot.png',sep = "/"),plot = SC5314.plot,device = "png",width = 18,height = 10, units = 'in')




##########################################################
library(stringr)
library(ontologyIndex)
candida_slim <- get_OBO(paste(data_dir, "candida_slim.obo",sep = "/"),
                        extract_tags = "minimal",
                        merge_equivalent_terms = TRUE)
candida_slim_names <- candida_slim$name

gene_assotiation <- read.csv(paste(data_dir, "gene_association.tsv",sep = "/"), header = F, sep = "\t")
gene_assotiation <- gene_assotiation %>%
  dplyr::select(c(11,5,9,13)) %>%
  setNames(c("geneID","ID","Category","taxon")) %>%
  dplyr::mutate(geneID = str_split_i(geneID, "\\|", 1),
                ID = as.factor(ID),
                Category = case_when(Category == "C" ~ "CC",
                                     Category == "F" ~ "MF",
                                     Category == "P" ~ "BP",),
                Category = as.factor(Category),
                taxon = str_split_i(taxon, ":", 2)) %>%
  dplyr::filter(taxon == "237561")

gene_association_n <- gene_assotiation %>% group_by(ID) %>% summarise(n = n())

slim_gene_assotiation <- read.csv(paste(data_dir, "gene_association_GO_slim.tsv",sep = "/"), header = F, sep = "\t")
slim_gene_assotiation <- slim_gene_assotiation %>%
  dplyr::select(c(11,5,9,13)) %>%
  setNames(c("geneID","ID","Category","taxon")) %>%
  dplyr::mutate(geneID = str_split_i(geneID, "\\|", 1),
                ID = as.factor(ID),
                Category = case_when(Category == "C" ~ "CC",
                                     Category == "F" ~ "MF",
                                     Category == "P" ~ "BP",),
                Category = as.factor(Category),
                taxon = str_split_i(taxon, ":", 2)) %>%
  dplyr::filter(taxon == "5476")
slim_terms <- slim_gene_assotiation %>% pull("ID") %>% levels(.)
slim_gene_assotiation$term <- as.factor(candida_slim_names[match(slim_gene_assotiation$ID, names(candida_slim_names))])
