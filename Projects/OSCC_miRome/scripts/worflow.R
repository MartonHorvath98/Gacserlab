#RNA-Seq data (2019/.../...)    #Horváth Márton
#
#this pipeline goes along the processing of RNA-Seq read starting with
#the .bam files obtained after the alignment step.

.libPaths("C:/Program Files/R/R-4.1.1/library")
setwd("G:/Saját/Labor/Projektek/HSC2_mRNA-Seq/")

library(GenomicFeatures)

#With the help of the "genomicFeatures" package we import a reference transcriptome
#into R that is in a .gtf/.gff3 extension... afterwards, we obtain the coordinates 
#of the exons from the txdb object that contains the annotated transcriptome...
human.txdb <- makeTxDbFromGFF("./grch38/Homo_sapiens.GRCh38.96.gff3", dataSource = "Ensembl", organism = "Homo sapiens")
human.genes <- exonsBy(human.txdb, by = "gene")



#we assign the path to the .bam files and save it in the form of a list, which will 
#later facilitate the work with them...
bam.files <- list.files("./hisat/", pattern = ".bam$", full.names = T)



# With the "Rsamtools" package we convert the files path list to a BamFileList object
# that we can directly use to count reads later on... to help the stability of the 
# processing with a limitid computing capacity (e.g. when using a laptop during the
# calculation) we can set a yield size to limit the number of reads processed contem-
# porarly
library(Rsamtools)
bam.list <- BamFileList(bam.files, yieldSize = 2000000)

library(GenomeInfoDb)
seqinfo(bam.list[1])

# The make calculation more easily processable we can use the 'serialParam' function to
#parameterize serial evaluation, primarily to facilitate easy transition from parallel 
# to serial code
library(BiocParallel)
register(SerialParam())
# Now, with the 'summarizeOverlaps' function of the "GenomicAlignments" package we 
# finally calculate readcounts - the used parameters are:
#
# - features           -> reference exons the align the reads to
# - reads              -> the BamFileList object (our reads)
# - mode               -> = "Union" (default) reads are counted when align with exactly one
#                         exon
# - ingore.strand      -> = FALSE, considering strand direction
# - singleEnd          -> = FALSE, would be TRUE (default), if the sequencing is single-ended
# - fragments          -> considering every read (singletons, etc.), during counting
# - preprocessed.reads -> ='invertStrand', strandspecificity (=firststrand, RF)
library(GenomicAlignments)
human_reads <- summarizeOverlaps(features = human.genes, 
                                 reads = bam.list, 
                                 mode = "Union", 
                                 ignore.strand = F, 
                                 singleEnd = F, 
                                 fragments = T, 
                                 preprocess.reads = invertStrand)

summarizeOverlaps(BPPARAM = SnowParam())

print(head(assay(human_reads), 20))
write.table(assay(human_reads), "./RStudio/total_readcounts.csv", sep = ",", na = "NA", dec = ".", row.names = T, col.names = T)


library(stringr)
bam.names <- as.list(strsplit(bam.files, "/")) #get samplenames
bam.names <- lapply(bam.names, function(x) x[[3]]) #exclude path
bam.time <- as.factor(sapply(bam.names, function(x) substr(x, 1, str_locate(x,"_") - 1))) #separate time
bam.treatment <- as.factor(sapply(bam.names, function(x) substr(x, str_locate(x,"_") + 1, str_locate(x, ".bam") - 2))) #separate treatment
bam.run <- as.factor(sapply(bam.names, function(x) paste("run", (substr(x, str_locate(x, ".bam") - 1, str_locate(x, ".bam") - 1))))) #separate run
bam.condition <- as.factor(paste(bam.treatment, bam.time, sep = "_")) #merge conditions -> treatment + time

coldata <- data.frame(samplenames = as.character(bam.names), time = bam.time, treatment = bam.treatment, run = bam.run, condition = bam.condition) #make metadate table


library(DESeq2)
library(edgeR)
#run the actual differencial expression analysis with several steps (1. leolvasási mélység megállapítása és korrigálása, 
#2. negatív binomiális eloszlási tényező kalkulálása, 3. szórás számítása (max. likelyhood módszerrel), 4. statisztikai teszt
#(Wald teszt))
deseq <- DESeqDataSetFromMatrix(assay(human_reads), colData = coldata, design = ~condition)
dds <- DESeq(deseq)



subset <- rowSums(cpm(counts(dds), normalized = T) >1 ) >= 3
dds <- dds[subset,]
#subset_coll <- rowSums(cpm(counts(dds_coll), normalized = T) >1 ) >= 3
#dds_coll <- dds_coll[subset,]
counts <- as.data.frame(assay(dds))

rld <- rlog(dds, blind = F)
#rld_super_coll <- rlog(dds_super_coll, blind = FALSE)

library(fdrtool)
results_Ca_1h <- results(dds, contrast = c("condition", "CA11_H1", "C_H1"), independentFiltering = T,
                         pAdjustMethod = "BH", alpha = 0.05)
#fdr_results_Ca_1h <- fdrtool(results_Ca_1h$stat, statistic = "normal",)
#results_Ca_1h$padj <- p.adjust(fdr_results_Ca_1h$pval, method = "BH")
summary(results_Ca_1h$padj)

Ca1h_res <- data.frame(results_Ca_1h, geneID = rownames(assay(dds)), row.names = rownames(assay(dds)))

library(openxlsx)
write.xlsx(Ca1h_res, "./RStudio/results/Ca1h_result.xlsx", asTable = T, rowNames = T, colNames = T)

Ca1h_sig.down <- data.frame(results_Ca_1h[which(results_Ca_1h$log2FoldChange < -1.5 & results_Ca_1h$pvalue < 0.05),])
write.xlsx(Ca1h_sig.down, "./RStudio/results/Ca1h_sig.down.xlsx", asTable = T, rowNames = T, colNames = T, sheetName = "Sig.down. genes")

Ca1h_sig.up <- data.frame(results_Ca_1h[which(results_Ca_1h$log2FoldChange > 1.5 & results_Ca_1h$pvalue < 0.05),])
write.xlsx(Ca1h_sig.up, "./RStudio/results/Ca1h_sig.up.xlsx", asTable = T, rowNames = T, colNames = T, sheetName = "Sig. up genes")


results_Ca_6h <- results(dds, contrast = c("condition", "CA11_H6", "C_H6"), independentFiltering = T,
                         pAdjustMethod = "BH", alpha = 0.05)
#fdr_results_Ca_6h <- fdrtool(results_Ca_6h$stat, statistic = "normal",)
#results_Ca_6h$padj <- p.adjust(fdr_results_Ca_6h$pval, method = "BH")
hist(results_Ca_6h$pvalue)

Ca6h_res <- data.frame(results_Ca_6h, geneID = rownames(assay(dds)), row.names = rownames(assay(dds)))
write.xlsx(Ca6h_res, "./RStudio/results/Ca6h_result.xlsx", asTable = T, rowNames = T, colNames = T)


Ca6h_sig.down <- data.frame(results_Ca_6h[which(results_Ca_6h$log2FoldChange < -1.5 & results_Ca_6h$pvalue < 0.05),])
write.xlsx(Ca6h_sig.down, "./RStudio/results/Ca6h_sig.down.xlsx", asTable = T, rowNames = T, colNames = T, sheetName = "Sig.down. genes")

Ca6h_sig.up <- data.frame(results_Ca_6h[which(results_Ca_6h$log2FoldChange > 1.5 & results_Ca_6h$pvalue < 0.05),])
write.xlsx(Ca6h_sig.up, "./RStudio/results/Ca6h_sig.up.xlsx", asTable = T, rowNames = T, colNames = T, sheetName = "Sig. up. genes")

Ca_total_results <- merge(Ca1h_res, Ca6h_res, by = 'geneID', suffixes = c(" (1h)", " (6h)"), all = T)
write.xlsx(Ca_total_results, "./RStudio/results/Ca_total_results.xlsx", asTable = T, rowNames = T, colNames = T)

################## CP LOW ######################
results_Cplow_6h <- results(dds, contrast = c("condition", "CP11_H6", "C_H6"), independentFiltering = T,
                         pAdjustMethod = "BH", alpha = 0.05)
#fdr_results_Cplow_6h <- fdrtool(results_Cplow_6h$stat, statistic = "normal",)
#results_Cplow_6h$padj <- p.adjust(fdr_results_Cplow_6h$pval, method = "BH")
hist(results_Cplow_6h$pvalue)

Cplow6h_res <- data.frame(results_Cplow_6h, geneID = rownames(assay(dds)), row.names = rownames(assay(dds)))
write.xlsx(Ca1h_res, "./RStudio/results/Cp-low-6h_result.xlsx", asTable = T, rowNames = T, colNames = T)

Cplow6h_sig.down <- data.frame(results_Cplow_6h[which(results_Cplow_6h$log2FoldChange < -1.5 & results_Cplow_6h$pvalue < 0.05),])
write.xlsx(Cplow6h_sig.down, "./RStudio/results/Cp-low-6h_sig.down.xlsx", asTable = T, rowNames = T, colNames = T, sheetName = "Sig.down. genes")

Cplow6h_sig.up <- data.frame(results_Cplow_6h[which(results_Cplow_6h$log2FoldChange > 1.5 & results_Cplow_6h$pvalue < 0.05),])
write.xlsx(Cplow6h_sig.up, "./RStudio/results/Cp-low-6h_sig.up.xlsx", asTable = T, rowNames = T, colNames = T, sheetName = "Sig. up genes")
################## CP HIGH #####################
results_Cphigh_6h <- results(dds, contrast = c("condition", "CP51_H6", "C_H6"), independentFiltering = T,
                            pAdjustMethod = "BH", alpha = 0.05)
#fdr_results_Cphigh_6h <- fdrtool(results_Cphigh_6h$stat, statistic = "normal",)
#results_Cphigh_6h$padj <- p.adjust(fdr_results_Cphigh_6h$pval, method = "BH")
hist(results_Cphigh_6h$pvalue)

Cphigh6h_res <- data.frame(results_Cphigh_6h, geneID = rownames(assay(dds)), row.names = rownames(assay(dds)))
write.xlsx(Ca1h_res, "./RStudio/results/Cp-high-6h_result.xlsx", asTable = T, rowNames = T, colNames = T)

Cphigh6h_sig.down <- data.frame(results_Cphigh_6h[which(results_Cphigh_6h$log2FoldChange < -1.5 & results_Cphigh_6h$pvalue < 0.05),])
write.xlsx(Cphigh6h_sig.down, "./RStudio/results/Cp-high-6h_sig.down.xlsx", asTable = T, rowNames = T, colNames = T, sheetName = "Sig.down. genes")

Cphigh6h_sig.up <- data.frame(results_Cphigh_6h[which(results_Cphigh_6h$log2FoldChange > 1.5 & results_Cphigh_6h$pvalue < 0.05),])
write.xlsx(Cphigh6h_sig.up, "./RStudio/results/Cp-high-6h_sig.up.xlsx", asTable = T, rowNames = T, colNames = T, sheetName = "Sig. up genes")

Cp_total_results <- merge(Cplow6h_res, Cphigh6h_res, by = 'geneID', suffixes = c(" (MOI 1:1)", " (MOI 5:1)"), all = T)
write.xlsx(Cp_total_results, "./RStudio/results/Cp_total_results.xlsx", asTable = T, rowNames = T, colNames = T)


##########################
#### IPA symbol table ####
##########################
library(data.table)
genes <- data.frame('geneID' = rownames(assay(dds)))

library(org.Hs.eg.db)
library(DOSE)
library(clusterProfiler)


ensembl <- as.character(as.character(genes$geneID))
symbols <- mapIds(org.Hs.eg.db, keys = ensembl, keytype = 'ENSEMBL', column = 'SYMBOL')
entrezID <- mapIds(org.Hs.eg.db, keys = ensembl, keytype = 'ENSEMBL', column = 'ENTREZID')
IDtable <- data.frame( 'ensemblId' = ensembl, 'geneID' = symbols, 'entrezId' = entrezID)

Ca1h_res$geneID <- mapIds(org.Hs.eg.db, keys = row.names(Ca1h_res), keytype = 'ENSEMBL', column = 'SYMBOL')
Ca6h_res$geneID <- mapIds(org.Hs.eg.db, keys = row.names(Ca6h_res), keytype = 'ENSEMBL', column = 'SYMBOL')
Cplow6h_res$geneID <- mapIds(org.Hs.eg.db, keys = row.names(Cplow6h_res), keytype = 'ENSEMBL', column = 'SYMBOL')
Cphigh6h_res$geneID <- mapIds(org.Hs.eg.db, keys = row.names(Cphigh6h_res), keytype = 'ENSEMBL', column = 'SYMBOL')



##########################
#### IPA C.albi table ####
##########################

IPA_ready_data_Ca <- merge(IDtable, Ca_total_results, by.x = "ensemblId", by.y = 'geneID')
write.xlsx(IPA_ready_data_Ca, "./RStudio/results/IPA_data_c.albi.xlsx")


Ca1h_genelist <- IPA_ready_data_Ca$`log2FoldChange (1h)`
names(Ca1h_genelist) <- as.character(IPA_ready_data_Ca$entrezId)
Ca1h_genelist <- sort(Ca1h_genelist, decreasing = T)
Ca1h_genes <- na.omit(names(Ca1h_genelist)[abs(Ca1h_genelist) > 1.5])


Ca1h_KEGG <- enrichKEGG(Ca1h_genes,
                        organism = 'hsa',
                        keyType = 'kegg',
                        universe =names(Ca1h_genelist),
                        pvalueCutoff = 0.05)

Ca1h_KEGG <- setReadable(Ca1h_KEGG, 'org.Hs.eg.db', keyType = 'ENTREZID')
Ca1h_KEGG_df <- as.data.frame(Ca1h_KEGG@result)
Ca1h_KEGG_df <- Ca1h_KEGG_df[order(Ca1h_KEGG_df$pvalue, decreasing = F),]

write.xlsx(head(Ca1h_KEGG_df,20), "./RStudio/results/Ca1h_KEGG_top_results.xlsx")

Ca1h_GO <- enrichGO(Ca1h_genes,
                    'org.Hs.eg.db',
                    keyType = 'ENTREZID',
                    universe = na.omit(names(Ca1h_genelist)),
                    ont = 'ALL',
                    readable = T)

Ca1h_GO_df <- (as.data.frame(Ca1h_GO@result))
Ca1h_GO_df <- Ca1h_GO_df[order(Ca1h_GO_df$pvalue, decreasing = F),]

write.xlsx(head(Ca1h_GO_df, 20), "./RStudio/results/Ca1h_GO_top_results.xlsx")


Ca6h_genelist <- IPA_ready_data_Ca$`log2FoldChange (6h)`
names(Ca6h_genelist) <- as.character(IPA_ready_data_Ca$entrezId)
Ca6h_genelist <- sort(Ca6h_genelist, decreasing = T)
Ca6h_genes <- na.omit(names(Ca6h_genelist)[abs(Ca6h_genelist) > 1])


Ca6h_KEGG <- enrichKEGG(Ca6h_genes,
                        organism = 'hsa',
                        keyType = 'kegg',
                        universe = na.omit(names(Ca6h_genelist)),
                        pvalueCutoff = 0.05)

Ca6h_KEGG <- setReadable(Ca6h_KEGG, 'org.Hs.eg.db', keyType = 'ENTREZID')
Ca6h_KEGG_df <- (as.data.frame(Ca6h_KEGG@result))
Ca6h_KEGG_df <- Ca6h_KEGG_df[order(Ca6h_KEGG_df$pvalue, decreasing = F),]

write.xlsx(head(Ca6h_KEGG_df, 20), "./RStudio/results/Ca6h_KEGG_top_results.xlsx")

Ca6h_GO <- enrichGO(Ca6h_genes,
                    'org.Hs.eg.db',
                    keyType = 'ENTREZID',
                    universe = na.omit(names(Ca6h_genelist)),
                    ont = 'ALL',
                    readable = T)

Ca6h_GO_df <- (as.data.frame(Ca6h_GO@result))
Ca6h_GO_df <- Ca6h_GO_df[order(Ca6h_GO_df$pvalue, decreasing = F),]

write.xlsx(head(Ca6h_GO_df, 20), "./RStudio/results/Ca6h_GO_top_results.xlsx")

##########################
#### IPA C.para table ####
##########################

IPA_ready_data_Cp <- merge(IDtable, Cp_total_results, by.x = "ensemblId", by.y = 'geneID')
write.xlsx(IPA_ready_data_Cp, "./RStudio/results/IPA_data_c.para.xlsx")


Cplow_genelist <- IPA_ready_data_Cp$`log2FoldChange (MOI 1:1)`
names(Cplow_genelist) <- as.character(IPA_ready_data_Cp$entrezId)
Cplow_genelist <- sort(Cplow_genelist, decreasing = T)
Cplow_genes <- na.omit(names(Cplow_genelist)[abs(Cplow_genelist) > 1.5])


Cplow_KEGG <- enrichKEGG(Cplow_genes,
                        organism = 'hsa',
                        keyType = 'kegg',
                        universe =names(Cplow_genelist),
                        pvalueCutoff = 0.05)

Cplow_KEGG <- setReadable(Cplow_KEGG, 'org.Hs.eg.db', keyType = 'ENTREZID')
Cplow_KEGG_df <- as.data.frame(Cplow_KEGG@result)
Cplow_KEGG_df <- Ca1h_KEGG_df[order(Cplow_KEGG_df$pvalue, decreasing = F),]

write.xlsx(head(Cplow_KEGG_df,20), "./RStudio/results/Cplow_KEGG_top_results.xlsx")

Cplow_GO <- enrichGO(Cplow_genes,
                       'org.Hs.eg.db',
                       keyType = 'ENTREZID',
                       universe = na.omit(names(Cplow_genelist)),
                       ont = 'ALL',
                       readable = T)

Cplow_GO_df <- (as.data.frame(Cplow_GO@result))
Cplow_GO_df <- Cplow_GO_df[order(Cplow_GO_df$pvalue, decreasing = F),]

write.xlsx(head(Cplow_GO_df, 20), "./RStudio/results/Cplow_GO_top_results.xlsx")


CPhigh_genelist <- IPA_ready_data_Cp$`log2FoldChange (MOI 5:1)`
names(CPhigh_genelist) <- as.character(IPA_ready_data_Cp$entrezId)
CPhigh_genelist <- sort(CPhigh_genelist, decreasing = T)
CPhigh_genes <- na.omit(names(CPhigh_genelist)[abs(CPhigh_genelist) > 1])


CPhigh_KEGG <- enrichKEGG(CPhigh_genes,
                        organism = 'hsa',
                        keyType = 'kegg',
                        universe = na.omit(names(CPhigh_genelist)),
                        pvalueCutoff = 0.05)

CPhigh_KEGG <- setReadable(CPhigh_KEGG, 'org.Hs.eg.db', keyType = 'ENTREZID')
CPhigh_KEGG_df <- (as.data.frame(CPhigh_KEGG@result))
CPhigh_KEGG_df <- CPhigh_KEGG_df[order(CPhigh_KEGG_df$pvalue, decreasing = F),]

write.xlsx(head(CPhigh_KEGG_df, 20), "./RStudio/results/CPhigh_KEGG_top_results.xlsx")

CPhigh_GO <- enrichGO(CPhigh_genes,
                    'org.Hs.eg.db',
                    keyType = 'ENTREZID',
                    universe = na.omit(names(CPhigh_genelist)),
                    ont = 'ALL',
                    readable = T)

CPhigh_GO_df <- (as.data.frame(CPhigh_GO@result))
CPhigh_GO_df <- CPhigh_GO_df[order(CPhigh_GO_df$pvalue, decreasing = F),]

write.xlsx(head(CPhigh_GO_df, 20), "./RStudio/results/CPhigh_GO_top_results.xlsx")




library(ggplot2)

Ca1h_KEGG_df$GeneRatio <- sapply(stringr::str_split(Ca1h_KEGG_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
Ca1h_GO_df$GeneRatio <- sapply(stringr::str_split(Ca1h_GO_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))


Ca6h_KEGG_df$GeneRatio <- sapply(stringr::str_split(Ca6h_KEGG_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
Ca6h_GO_df$GeneRatio <- sapply(stringr::str_split(Ca6h_GO_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))


Cplow_KEGG_df$GeneRatio <- sapply(stringr::str_split(Cplow_KEGG_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
Cplow_GO_df$GeneRatio <- sapply(stringr::str_split(Cplow_GO_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))

CPhigh_KEGG_df$GeneRatio <- sapply(stringr::str_split(CPhigh_KEGG_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
CPhigh_GO_df$GeneRatio <- sapply(stringr::str_split(CPhigh_GO_df$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))



library(dplyr)
library(forcats)
library(numbers)
library(plyr)
##########################
### Dotplots C.albi 1h ###
##########################

##### C.albi 1h KEGG #####
attach(Ca1h_KEGG_df)
limits = c(min(Count),max(Count))
modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
breaks = if (max(Count) > 100) {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
} else if (max(Count) > 10 & max(Count) < 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
} else if (max(Count) > 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
} else {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
}
detach(Ca1h_KEGG_df)

png(filename = "./RStudio/pictures/Ca1h_KEGG.png", width = 2260, height = 1580, units = 'px', res = 300)
Ca1h_KEGG_df %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  .[1:min(10, (dim(Ca1h_KEGG_df)[1])),] %>%
  ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
  geom_point(color = 'red') + 
  scale_size(range = c(2,9), limits = limits, breaks =  breaks) + 
  theme_light() + 
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
  labs(size = "Count", y = "")
dev.off()

###### C.albi 1h GO ######
attach(Ca1h_GO_df)
limits = c(min(Count),max(Count))
modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
breaks = if (max(Count) > 100) {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
} else if (max(Count) > 10 & max(Count) < 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
} else if (max(Count) > 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
} else {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
}
detach(Ca1h_GO_df)

png(filename = "./RStudio/pictures/Ca1h_GO.png", width = 2260, height = 1580, units = 'px', res = 300)
Ca1h_GO_df %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  .[1:min(10, (dim(Ca1h_GO_df)[1])),] %>%
  ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
  geom_point(color = 'red') + 
  scale_size(range = c(2,9), limits = limits, breaks =  breaks) + 
  theme_light() + 
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
  labs(size = "Count", y = "")
dev.off()

##### C.albi 6h KEGG #####
attach(Ca6h_KEGG_df)
limits = c(min(Count),max(Count))
modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
breaks = if (max(Count) > 100) {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
} else if (max(Count) > 10 & max(Count) < 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
} else if (max(Count) > 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
} else {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
}
detach(Ca6h_KEGG_df)

png(filename = "./RStudio/pictures/Ca6h_KEGG.png", width = 2260, height = 1580, units = 'px', res = 300)
Ca6h_KEGG_df %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  .[1:min(10, (dim(Ca6h_KEGG_df)[1])),] %>%
  ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
  geom_point(color = 'red') + 
  scale_size(range = c(2,9), limits = limits, breaks =  breaks) + 
  theme_light() + 
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
  labs(size = "Count", y = "")
dev.off()

##### C.albi 6h GO #####
attach(Ca6h_GO_df)
limits = c(min(Count),max(Count))
modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
breaks = if (max(Count) > 100) {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
} else if (max(Count) > 10 & max(Count) < 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
} else if (max(Count) > 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
} else {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
}
detach(Ca6h_GO_df)

png(filename = "./RStudio/pictures/Ca6h_GO.png", width = 2260, height = 1580, units = 'px', res = 300)
Ca6h_GO_df %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  .[1:min(10, (dim(Ca6h_GO_df)[1])),] %>%
  ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
  geom_point(color = 'red') + 
  scale_size(range = c(2,9), limits = limits, breaks =  breaks) + 
  theme_light() + 
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
  labs(size = "Count", y = "")
dev.off()

##### C.para MOI 5:1 KEGG #####
attach(CPhigh_KEGG_df)
limits = c(min(Count),max(Count))
modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
breaks = if (max(Count) > 100) {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
} else if (max(Count) > 10 & max(Count) < 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
} else if (max(Count) > 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
} else {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
}
detach(CPhigh_KEGG_df)

png(filename = "./RStudio/pictures/Cphigh_KEGG.png", width = 2260, height = 1580, units = 'px', res = 300)
CPhigh_KEGG_df %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  .[1:min(10, (dim(CPhigh_KEGG_df)[1])),] %>%
  ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
  geom_point(color = 'red') + 
  scale_size(range = c(2,9), limits = limits, breaks =  breaks) + 
  theme_light() + 
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
  labs(size = "Count", y = "")
dev.off()


##### C.para MOI 5:1 GO #####
attach(CPhigh_GO_df)
limits = c(min(Count),max(Count))
modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
breaks = if (max(Count) > 100) {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
} else if (max(Count) > 10 & max(Count) < 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
} else if (max(Count) > 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
} else {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
}
detach(CPhigh_GO_df)

png(filename = "./RStudio/pictures/Cphigh_GO.png", width = 2260, height = 1580, units = 'px', res = 300)
CPhigh_GO_df %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  .[1:min(10, (dim(CPhigh_GO_df)[1])),] %>%
  ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
  geom_point(color = 'red') + 
  scale_size(range = c(2,9), limits = limits, breaks =  breaks) + 
  theme_light() + 
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
  labs(size = "Count", y = "")
dev.off()

##### C.para MOI 1:1 KEGG #####
attach(Cplow_KEGG_df)
limits = c(min(Count),max(Count))
modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
breaks = if (max(Count) > 100) {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
} else if (max(Count) > 10 & max(Count) < 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
} else if (max(Count) > 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
} else {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
}
detach(Cplow_KEGG_df)

png(filename = "./RStudio/pictures/Cplow_KEGG.png", width = 2260, height = 1580, units = 'px', res = 300)
Cplow_KEGG_df %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  .[1:min(10, (dim(Cplow_KEGG_df)[1])),] %>%
  ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
  geom_point(color = 'red') + 
  scale_size(range = c(2,9), limits = limits, breaks =  breaks) + 
  theme_light() + 
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
  labs(size = "Count", y = "")
dev.off()

##### C.para MOI 1:1 GO #####
attach(Cplow_GO_df)
limits = c(min(Count),max(Count))
modulus = mod(seq(limits[1], limits[2]), ceiling(length(seq(limits[1], limits[2])) / 10))
breaks = if (max(Count) > 100) {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 10, f= ceiling))
} else if (max(Count) > 10 & max(Count) < 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 2, f= ceiling))
} else if (max(Count) > 20){
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 5, f= ceiling))
} else {
  unique(round_any(seq(limits[1],limits[2])[which(modulus == 0)], 1, f= ceiling))
}
detach(Cplow_GO_df)

png(filename = "./RStudio/pictures/Cplow_GO.png", width = 2260, height = 1580, units = 'px', res = 300)
Cplow_GO_df %>%
  mutate(Description = fct_reorder(Description, GeneRatio)) %>%
  .[1:min(10, (dim(Cplow_GO_df)[1])),] %>%
  ggplot(aes(x = GeneRatio, y = Description, size = as.numeric(Count))) + 
  geom_point(color = 'red') + 
  scale_size(range = c(2,9), limits = limits, breaks =  breaks) + 
  theme_light() + 
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) + 
  labs(size = "Count", y = "")
dev.off()




Ca_miR_targets <- read.xlsx("G:/mirna_human/szeged_c.albi/IPA/Ca6h_miR-targets.xlsx")

Ca_miR_target_LFC <- Ca_miR_targets[,2]
Ca_miR_target_names <- names(Ca_miR_target_LFC) <- as.character(Ca_miR_targets[,4])
miR_entrezIDs <- mapIds(org.Hs.eg.db, keys = Ca_miR_target_names, keytype = 'SYMBOL', column = 'ENTREZID')
Ca_miR_targets[,4] <- miR_entrezIDs

Ca_miR_target_LFC <- sort(Ca_miR_target_LFC, decreasing = T)
Ca_miR_target_names <- names(Ca_miR_target_LFC)[abs(Ca_miR_target_LFC) > 1.5]

Ca_miR_target_GO6h <- gseGO(geneList     = Ca_miR_target_LFC,
                            OrgDb        = org.Hs.eg.db,
                            ont          = "ALL",
                            #keytype      = "symbol",
                            nPerm        = 1000,
                            minGSSize    = 10,
                            maxGSSize    = 500,
                            pvalueCutoff = 0.05,
                            verbose      = T)

Ca_miR_target_GO6h <- enrichGO(Ca_miR_target_names,
                               'org.Hs.eg.db',
                               keyType = 'ENTREZID',
                               universe = names(Ca6h_genelist),
                               ont = 'ALL',
                               readable = T)


Ca_miR_target_GO6h <- setReadable(Ca_miR_target_GO6h, 'org.Hs.eg.db', 'ENTREZID')
Ca_miR_target_GO6h_df <- (as.data.frame(Ca_miR_target_GO6h))
Ca_miR_target_GO6h_df <- Ca_miR_target_GO6h_df[order(Ca_miR_target_GO6h_df$pvalue, decreasing = F),]

write.xlsx(head(Ca_miR_target_GO6h_df, 20), "G:/rna_human/szeged_c.albi/Ca_miR_target_GO6h_top_results.xlsx")


dotplot(Ca_miR_target_GO6h,
        showCategory =  20,
        color = 'p.adjust',
        font.size = 16)

library(ViSEAGO)


EntrezGene <- ViSEAGO::EntrezGene2GO()

Ca_miR_2GO <- annotate("9606",
                       EntrezGene)

Ca_miR_GO <- ViSEAGO::create_topGOdata(geneSel = Ca_miR_target_names,
                              allGenes = names(Ca6h_genelist),
                              geneList = NULL,
                              gene2GO = Ca_miR_2GO,
                              ont = 'BP',
                              nodeSize = 5
                              )
Ca_miR_classic <- topGO::runTest(
  Ca_miR_GO,
  algorithm = 'classic',
  statistic = 'fisher'
)
BP_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("Ca_miR_GO","Ca_miR_classic")
  )
)
ViSEAGO::show_table(BP_sResults)
                    #"G:/rna_human/szeged_c.albi/Ca_miR_target_ViSEAGO.xlsx")
ViSEAGO::GOcount(BP_sResults)

myGOs <- ViSEAGO::build_GO_SS(
  gene2GO = Ca_miR_2GO,
  enrich_GO_terms = BP_sResults
)
myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Wang"
)
ViSEAGO::MDSplot(myGOs,
                 file = "G:/rna_human/szeged_c.albi/Ca_miR_target_ViSEAGO.png")

Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=T,
  showGOlabels=T,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=0,
        minClusterSize =5
      )
    )
  ),
  samples.tree=NULL
)
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms")
ViSEAGO::show_table(
  Wang_clusters_wardD2,
  file = "G:/rna_human/szeged_c.albi/Ca_miR_target_ViSEAGO_heatmap.xls")

# calculate semantic similarites between clusters of GO terms
Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
  Wang_clusters_wardD2,
  distance=c("max", "avg","rcmax", "BMA")
)
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOclusters")
  #file = "G:/rna_human/szeged_c.albi/Ca_miR_target_ViSEAGO_clusters.html")

Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
  Wang_clusters_wardD2,
  tree=list(
    distance="BMA",
    aggreg.method="ward.D2"
  )
)
# sisplay the GOClusters heatmap
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOclusters"
)


