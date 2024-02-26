# "devtools" package to be able to download packages from different sources
if (!("devtools" %in% installed.packages())) {
  BiocManager::install("devtools", update = FALSE)
}
# "stringr" to manipulate strings
if (!("stringr" %in% installed.packages())) {
  BiocManager::install("stringr", update = FALSE)
}
# "dplyr" for the %>% pipe operator and data wrangling
if (!("dplyr" %in% installed.packages())) {
  BiocManager::install("dplyr", update = FALSE)
}
# "genomicFeatures" package to load reference transcriptome
if (!("GenomicFeatures" %in% installed.packages())) {
  BiocManager::install("GenomicFeatures", update = FALSE)
}
# "Rsamtools" package to load and work with .bam files
if (!("Rsamtools" %in% installed.packages())) {
  BiocManager::install("Rsamtools", update = FALSE)
}
# "GenomeInfoDb" package to work with annotation data bases
if (!("GenomeInfoDb" %in% installed.packages())) {
  BiocManager::install("GenomeInfoDb", update = FALSE)
}
# "TCGAbiolinks" package to work with annotation data bases
if (!("TCGAbiolinks" %in% installed.packages())) {
  BiocManager::install("TCGAbiolinks", update = FALSE)
}
# "BiocParallel" package to enable parallel computing
if (!("BiocParallel" %in% installed.packages())) {
  BiocManager::install("BiocParallel", update = FALSE)
}
# "GenomicAlignments" package to calculate readcounts
if (!("GenomicAlignments" %in% installed.packages())) {
  BiocManager::install("GenomicAlignments", update = FALSE)
}
#summarizeOverlaps(function) - parameters: 
#   - features                            -> reference exons we align the reads
#   - reads                               -> the BamFileList object (our reads)
#   - mode [= "union"]                    -> (default) reads are counted when  
#                                            align with exactly one exon
#   - ingore.strand [= FALSE]             -> considering strand direction
#   - singleEnd [= FALSE]                 -> would be TRUE (default), if the 
#                                            sequencing is single-ended
# - fragments                             -> considering every read (singletons, 
#                                            etc.), during counting
# - preprocessed.reads [='invertStrand']  -> strandspecificity 
#                                            (= leading revers-strand, RF)

# "DESeq2" package for differential expression analysis
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}
# "ggplot2" package for visualisation of results
if (!("ggplot2" %in% installed.packages())) {
  BiocManager::install("ggplot2", update = FALSE)
}
# "ggbiplot" package for visualisation of pca
if (!("ggbiplot" %in% installed.packages())) {
  devtools::install_github("vqv/ggbiplot", upgrade = FALSE)
}
# "ggrepel" package for labels on plots
if (!("ggrepel" %in% installed.packages())) {
  BiocManager::install("ggrepel", update = FALSE)
}
# "fdrtool" package to positive error rate adjustment
if (!("fdrtool" %in% installed.packages())) {
  BiocManager::install("fdrtool", update = FALSE)
}
# "apeglm" package for lfc normalization
if (!("ashr" %in% installed.packages())) {
  BiocManager::install("ashr", update = FALSE)
}
# "org.Hs.eg.db" package for annotation of genes
if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", update = FALSE)
}
# "openxlsx" package for saving data tables as excel files
if (!("openxlsx" %in% installed.packages())) {
  BiocManager::install("openxlsx", update = FALSE)
}
# "clusterProfiler" package for overexpression and gene set enrichment analyses
if (!("clusterProfiler" %in% installed.packages())) {
  BiocManager::install("clusterProfiler", update = FALSE)
}
# "msigdbr" package for pathway and hallmark gene set data
if (!("msigdbr" %in% installed.packages())) {
  BiocManager::install("msigdbr", update = FALSE)
}
# "numbers" & "forcats" package for mathematical formulas
if (!("numbers" %in% installed.packages())) {
  BiocManager::install("numbers", update = FALSE)
}
if (!("forcats" %in% installed.packages())) {
  BiocManager::install("forcats", update = FALSE)
}
# "viridis" package for the viridis colour scheme
if (!("viridis" %in% installed.packages())) {
  BiocManager::install("viridis", update = FALSE)
}
# "ComplexUpset" package to make upset plot of shared genes
if (!("ComplexUpset" %in% installed.packages())) {
  BiocManager::install("ComplexUpset", update = FALSE)
}
# "VennDiagram" package to make venn diagrams
if (!("VennDiagram" %in% installed.packages())) {
  BiocManager::install("VennDiagram", update = FALSE)
}
# "GOplot" package to make venn diagrams
if (!("GOplot" %in% installed.packages())) {
  BiocManager::install("GOplot", update = FALSE)
}
# "GEOquery" package to make venn diagrams
if (!("GEOquery" %in% installed.packages())) {
  BiocManager::install("GEOquery", update = FALSE)
}
# dependency of WGCNA that has trouble loading sometimes
if (!("impute" %in% installed.packages())) {
  BiocManager::install("impute", update = F)
}
# "WGCNA" package to make weighted co-variance analysis for DE clustering
if (!("WGCNA" %in% installed.packages())) {
  BiocManager::install("WGCNA", update = F)
}

if (!("ggforce" %in% installed.packages())) {
  BiocManager::install("ggforce", update = F)
}

if (!("ComplexHeatmap" %in% installed.packages())) {
  BiocManager::install("ComplexHeatmap", update = F)
}

if(!"RCy3" %in% installed.packages()){
  BiocManager::install("RCy3",update = F)
}

if(!"gProfileR" %in% installed.packages()){
  BiocManager::install("gProfileR",update = F)
}

if (!"ggforest" %in% installed.packages()){
  BiocManager::install("ggforest", update = F)
}

