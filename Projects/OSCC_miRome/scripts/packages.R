# GenomicFeatures
if (!("GenomicFeatures" %in% installed.packages())) {
  BiocManager::install("GenomicFeatures", update = FALSE)
}
library(GenomicFeatures)
# Rsamtools
if (!("Rsamtools" %in% installed.packages())) {
  BiocManager::install("Rsamtools", update = FALSE)
}
library(Rsamtools)
# GenomeInfoDb
if (!("GenomeInfoDb" %in% installed.packages())) {
  BiocManager::install("GenomeInfoDb", update = FALSE)
}
library(GenomeInfoDb)
# GenomicAlignments
if (!("GenomicAlignments" %in% installed.packages())) {
  BiocManager::install("GenomicAlignments", update = FALSE)
}
library(GenomicAlignments)
# BiocParallel
if (!("BiocParallel" %in% installed.packages())) {
  BiocManager::install("BiocParallel", update = FALSE)
}
library(BiocParallel)
# stringr
if (!("stringr" %in% installed.packages())) {
  BiocManager::install("stringr", update = FALSE)
}
library(stringr)
# dplyr
if (!("dplyr" %in% installed.packages())) {
  BiocManager::install("dplyr", update = FALSE)
}
library(dplyr)
# openxlsx
if (!("openxlsx" %in% installed.packages())) {
  BiocManager::install("openxlsx", update = FALSE)
}
library(openxlsx)
# DESeq2
if (!("DESeq2" %in% installed.packages())) {
  BiocManager::install("DESeq2", update = FALSE)
}
library(DESeq2)
# edgeR
if (!("edgeR" %in% installed.packages())) {
  BiocManager::install("edgeR", update = FALSE)
}
library(edgeR)
# ggplot2
if (!("ggplot2" %in% installed.packages())) {
  BiocManager::install("ggplot2", update = FALSE)
}
library(ggplot2)
# scales
if (!("scales" %in% installed.packages())) {
  BiocManager::install("scales", update = FALSE)
}
library(scales)
# ggbiplot
if (!("ggbiplot" %in% installed.packages())) {
  devtools::install_github("vqv/ggbiplot")
}
library(ggbiplot)
# ggrepel
if (!("ggrepel" %in% installed.packages())) {
  BiocManager::install("ggrepel", update = FALSE)
}
library(ggrepel)
# fdrtool
if (!("fdrtool" %in% installed.packages())) {
  BiocManager::install("fdrtool", update = FALSE)
}
library(fdrtool)
# ashr
if (!("ashr" %in% installed.packages())) {
  BiocManager::install("ashr", update = FALSE)
}
library(ashr)
# org.Hs.eg.db
if (!("org.Hs.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Hs.eg.db", force = T)
}
library(org.Hs.eg.db)
# clusterProfiler
if (!("clusterProfiler" %in% installed.packages())) {
  BiocManager::install("clusterProfiler", force = T)
}
library(clusterProfiler)
# KEGGgraph
if (!("KEGGgraph" %in% installed.packages())) {
  BiocManager::install("KEGGgraph", force = T)
}
library(KEGGgraph)
# Rgraphviz
if (!("Rgraphviz" %in% installed.packages())) {
  BiocManager::install("Rgraphviz", force = T)
}
library(Rgraphviz)
# createKEGGdb
if (!("createKEGGdb" %in% installed.packages())) {
  devtools::install_github("YuLab-SMU/createKEGGdb")
}
library(createKEGGdb)
# forcats
if (!("forcats" %in% installed.packages())) {
  BiocManager::install("forcats", force = T)
}
library(forcats)
# reshape
if (!("reshape" %in% installed.packages())) {
  BiocManager::install("reshape", force = T)
}
library(reshape)
# numbers
if (!("numbers" %in% installed.packages())) {
  BiocManager::install("numbers", force = T)
}
library(numbers)
# plyr
if (!("plyr" %in% installed.packages())) {
  BiocManager::install("plyr", update = FALSE)
}
library(plyr)
# tibble
if (!("tibble" %in% installed.packages())) {
  BiocManager::install("tibble", update = FALSE)
}
library(tibble)
# VennDiagram
if (!("VennDiagram" %in% installed.packages())) {
  BiocManager::install("VennDiagram", update = FALSE)
}
library(VennDiagram)
# GOplot
if (!("GOplot" %in% installed.packages())) {
  BiocManager::install("GOplot", update = FALSE)
}
library(GOplot)
# ComplexUpset
if (!("ComplexUpset" %in% installed.packages())) {
  BiocManager::install("ComplexUpset", update = FALSE)
}
library(ComplexUpset)
# conicfit
if (!("conicfit" %in% installed.packages())) {
  BiocManager::install("conicfit", update = FALSE)
}
library(conicfit)
# metan
if (!("pheatmap" %in% installed.packages())) {
  BiocManager::install("pheatmap", update = FALSE)
}
library(pheatmap)
# purrr
if (!("purrr" %in% installed.packages())) {
  BiocManager::install("purrr", update = FALSE)
}
library(purrr)
# rrvgo
if (!("rrvgo" %in% installed.packages())) {
  BiocManager::install("rrvgo", update = FALSE)
}
library(rrvgo)
# GOSemSim
if (!("GOSemSim" %in% installed.packages())) {
  BiocManager::install("GOSemSim", update = FALSE)
}
library(GOSemSim)
#igraph
if (!("igraph" %in% installed.packages())) {
  BiocManager::install("igraph", update = FALSE)
}
library(igraph)
#arrow
if (!("arrow" %in% installed.packages())) {
  BiocManager::install("arrow", update = FALSE)
}
library(arrow)
#circlize
if (!("circlize" %in% installed.packages())) {
  BiocManager::install("circlize", update = FALSE)
}
library(circlize)
