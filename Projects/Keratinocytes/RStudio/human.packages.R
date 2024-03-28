################################################################################
# Load required packages                                                       #
################################################################################

# 1.) Loader packages ##########################################################
if (!require("BiocManager")) install.packages("BiocManager")
suppressPackageStartupMessages(library(BiocManager))
# Provides tools for managing Bioconductor packages
if (!require("openxlsx")) install.packages("openxlsx")
suppressPackageStartupMessages(library(openxlsx)) 
# Simplifies the creation, writing, and editing of '.xlsx' files
if (!require("devtools")) install.packages("devtools")
suppressPackageStartupMessages(library(devtools)) 
# Collection of package development tools.

# 2.) Utility packages #########################################################

if (!require("crayon")) install.packages("crayon")
suppressPackageStartupMessages(library(crayon))
# Provides a set of functions to style console output
if (!require("dplyr")) install.packages("dplyr")
suppressPackageStartupMessages(library(dplyr))
# tool for working with data frame like objects
if (!require("tidyr")) install.packages("tidyr")
suppressPackageStartupMessages(library(tidyr))
# contains tools for changing the shape (pivoting) and hierarchy (nesting and 
# 'unnesting') of a dataset
if (!require("tibble")) install.packages("tibble")
suppressPackageStartupMessages(library(tibble))
# provides a 'tbl_df' class that offers better formatting than the default data 
# frame
if (!require("stringr")) install.packages("stringr")
suppressPackageStartupMessages(library(stringr))
# provides a consistent, simple and easy to use set of wrappers around the
# 'stringi' package 
if (!require("BiocParallel")) install.packages("BiocParallel")
suppressPackageStartupMessages(library(BiocParallel)) 
# Provides a framework for parallel evaluation of R expressions using various
# back-ends


# 3.) Annotation and aligner packages ##########################################

if (!require("Rsamtools")) install.packages("Rsamtools")
suppressPackageStartupMessages(library(Rsamtools))
# Provides an interface to the 'samtools', 'bcftools', and 'tabix' utilities
if (!require("GenomicAlignments")) install.packages("GenomicAlignments")
suppressPackageStartupMessages(library(GenomicAlignments))
# Provides efficient containers for storing and manipulating short genomic
# alignments
if (!require("GenomicFeatures")) install.packages("GenomicFeatures")
suppressPackageStartupMessages(library(GenomicFeatures))
# Provides tools for making and manipulating transcript centric annotations
if (!require("GenomicRanges")) install.packages("GenomicRanges")
suppressPackageStartupMessages(library(GenomicRanges))
# Provides efficient containers for storing and manipulating genomic intervals
if (!require("GenomeInfoDb")) install.packages("GenomeInfoDb")
suppressPackageStartupMessages(library(GenomeInfoDb))
# Provides infrastructure for manipulating and summarizing genome-wide data

# 4.) Differential expression analysis packages ################################

if (!require("DESeq2")) install.packages("DESeq2")
suppressPackageStartupMessages(library(DESeq2)) 
# provides methods to test for differential expression by estimating variance-mean
# dependence in count data from high-throughput sequencing assays and test for
# differential expression based on a model using the negative binomial distribution
if (!require("edgeR")) BiocManager::install("edgeR", update = F)
suppressPackageStartupMessages(library(edgeR))
# provides functions for analyzing read counts from RNA-Seq experiments
if (!require("fdrtool")) install.packages("fdrtool", repos = c("http://R-Forge.R-project.org"))
suppressPackageStartupMessages(library(fdrtool))
# provides a method to estimate the proportion of true null hypotheses among a 
# set of hypotheses rejected at a given significance level
if (!require("AnnotationDbi")) install.packages("AnnotationDbi")
suppressPackageStartupMessages(library(AnnotationDbi))
# provides a set of database interfaces for common annotation packages
if (!require("org.Hs.eg.db")) BiocManager::install("org.Hs.eg.db", update = F)
suppressPackageStartupMessages(library(org.Hs.eg.db))
# provides mappings between Entrez Gene identifiers and other gene identifiers
if (!require("genefilter")) BiocManager::install("genefilter",version = "3.8", update = F)
suppressPackageStartupMessages(library(genefilter))
# provides methods for filtering genes based on different criteria 

# 5.) Pathway analysis packages ################################################
if (!require("dendextend")) install.packages("dendextend")
suppressPackageStartupMessages(library(dendextend))
# provides a set of functions to manipulate dendrograms
if (!require("dendsort")) install.packages("dendsort")
suppressPackageStartupMessages(library(dendsort))
# provides a set of functions to sort dendrograms
if (!require("cluster")) install.packages("cluster")
suppressPackageStartupMessages(library(cluster))
# provides a set of functions to cluster data

# 6.) Visualization packages ###################################################

if (!require("ggplot2")) install.packages("ggplot2")
suppressPackageStartupMessages(library(ggplot2))
# provides a system for 'declaratively' creating graphics, based on "The Grammar
# of Graphics"
if (!require("ggbiplot")) install.packages("ggbiplot")
suppressPackageStartupMessages(library(ggbiplot))
# provides a ggplot2 extension for creating biplots
if (!require("ggrepel")) install.packages("ggrepel")
suppressPackageStartupMessages(library(ggrepel))
# provides geom for ggplot2 to repel overlapping text labels
if (!require("pheatmap")) install.packages("pheatmap")
suppressPackageStartupMessages(library(pheatmap))
# provides a function to create pretty heatmaps
if (!require("RColorBrewer")) install.packages("RColorBrewer")
suppressPackageStartupMessages(library(RColorBrewer))
# provides a set of color palettes for creating pretty graphics
if (!require("ComplexHeatmap")) install.packages("ComplexHeatmap")
suppressPackageStartupMessages(library(ComplexHeatmap))
# provides a function to create complex heatmaps with multiple annotations
if (!require("circlize")) install.packages("circlize")
suppressPackageStartupMessages(library(circlize))
# provides a function to create circular plots for genomic data
if (!require("factoextra")) install.packages("factoextra")
suppressPackageStartupMessages(library(factoextra))
# provides a set of functions to extract and visualize the results of multivariate
# data analyses
if (!require("scales")) install.packages("scales")
suppressPackageStartupMessages(library(scales))
# provides functions to scale data
if (!require("viridis")) install.packages("viridis")
suppressPackageStartupMessages(library(viridis))
# provides a set of color palettes for creating pretty graphics
if (!require("VennDiagram")) install.packages("VennDiagram")
suppressPackageStartupMessages(library(VennDiagram))
# provides a function to create Venn diagrams
if (!require("ggvenn")) install.packages("ggvenn")
suppressPackageStartupMessages(library(ggvenn))
# provides a function to create Venn diagrams
if (!require("GOplot")) install.packages("GOplot", repos = c("http://R-Forge.R-project.org"))
suppressPackageStartupMessages(library(GOplot))
# provides a function to create gene ontology plots