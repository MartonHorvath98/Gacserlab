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

# 3.) Statistical packages #####################################################

if (!require("forcats")) install.packages("forcats")
suppressPackageStartupMessages(library(forcats))
# provides a suite of tools that solve common problems with factors
if (!require("numbers")) install.packages("numbers")
suppressPackageStartupMessages(library(numbers))
# provides a consistent way to work with numbers in R
if (!require("countToFPKM")) install.packages("countToFPKM")
suppressPackageStartupMessages(library(countToFPKM))
# provides a function to convert read counts to FPKM values

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


# 5.) Visualization packages ###################################################

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
if (!require("FactomineR")) install.packages("FactomineR")
suppressPackageStartupMessages(library(FactomineR))
# provides a set of functions for multivariate exploratory data analysis


