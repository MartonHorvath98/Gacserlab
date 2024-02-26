library(msigdbr)
library(dplyr)
EMT_candida <- shared_set$geneID
EMT_hallmark <- msigdbr_df %>%
  dplyr::filter(
    gs_cat == "H", # only canonical representations (compiled by experts)
    stringr::str_detect(gs_name, "EPITHELIAL_MESENCHYMAL") # KEGG pathways
  ) %>%
  dplyr::pull(gene_symbol)

pEMT <- read.csv('./data/pEMT_genes.txt',header = F)
pEMT_genes <- pEMT$V1

library(VennDiagram)
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(
    x,
    filename = NULL,
    ...
  )
  grid.draw(venn_object)
}

EMT_venn_base <- list(
  'Candida-driven EMT' = EMT_candida,
  'EMT hallmark' = EMT_hallmark,
  'partial-EMT' = pEMT_genes
) 

library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

display_venn(EMT_venn_base,
             # Circles
             lwd = 2,
             lty = 'solid',
             fill = myCol,
             
             # Numbers
             cex = 1.5,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 2,
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 0),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)
