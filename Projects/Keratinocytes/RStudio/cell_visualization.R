################################################################################
# 1. Create a PCA plot for each cell type                                      #
################################################################################
# 1.1 HaCat PCA
hacat.pca <- make_pca(rld = hacat.deseq$rld, 
                      group = hacat.deseq$rld$treatment,
                      labs = c("Control", "C.albicans (SC5314)",
                               "C.parapsilosis (CLIB214)"),
                      cols = c("C.para" = "salmon",
                               "C.albi" = "darkred",
                               "ctrl" = "turquoise"))
hacat.pca <- hacat.pca + 
  theme(legend.position = "bottom")
ggsave("./Results/Rstudio/pictures/hacat_pca.png", plot = hacat.pca,
       width = 10, height = 10, units = 'in')

# 1.2 HPV-KER PCA
hpvker.pca <- make_pca(rld = hpvker.deseq$rld, 
                       group = hpvker.deseq$rld$treatment,
                       labs = c("Control", "C.albicans (SC5314)",
                                "C.parapsilosis (CLIB214)"),
                       cols = c("C.para" = "salmon",
                                "C.albi" = "darkred",
                                "ctrl" = "turquoise"))
hpvker.pca <- hpvker.pca + 
  theme(legend.position = "bottom")
ggsave("./Results/Rstudio/pictures/hpvker_pca.png", plot = hpvker.pca,
       width = 10, height = 10, units = 'in')

# 1.3 Total PCA - compare the two cell types
cells.pca <- make_pca(rld = human.deseq$rld, 
                      group = human.deseq$rld$condition,
                      labs = c("C.albicans (HaCat)", "C.parapsilosis (HaCat)", "Control (HaCat)",
                               "C.albicans (HPV-KER)", "C.parapsilosis (HPV-KER)", "Control (HPV-KER)"),
                      cols = c("HaCat C.para" = "salmon",
                               "HaCat C.albi" = "darkred",
                               "HaCat ctrl" = "orange",
                               "HPVKER C.para" = "blue",
                               "HPVKER C.albi" = "purple",
                               "HPVKER ctrl" = "turquoise"))
ggsave("./Results/Rstudio/pictures/total_pca.png", plot = cells.pca,
       width = 12, height = 12, units = 'in')

