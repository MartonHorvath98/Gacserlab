library(ggplot2)
#Increased text size for descriptions and facet titles for better readibility
(new_hacat_Calbi_go <- human.GOs_plots$`CTRL vs C.albi (HaCat)` + 
  theme(axis.text.y = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.title.y = element_blank()))

ggsave(filename = "./hacat_vs_calbi_GO.png",
       plot = new_hacat_Calbi_go,
       device = "png",
       width = 16, height = 12, units = 'in')

(new_hacat_Cpara_go <- human.GOs_plots$`CTRL vs C.para (HaCat)` + 
    theme(axis.text.y = element_text(size = 16),
          strip.text.y = element_text(size = 16),
          axis.title.y = element_blank()))

ggsave(filename = "./hacat_vs_cpara_GO.png",
       plot = new_hacat_Cpara_go,
       device = "png",
       width = 16, height = 12, units = 'in')

(new_hpvker_Calbi_go <- human.GOs_plots$`CTRL vs C.albi (HPV-KER)` + 
    theme(axis.text.y = element_text(size = 16),
          strip.text.y = element_text(size = 16),
          axis.title.y = element_blank()))

ggsave(filename = "./hpvker_vs_calbi_GO.png",
       plot = new_hpvker_Calbi_go,
       device = "png",
       width = 16, height = 12, units = 'in')

(new_hpvker_Cpara_go <- human.GOs_plots$`CTRL vs C.para (HPV-KER)` + 
    theme(axis.text.y = element_text(size = 16),
          strip.text.y = element_text(size = 16),
          axis.title.y = element_blank()))

ggsave(filename = "./hpvker_vs_cpara_GO.png",
       plot = new_hpvker_Cpara_go,
       device = "png",
       width = 16, height = 12, units = 'in')

################################################################
# Improved GOSlim set for biologically more relevant answers
library(openxlsx)
library(GOplot)

SC5314_GO_mod <- list(
  BP = read.xlsx(paste(results_dir,"SC5314vsHK_BP_mod.xlsx", sep = "/"), colNames = T),
  MF = read.xlsx(paste(results_dir,"SC5314vsHK_MF_mod.xlsx", sep = "/"), colNames = T),
  CC = read.xlsx(paste(results_dir,"SC5314vsHK_CC.xlsx", sep = "/"), colNames = T)
)

for (i in names (SC5314_GO_mod)){
  SC5314_GO_mod[[i]] <- SC5314_GO_mod[[i]] %>%
    dplyr::mutate(Category = i)
}

SC5314_mod.df <- do.call(rbind.fill, append(SC5314_GO_mod, list(KEGG=SC5314_KEGG_df)))
SC5314_mod.df$Category <- factor(SC5314_mod.df$Category, levels=c("BP", "MF", "CC", "KEGG")) 

SC5314_mod.df <- make_GObase(SC5314_mod.df, SC5314.res$sig_df)
SC5314_mod.df$zscore <- SC5314.df$zscore[match(SC5314_mod.df$id, SC5314.df$id)]

(SC5314_mod.plot <- compound_GOplot(SC5314_mod.df))
ggsave('SC5314_GO_new_plot.png',
       plot = SC5314_mod.plot,device = "png",
       width = 18,height = 12, units = 'in')

##########################
WO1_GO_mod <- list(
  BP = read.xlsx(paste(results_dir,"WO1vsHK_BP_mod.xlsx", sep = "/"), colNames = T),
  MF = read.xlsx(paste(results_dir,"WO1vsHK_MF_mod.xlsx", sep = "/"), colNames = T),
  CC = read.xlsx(paste(results_dir,"WO1vsHK_CC_mod.xlsx", sep = "/"), colNames = T)
)

for (i in names (WO1_GO_mod)){
  WO1_GO_mod[[i]] <- WO1_GO_mod[[i]] %>%
    dplyr::mutate(Category = i)
}

WO1_mod.df <- do.call(rbind.fill, append(WO1_GO_mod, list(KEGG=WO1_KEGG_df)))
WO1_mod.df$Category <- factor(WO1_mod.df$Category, levels=c("BP", "MF", "CC", "KEGG")) 

WO1_mod.df <- make_GObase(WO1_mod.df, WO1.res$sig_df)
WO1_mod.df$zscore <- WO1.df$zscore[match(WO1_mod.df$id, WO1.df$id)]

(WO1_mod.plot <- compound_GOplot(WO1_mod.df))
ggsave('WO1_GO_new_plot.png',
       plot = WO1_mod.plot,device = "png",
       width = 18,height = 12, units = 'in')


###########################################
#modifyin fungi GO plots
library(ggplot2)
library(ggrepel)
(WO1_plotBP_new <- make_FGplot(WO1_GO$BP, "GO") + 
  # geom_text_repel(aes(label = Description), arrow = NULL,
  #                 size = 8, box.padding = unit(.5, "cm"),
  #                 point.padding = unit(1, "cm"), force = 1,
  #                 min.segment.length = 1, segment.color = "transparent",
  #                 direction = "x") + 
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 16)))
ggsave(paste(plots_dir,"WO1vsHK_BP_new_dotplot.png", sep = "/"), WO1_plotBP_new,
       width = 12, height = 10, units = 'in')

(WO1_plotMF_new <- make_FGplot(WO1_GO$MF, "GO") + 
    # geom_text_repel(aes(label = stringr::str_wrap(Description, 50)),
    #                 arrow = NULL,
    #                 size = 8, box.padding = unit(.5, "cm"),
    #                 point.padding = unit(1, "cm"), force = 1,
    #                 min.segment.length = 1, segment.color = "transparent",
    #                 direction = "x") + 
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.text = element_text(size = 16)))
ggsave(paste(plots_dir,"WO1vsHK_MF_new_dotplot.png", sep = "/"), WO1_plotMF_new,
       width = 12, height = 10, units = 'in')

(WO1_plotCC_new <- make_FGplot(WO1_GO$CC, "GO") + 
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.text = element_text(size = 16)))
ggsave(paste(plots_dir,"WO1vsHK_CC_new_dotplot.png", sep = "/"), WO1_plotCC_new,
       width = 12, height = 10, units = 'in')

#SC5314
(SC5314_plotBP_new <- make_FGplot(SC5314_GO$BP, "GO") + 
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.text = element_text(size = 16)))
ggsave(paste(plots_dir,"SC5314vsHK_BP_new_dotplot.png", sep = "/"), SC5314_plotBP_new,
       width = 12, height = 10, units = 'in')

(SC5314_plotMF_new <- make_FGplot(SC5314_GO$MF, "GO") + 
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.text = element_text(size = 16)))
ggsave(paste(plots_dir,"SC5314vsHK_MF_new_dotplot.png", sep = "/"), SC5314_plotMF_new,
       width = 12, height = 10, units = 'in')

(SC5314_plotCC_new <- make_FGplot(SC5314_GO$CC, "GO") + 
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 16),
          axis.text = element_text(size = 16)))
ggsave(paste(plots_dir,"SC5314vsHK_CC_new_dotplot.png", sep = "/"), SC5314_plotCC_new,
       width = 12, height = 10, units = 'in')


