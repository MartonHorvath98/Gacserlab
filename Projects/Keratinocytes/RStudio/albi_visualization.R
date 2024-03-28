################################################################################
# 1. Create a PCA plot                                                         #
################################################################################
albi.pca <- make_pca(albi.deseq$rld, group = experiment,
                     labs = c("SC5314(ctrl)","SC5314 vs HaCaT",
                              "WO1 (ctrl)","WO1 vs HaCaT"), 
                     cols = c("salmon","darkred","steelblue","blue"))
ggsave(paste(file.path(folder, date, plots_dir),"albicans_pca.png",sep="/"),
       plot = albi.pca, width = 8, height = 8, units = 'in')
