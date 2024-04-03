##############################################
# 12.) Compare HSC2 & OKF6 cell responses    #
##############################################

# [1] "Angiogenesis"         "Cancer.Metabolism"    "ECM.Layers"           "ECM.Remodeling"      
# [5] "EMT"                  "Hypoxia"              "Metastasis"           "Transcription.Factor"
# [9] "Tumor.Growth"         "Tumor.Invasion" 
#Angiogenesis
Angiogenesis_shared <- make_quad_venn(Angiogenesis_vennbase, OKF6_Angiogenesis_vennbase,
                                      c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                        "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"),
                                      title = "Angiogenesis", dir = plots_dir)
Angiogenesis_shared$df <- make_venn_df(Angiogenesis_shared$df,
                                       c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                         "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"))

openxlsx::write.xlsx(Angiogenesis_shared$df, paste(results_dir, "Angiogenesis_shared.xlsx",sep = "/"))  

# png(paste(plots_dir,"Angiogenesis_shared.png", sep = "/"),
#     width = 16, height = 10, units = 'in', res = 300)
# Angiogenesis_shared$plot
# dev.off()

#Cancer.Metabolism
Cancer_metabolism_shared <- make_quad_venn(Cancer_metabolism_vennbase, OKF6_metabolism_vennbase,
                                           c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                             "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"),
                                           title = "Metabolism", dir = plots_dir)

Cancer_metabolism_shared$df <- make_venn_df(Cancer_metabolism_shared$df,
                                            c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                              "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"))
openxlsx::write.xlsx(Cancer_metabolism_shared$df, paste(results_dir, "Cancer_metabolism_shared.xlsx",sep = "/"))  


#ECM.Remodeling
ECM_remodeling_shared <- make_quad_venn(ECM_remodeling_vennbase, OKF6_ECM_remodeling_vennbase,
                                        c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                          "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"),
                                        title = "ECM_remodeling", dir = plots_dir)
ECM_remodeling_shared$df <- make_venn_df(ECM_remodeling_shared$df,
                                         c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                           "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"))
openxlsx::write.xlsx(ECM_remodeling_shared$df, paste(results_dir, "ECM_remodeling_shared.xlsx",sep = "/"))  

#EMT
EMT_shared <- make_quad_venn(EMT_vennbase, OKF6_EMT_vennbase,
                             c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                               "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"),
                             title = "EMT", dir = plots_dir)

EMT_shared$df <- make_venn_df(EMT_shared$df,
                              c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"))
openxlsx::write.xlsx(EMT_shared$df, paste(results_dir, "EMT_shared.xlsx",sep = "/"))  

#Hypoxia
Hypoxia_shared <- make_quad_venn(Hypoxia_vennbase, OKF6_Hypoxia_vennbase,
                                 c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                   "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"),
                                 title = "Hypoxia", dir = plots_dir)
Hypoxia_shared$df <- make_venn_df(Hypoxia_shared$df,
                                  c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                    "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"))
openxlsx::write.xlsx(Hypoxia_shared$df, paste(results_dir, "Hypoxia_shared.xlsx",sep = "/"))  
#Metastasis
Metastasis_shared <- make_quad_venn(Metastasis_vennbase, OKF6_Metastasis_vennbase,
                                    c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                      "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"),
                                    title = "Metastasis", dir = plots_dir)
Metastasis_shared$df <- make_venn_df(Metastasis_shared$df,
                                     c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                       "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"))
openxlsx::write.xlsx(Metastasis_shared$df, paste(results_dir, "Metastasis_shared.xlsx",sep = "/"))  

#Transcription.Factor
Cancer_TFs_shared <- make_quad_venn(Transcription_vennbase, OKF6_Transcription_vennbase,
                                    c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                      "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"),
                                    title = "Cancer_TFs", dir = plots_dir)
Cancer_TFs_shared$df <- make_venn_df(Cancer_TFs_shared$df,
                                     c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                       "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"))
openxlsx::write.xlsx(Cancer_TFs_shared$df, paste(results_dir, "Cancer_TFs_shared.xlsx",sep = "/")) 
#Tumor.Growth
Tumor_Growth_shared <- make_quad_venn(Tumor_growth_vennbase, OKF6_Tumor_growth_vennbase,
                                      c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                        "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"),
                                      title = "Tumor_growth", dir = plots_dir)

Tumor_Growth_shared$df <- make_venn_df(Tumor_Growth_shared$df,
                                       c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                         "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"))
openxlsx::write.xlsx(Tumor_Growth_shared$df, paste(results_dir, "Tumor_growth_shared.xlsx",sep = "/"))
#Tumor.Invasion
Tumor_Invasion_shared <- make_quad_venn(Tumor_invasion_venn, OKF6_Tumor_invasion_venn,
                                        c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                          "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"),
                                        title = "Tumor_invasion", dir = plots_dir)

Tumor_Invasion_shared$df <- make_venn_df(Tumor_Invasion_shared$df,
                                         c("C.albi 1:1 (HSC2)", "C.para 5:1 (HSC2)",
                                           "C.albi 1:1 (OKF6)", "C.para 5:1 (OKF6)"))
openxlsx::write.xlsx(Tumor_Invasion_shared$df, paste(results_dir, "Tumor_Invasion_shared.xlsx",sep = "/"))
