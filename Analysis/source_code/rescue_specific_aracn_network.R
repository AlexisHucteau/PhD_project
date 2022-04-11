library(dplyr)

R_NR_msviper[["mrs"]][["regulon"]]

regulon_aml_patient <- lapply(names(R_NR_msviper[["mrs"]][["regulon"]]), function(TF){
  nb <- length(R_NR_msviper[["mrs"]][["regulon"]][[TF]]$tfmode)
  data.frame(TFs = rep(TF, nb),
             target = names(R_NR_msviper[["mrs"]][["regulon"]][[TF]]$tfmode),
             mor = R_NR_msviper[["mrs"]][["regulon"]][[TF]]$tfmode)
}) %>% data.table::rbindlist()

write.csv(regulon_aml_patient, "GitHub/PhD_project/Analysis/In_silico/Results_koichi_multiomic/regulons.csv")

regulon_aml_patient_to_cytoscape <- regulon_aml_patient
regulon_aml_patient_to_cytoscape$shared_name <- paste0(regulon_aml_patient_to_cytoscape$TFs, " (interacts with) ",regulon_aml_patient_to_cytoscape$target)
write.csv(regulon_aml_patient_to_cytoscape, "GitHub/PhD_project/Analysis/In_silico/Results_koichi_multiomic/regulons2cytos.csv")
