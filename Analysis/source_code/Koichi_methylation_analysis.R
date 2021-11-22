library(dplyr)
library(ChAMP)

load("~/DATA/BMIQ_Koichi.RData")
Clinical_patient_data <- read.csv("~/GitHub/Koichi_gene_expression_git/Koichi_gene_expression_analyses/DATA/Clinical_patient_data.csv")
Clinical_patient_data <- dplyr::filter(Clinical_patient_data, Baseline_Sample %in% colnames(BMIQ_Koichi) | Post_treatment_sample %in% colnames(BMIQ_Koichi))[,c(5,7, 8)] %>% unique()
colnames(Clinical_patient_data)[c(2,3)] <- rep("Sample", 2)
Clinical_patient_data <- rbind(Clinical_patient_data[,c(1,2)], Clinical_patient_data[,c(1,3)])
Clinical_patient_data$Best_response[c(1:64)] <- paste(Clinical_patient_data$Best_response[c(1:64)], "B", sep = ".")
Clinical_patient_data$Best_response[c(65:128)] <- paste(Clinical_patient_data$Best_response[c(65:128)], "REL", sep = ".")
Phenotype <- Clinical_patient_data$Best_response
names(Phenotype) <- Clinical_patient_data$Sample

Factor_Koichi_response_methylation <- Phenotype[colnames(BMIQ_Koichi)]

Factor_Koichi_response_methylation <- ifelse(Factor_Koichi_response_methylation %in% c("CR.B", "CRi.B"), "R", 
                                             ifelse(Factor_Koichi_response_methylation %in% c("PD.B", "SD.B"), "NR", "Other"))

Focus_samples <- Factor_Koichi_response_methylation %in% c("R", "NR")

DMR_Res_vs_NonRes_Koichi <- champ.DMR(as.matrix(BMIQ_Koichi[,Focus_samples]), 
                                      pheno = Factor_Koichi_response_methylation[Focus_samples], 
                                      cores = 6, 
                                      arraytype = "EPIC")

DMP_Res_vs_NonRes_Koichi <- champ.DMP(as.matrix(BMIQ_Koichi[,Focus_samples]), 
                                      pheno = Factor_Koichi_response_methylation[Focus_samples], 
                                      arraytype = "EPIC")

Block_Res_vs_NonRes_Koichi <- champ.Block(beta = as.matrix(BMIQ_Koichi[,Focus_samples]), 
            pheno = Factor_Koichi_response_methylation[Focus_samples], 
            arraytype = "EPIC", cores = 6)

DNAmet_GSEA <- champ.GSEA(beta = as.matrix(BMIQ_Koichi[,Focus_samples]), 
           DMP = DMP_Res_vs_NonRes_Koichi, 
           DMR = DMR_Res_vs_NonRes_Koichi, 
           pheno = Factor_Koichi_response_methylation[Focus_samples],
           arraytype = "EPIC")


DMP.GUI(DMP = DMP_Res_vs_NonRes_Koichi$R_to_NR, 
        beta = as.matrix(BMIQ_Koichi[,Focus_samples]), 
        pheno = Factor_Koichi_response_methylation[Focus_samples]
)

DMR.GUI(DMR = DMR_Res_vs_NonRes_Koichi, 
        beta = as.matrix(BMIQ_Koichi[,Focus_samples]), 
        pheno = Factor_Koichi_response_methylation[Focus_samples], 
        arraytype = "EPIC"
)

Block.GUI(Block = Block_Res_vs_NonRes_Koichi, 
          beta = as.matrix(BMIQ_Koichi[,Focus_samples]), 
          pheno = Factor_Koichi_response_methylation[Focus_samples], 
          arraytype = "EPIC")
