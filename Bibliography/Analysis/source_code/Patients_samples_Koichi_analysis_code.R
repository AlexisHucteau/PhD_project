library(FactoMineR)
library(igraph)
library(limma)
library(viper)
library(stringr)
library(data.table)
library(dplyr)

"%ni%" <- Negate("%in%")
RNAseq_diff_gene_expression_analysis <- list()
RNAseq <- read.csv("~/GitHub/Koichi_gene_expression_git/Koichi_gene_expression_analyses/DATA/RNAseq_parsed.csv", row.names = 1, header = T, check.names = F)
Clinical_patient_data <- read.csv("~/GitHub/Koichi_gene_expression_git/Koichi_gene_expression_analyses/DATA/Clinical_patient_data.csv")


PCA_rnaseq <- PCA(t(RNAseq))

Make_factor <- function(Samplesheet = Clinical_patient_data, 
                        Samples_names, 
                        Mutations_to_ignore = 0, 
                        Clinical_outcome_A, 
                        Clinical_name_A,
                        Clinical_outcome_B,  
                        Clinical_name_B,
                        Clinical_outcome_C, 
                        Clinical_name_C){
  # Function made for Clinical_patient_data 
  # Create a factor that can be used for Differential_analysis function
  # Samplesheet = Clinical_patient_data
  # Mutations_to_ignore: A vector of mutations that have to be taken into account (type 0 no mutations to ignore)
  # Clinical_outcome_A: A vector of best response corresponding to the phenotype A
  # Clinical_outcome_B: A vector of best response corresponding to the phenotype B
  # Clinical_outcome_C: A vector of best response corresponding to the phenotype C
  # Baseline_sample: A logical variable indicating whether Baseline samples are taken or not
  # Relapse_sample: A logical variable indicating whether Relapse samples are taken or not
  
  # Phenotype_A: The name of the first phenotype that have to be compared to
  # Phenotype _B: The name of the second phenotype that have to be compared to
  # Clinical_outcome_comparison: A logical variable indicating whether clinical outcome are taken into account
  # Baseline: 
  # Relapse: A logical variable indicating whether Relapse samples are taken or not
  if(typeof(Mutations_to_ignore) != "double"){
    Mutations_samples <- Samplesheet[which(duplicated(str_split(Samplesheet$mutations, pattern=","), Mutations_to_ignore)),] %>% 
      c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>% 
      na.omit()
    Mutations_factor <- factor(ifelse(Samples_names %in% Mutations_samples, "Mut", "WT"))
  }else{
    Mutations_factor <- factor(rep("", length(Samples_names)))
  }
  
  Clinical_outcome_A <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_A),] %>% 
    c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>% 
    na.omit()
  Clinical_outcome_B <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_B),] %>% 
    c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>% 
    na.omit()
  Clinical_outcome_C <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_C),] %>% 
    c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>% 
    na.omit()
  
  Clinical_outcome <- factor(ifelse(Samples_names %in% Clinical_outcome_A, Clinical_name_A,
                                    ifelse(Samples_names %in% Clinical_outcome_B, Clinical_name_B, 
                                           ifelse(Samples_names %in% Clinical_outcome_C, Clinical_name_C, ""))))
  Sample_timing <- factor(ifelse(Samples_names %in% Samplesheet$Baseline_RNAseq_data, "B", "REL"))
  if(typeof(Mutations_to_ignore) != "double"){
    Final_factor <- paste(Mutations_factor, Clinical_outcome, Sample_timing, sep = ".") %>% as.factor()
  }else{
    Final_factor <- paste(Clinical_outcome, Sample_timing, sep = ".") %>% as.factor()
  }
  
  return(Final_factor)
}

Factor_R_OR_NR_B <- Make_factor(Clinical_patient_data, 
                                colnames(RNAseq),
                                0, 
                                c("CR", "CRi"), 
                                "R", 
                                c("MLFS", "HI", "CRp", "PR"), 
                                "OR",
                                c("SD", "PD"), 
                                "NR")


Differential_analysis <- function(Focused_variable, DATA){
  design.pairs <- function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n - 1))
      for (j in (i + 1):n) {
        k <- k + 1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i], "-", levels[j],sep = "")
      }
    design
  }
  design <- model.matrix(~0 + Focused_variable)
  contr.matrix <- design.pairs(levels(factor(Focused_variable)))
  colnames(design) <- rownames(contr.matrix)   
  Fit <- lmFit(DATA, design) %>%
    contrasts.fit(., contr.matrix) %>%
    eBayes(., trend = TRUE)
  
  FitList <- list()
  for (i in 1:ncol(contr.matrix)) {
    FitList[[i]] <- topTable(Fit, coef = i, adjust.method = "BH", number = nrow(DATA)) %>%
      mutate(ID = rownames(.))
    
    message(paste0(i, " done"))
    
  }
  names(FitList) <- colnames(contr.matrix)
  return(FitList)
}

RNAseq_diff_gene_expression_analysis[["R_OR_NR_B"]] <- Differential_analysis(Factor_R_OR_NR_B, RNAseq)

####Miguel Functions for Viper

# Function to convert a dorothea table to viper regulons

# Function extracted from dorothea code
# https://github.com/saezlab/dorothea/blob/master/R/helpers.R#L17
dorothea2viper_regulons <- function(df) {
  regulon_list <- split(df, df$tf)
  viper_regulons <- lapply(regulon_list, function(regulon) {
    tfmode <- stats::setNames(regulon$mor, regulon$target)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })
  
  return(viper_regulons)
}


# Function to convert viper regulons to a dorothea table

# Function extracted from
# https://github.com/saezlab/ConservedFootprints/blob/master/src/dorothea_analysis.R#L126
viper_regulons2dorothea <- function(r) {
  res <- r %>%
    purrr::map_df(
      .f = function(i) {
        tf_target <- i$tfmode %>%
          tibble::enframe(name = "target", value = "mor") %>%
          mutate(likelihood = i$likelihood)
      },
      .id = "tf"
    )
  return(res)
}

# Function to convert dorothea database and gene expression to aracne regulons

# We need to create a network file where the first TF is the regulon with every target with a true confidence score
# in our case from Dorothea it always be 1
# Finally use it in aracne2regulon function from viper package
dorothea2aracne2viper_regulons <- function(dorothea, exprs_m) {
  dorothea_aggregation_tf <- dorothea %>%
    select(tf, target) %>%
    group_by(tf) %>%
    summarise(targets = str_c(target, collapse = ";"))
  tmp_file <- tempfile()
  for (i in 1:nrow(dorothea_aggregation_tf)) {
    tf_targets <- str_split(dorothea_aggregation_tf$targets[i], ";")[[1]]
    row <- c(dorothea_aggregation_tf$tf[i], unlist(mapply(c, tf_targets, rep(1, length(tf_targets)), SIMPLIFY = F)))
    cat(str_c(row, collapse = "\t"), "\n", file = tmp_file, append = T)
  }
  aracne_regulons <- aracne2regulon(tmp_file, exprs_m, format = "adj", verbose = F)
  file.remove(tmp_file)
  return(aracne_regulons)
}


# Main function to run msviper

run_msviper <- function(exprs_m, dorothea, use_aracne, ref, treat, ref_name, treat_name, minsize, ges.filter) {
  # First we need to generate the phenotype table (AnnotatedDataFrame)
  conditions <- rep("NA", ncol(exprs_m))
  conditions[ref] <- ref_name
  conditions[treat] <- treat_name
  names(conditions) <- colnames(exprs_m)
  conditions <- conditions[which(conditions != "NA")]
  
  phenotype <- data.frame(condition = factor(conditions))
  rownames(phenotype) <- names(conditions)
  
  phenoData <- new("AnnotatedDataFrame", data = phenotype)
  
  exprs_m <- exprs_m[,which(colnames(exprs_m) %in% rownames(phenotype))] %>% as.matrix()
  
  # Create Expression set from phenotyble table and expression matrix
  dset_viper <- ExpressionSet(assayData = exprs_m, phenoData = phenoData)
  dset_viper$sampleID <- factor(colnames(exprs_m))
  
  # Aracne can be used to estimate the mor instead using the -1, 1 from dorothea
  regulons <- NULL
  if (use_aracne) {
    regulons <- dorothea2aracne2viper_regulons(dorothea, dset_viper)
  } else {
    regulons <- dorothea2viper_regulons(dorothea)
  }
  
  # We need to create the statistics signature from the conditions
  signature <- rowTtest(dset_viper, "condition", treat_name, ref_name)
  statistics_signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
  # Generate the null model with bootstrapping (1000 iterations)
  nullmodel <- ttestNull(dset_viper, "condition", treat_name, ref_name, per = 1000, repos = T, verbose = F)
  # Run msviper using the statistics signature, the regulons converted from dorothea table, the null model the minSize of regulon and the ges.filter
  mrs <- msviper(ges = statistics_signature, regulon = regulons, nullmodel = nullmodel, minsize = minsize, ges.filter = ges.filter, verbose = F)
  # Convert the msviper regulons to dorothea
  dorothea_mrs_regulons <- viper_regulons2dorothea(mrs$regulon) %>%
    mutate(state = ifelse(mor > 0, "activation", "inhibition"))
  # Generate a table with the TFs, the regulon size, the NES score, the pval and the pval.fdr
  mrs_table <- tibble(TF = names(mrs$es$p.value), size = mrs$es$size, nes = mrs$es$nes, pval = mrs$es$p.value, pval.fdr = p.adjust(mrs$es$p.value, method = "fdr")) %>% arrange(pval)
  
  list(mrs_table = mrs_table, mrs = mrs, regulons = dorothea_mrs_regulons)
}


# Extra function to generate the cytoscape networok from msviper result


mrs2cytoscape <- function(mrs,full.path) {
  
  all_nodes <- unique(c(mrs$regulons$tf, mrs$regulons$target))
  tnodes <- tibble(TF = all_nodes)
  all_nodes_metadata <- right_join(mrs$mrs_table, tnodes, by = "TF")
  regulons_network <- graph.data.frame(mrs$regulons, directed = T, vertices = all_nodes_metadata)
  deleteAllNetworks()
  createNetworkFromIgraph(regulons_network, "regulons_network")
  # setVisualStyle(cytoscape_id_network, 'default')
  # setVisualStyle("default")
  
  my_style <- "my_style"
  
  
  
  
  
  createVisualStyle(my_style, list())
  setNodeColorDefault("#D3D3D3", style.name = my_style)
  blue_white_red <- c("#0000FF", "#FFFFFF", "#FF0000")
  setNodeColorMapping("nes", c(min(V(regulons_network)$nes, na.rm = T), mean(V(regulons_network)$nes, na.rm = T), max(V(regulons_network)$nes, na.rm = T)), blue_white_red, style.name = my_style)
  
  setEdgeTargetArrowShapeMapping("state", c("activation", "inhibition"), c("DELTA", "T"), style.name = my_style)
  
  setEdgeColorMapping("mor", c(min(E(regulons_network)$mor, na.rm = T), mean(E(regulons_network)$mor, na.rm = T), max(E(regulons_network)$mor, na.rm = T)), blue_white_red, style.name = my_style)
  setNodeLabelMapping('id'
  )
  setVisualStyle("my_style")
  
  createColumnFilter(filter.name='null', column='pval', 0.05, 'GREATER_THAN', network = regulons_network)
  applyFilter('null', hide=T, network = regulons_network)
  
  exportImage(full.path, 'SVG', zoom=200) 
  
}

data(dorothea_hs, package = "dorothea") 
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

ref_R_B <- Factor_R_OR_NR_B == "R.B"
ref_NR_B <- Factor_R_OR_NR_B == "NR.B"

R_NR_msviper <- run_msviper(RNAseq, regulons, use_aracne = T, ref_NR_B, ref_R_B,  "NR", "R", minsize = 4, ges.filter=T)

Focus_on_one_gene_not_TF <- function(RNAseq, Gene, Comparison_A, Comparison_A_name, Comparison_B, Comparison_B_name, phenotype){
  df <- RNAseq[rownames(RNAseq) == Gene, phenotype %in% c(Comparison_A, Comparison_B)]
  phenotype <- phenotype[phenotype %in% c(Comparison_A, Comparison_B)]
  pheno <- ifelse(phenotype == Comparison_A, Comparison_A_name, Comparison_B_name)
  df <- t(df) %>% as.data.frame()
  df$Phenotype <- pheno
  
  df[,1] <- as.numeric(df[,1])
  colnames(df)[1] <- "Gene"
  df
}

Make_gene_expr_boxplots <- function(RNAseq, Gene_to_focus, Comparison_A, Comparison_A_name, Comparison_B, Comparison_B_name, Phenotype){
  Phenotype_of_interest <- Phenotype
  
  Data_on_the_gene <- Focus_on_one_gene_not_TF(RNAseq, Gene_to_focus, Comparison_A, Comparison_A_name, Comparison_B, Comparison_B_name, Phenotype)
  ggplot(Data_on_the_gene, aes(x=Phenotype, y = Gene, fill=Phenotype))+
    geom_boxplot() +
    geom_jitter(Data_on_the_gene, inherit.aes = FALSE, mapping = aes(y = Gene, x = Phenotype), width = 0.25, alpha = 0.5, colour = "darkred")+
    ggtitle(paste("Expression", Gene_to_focus, Comparison_A_name, "vs", Comparison_B_name))+
    ylab("Gene expression")
}

# Network signatures

Prepare_features <- function(feature_data_frame, column_of_interest, type_of_data){
  if(type_of_data == "DEG"){
    res <- feature_data_frame[,column_of_interest]
    colnames(res)[1:3] <- c("Gene", "logFC", "P.Value")
  }else{
    res <- feature_data_frame[,column_of_interest]
    colnames(res)[1:3] <- c("Gene", "nes", "pval")
  }
  return(res)
}

Find_most_importants_genes <- function(network){
  res <- list()
  
  ranked_eigen_gene <- network$features[order(-network$features$Eigen_centrality),] %>% head(15) %>% .$Gene
  ranked_page_rank_gene <- network$features[order(-network$features$Page_rank),] %>% head(15) %>% .$Gene
  
  res$ranked_eigen_gene <- ranked_eigen_gene
  res$ranked_page_rank_gene <- ranked_page_rank_gene
  
  V_of_interest <- V(network$network) %>% .[which(names(.) %in% intersect(ranked_eigen_gene, ranked_page_rank_gene))]
  E_of_interest <- E(network$network)[from(V_of_interest) | to(V_of_interest)]
  
  filtered_graph <- subgraph.edges(network$network, E_of_interest)
  
  res$network <- filtered_graph
  
  return(res)
}

Prepare_Cytoscape_network <- function(Big_Network = igraph_PPI_TF_target_Network, DEG_analysis, TF_analysis, logFC_treshold = 0.75){
  DEG_of_interest <- DEG_analysis %>% dplyr::filter(abs(logFC) > logFC_treshold & P.Value < 0.1) %>% .$Gene
  TF_of_interest <- TF_analysis %>% dplyr::filter(pval < 0.1) %>% .$Gene
  
  V_of_interest <- V(Big_Network) %>% .[which(names(.) %in% unique(c(DEG_of_interest, TF_of_interest)))] 
  
  filtered_graph <- induced_subgraph(Big_Network, V_of_interest)
  
  
  eigen_centrality_result <- eigen_centrality(filtered_graph, directed = F)$vector
  
  page_rank_result <- igraph::page.rank(filtered_graph, directed = F)$vector
  
  features <- merge(DEG_analysis, TF_analysis, by = "Gene", all = T)
  features <- merge(features, eigen_centrality_result, by.x = "Gene", by.y = 0, all = T)
  colnames(features)[ncol(features)] <- "Eigen_centrality"
  features <- merge(features, page_rank_result, by.x = "Gene", by.y = 0, all = T)
  colnames(features)[ncol(features)] <- "Page_rank"
  
  set(features,which(is.na(features[["nes"]])),"nes",0)
  set(features,which(is.na(features[["pval"]])),"pval",1)
  set(features,which(is.na(features[["Eigen_centrality"]])),"Eigen_centrality",0)
  set(features,which(is.na(features[["Page_rank"]])),"Page_rank",0)
  
  features$TF <- ifelse(features$nes == 0, F, T)
  
  clustering_eigen <- cluster_leading_eigen(filtered_graph) %>% membership() %>% print() %>% data.frame()
  
  features <- merge(features, clustering_eigen, by.x = "Gene", by.y = 0, all = T)
  set(features,which(is.na(features[["."]])),".",999)
  colnames(features)[ncol(features)] <- "Cluster"
  
  res <- list("features" = features,
              "network" = filtered_graph
  )
  res$most_important_network <- Find_most_importants_genes(res)
  return(res)
}

All_workflow <- function(feature_DEG_df, column_DEG, feature_tf_df, column_TF, NET = igraph_PPI_TF_target_Network, logFC_treshold = 0.75){
  DEG <- Prepare_features(feature_DEG_df, column_DEG, "DEG")
  TF <- Prepare_features(feature_tf_df, column_TF, "TF")
  res <- Prepare_Cytoscape_network(NET, DEG, TF, logFC_treshold)
  return(res)
}

Combined_network <- read.csv("~/GitHub/Koichi_gene_expression_git/Koichi_gene_expression_analyses/Results/Tables/Combined_Networks.tsv", sep = "\t")
PPI_TF_target_Network <- graph_from_data_frame(Combined_network, directed = T)

NR_R_network <- All_workflow(RNAseq_diff_gene_expression_analysis$R_OR_NR_B$`NR.B-R.B`, c(7, 1, 4), R_NR_msviper$mrs_table, c(1,3,4), PPI_TF_target_Network)

write.csv(NR_R_network$features, "~/tmp/NR_R_network_features.csv", quote = F)
NR_R_network$network %>% igraph::as_data_frame() %>% write.csv("~/tmp/NR_R_network_network.csv", quote = F)


Do_cool_scatterplot <- function(Feature, title){
  Feature <- dplyr::filter(Feature, Eigen_centrality > 0.0005 & Page_rank != 0)
  DEG <- ifelse(Feature$logFC < 0, "DOWN", "UP")
  ggplot(Feature, aes(x = log(Page_rank), y = log(Eigen_centrality), label = Gene, colour = DEG))+
    geom_text(check_overlap = F, size = 2, nudge_x = 0.05, hjust = 0, outlier.size = 0)+
    geom_point(size = 0.5)+
    labs(title = paste0("Network-based node prioritization ", title))+
    xlab("Page Rank (log)")+
    ylab("Eigen Centrality (log)")+
    scale_colour_manual(values=c("#0000FF", "#FF0000"))
}

Variability_in_Patients <- data.frame(Variability = sapply(RNAseq, function(x){var(x)}),
                                      Pheno = Factor_R_OR_NR_B)



ggplot(Variability_in_Patients, aes(x= Pheno, y = Variability))+
  geom_boxplot()

gc()
