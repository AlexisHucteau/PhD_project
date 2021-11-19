library(hugene20sttranscriptcluster.db)
library(oligo)
library(dplyr)
library(FactoMineR)
library(ggsignif)
library(viper)
library(dorothea)
library(data.table)
library(stringr)
library(igraph)
library(ggplot2)

# Preparation of transcriptomes
Annotations <- data.frame(ACCNUM = sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse = ", "), 
                          SYMBOL = sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse = ", "), 
                          DESC = sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse = ", "))
SDRF <- read.csv("~/GitHub/Epigenomic_integration/Epigenomic_integration/DATA_RNAseq_Lucille/Samplesheet.sdrf.csv")
celFiles <- list.celfiles("~/GitHub/Epigenomic_integration/Epigenomic_integration/DATA_RNAseq_Lucille", full.names = T)
Transcriptomes <- read.celfiles(celFiles) %>% rma() %>% oligo::exprs() %>% as.matrix()
Transcriptomes <- merge(Annotations, Transcriptomes, by.x = 0, by.y = 0, all.y = T)
Transcriptomes <- tibble::column_to_rownames(Transcriptomes, var = "Row.names")
Annot <- Transcriptomes[,c(1:3)]
Transcriptomes <- Transcriptomes[,c(4:28)]
Phenotype <- paste(SDRF$Characteristics.cell.line., SDRF$Characteristics.genotype., SDRF$Characteristics.treatment., sep = ".")

PCA_analysis <- PCA(as.matrix(t(Transcriptomes)), graph = F)
# fviz_eig(PCA_analysis, addlabels = T)

DEGs <- function(data, Phenotype, Annotations){
  require(limma)
  require(dplyr)
  #require(dendextend)
  
  
  message("[===========================]")
  message("[<<<<<<< DGEs START >>>>>>>>>]")
  message("[<<<< Pairwise analysis >>>>>]")
  message("-----------------------------")
  
  # This function creates the pairs for the pairwise matrices
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
  
  # This function creates the pairs for the pairwise matrices
  
  design <- model.matrix(~0 + Phenotype)
  
  #Removing heteroscedascity from data
  
  contr.matrix <- design.pairs(levels(factor(Phenotype)))
  colnames(design) <- rownames(contr.matrix)   
  
  # Fitting linear models for comparisons of interest
  Fit <- lmFit(data, design) %>%
    contrasts.fit(., contr.matrix) %>%
    eBayes(., trend = TRUE)
  
  FitList <- list()
  for (i in 1:ncol(contr.matrix)) {
    FitList[[i]] <- topTable(Fit, coef = i, adjust.method = "BH", number = nrow(data)) %>%
      mutate(ID = rownames(.))
    FitList[[i]] <- merge(FitList[[i]], Annotations, by.x = 0, by.y = 0, all.x = T)
    
    message(paste0(i, " done"))
    
  }
  names(FitList) <- colnames(contr.matrix)
  return(FitList)
}

Differential_Cell_lines_analysis <- DEGs(Transcriptomes, Phenotype, Annot)


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

Transcriptomes_SYMBOL <- merge(Transcriptomes, Annot, by.x = 0, by.y = 0, all = T) %>%
  dplyr::filter(., SYMBOL != "NA") %>%
  split(., .$SYMBOL) %>% 
  lapply(., function(x){
    l <- length(x[1,])-3
    cnames <- colnames(x)[c(2:l)]
    df <- x[,c(2:l)] %>%
      as.matrix(.) %>%
      colMeans(.) %>% 
      data.frame(.) %>%
      t(.) %>%
      data.frame(.)
    colnames(df) <- cnames
    df
  }) %>%
  rbindlist(.) %>% 
  data.frame(.)

rownames(Transcriptomes_SYMBOL) <- unique(Annot$SYMBOL)[-1]
colnames(Transcriptomes_SYMBOL) <- str_remove_all(colnames(Transcriptomes_SYMBOL), "X")

HL60_Mut_IDHi_vs_no_treat <- run_msviper(Transcriptomes_SYMBOL, regulons, use_aracne = T, Phenotype == "HL60.Mut.AGI5198", Phenotype == "HL60.Mut.DMF",  "HL60_IDHi", "HL60_DMF", minsize = 4, ges.filter=T)
HL60_Mut_vs_HL60_WT <- run_msviper(Transcriptomes_SYMBOL, regulons, use_aracne = T, Phenotype == "HL60.Mut.None", Phenotype == "HL60.WT.None",  "HL60_Mut", "HL60_WT", minsize = 4, ges.filter=T)
MOLM14_Mut_IDHi_vs_no_treat <- run_msviper(Transcriptomes_SYMBOL, regulons, use_aracne = T, Phenotype == "MOLM14.Mut.AGI5198", Phenotype == "MOLM14.Mut.DMF",  "MOLM14_IDHi", "MOLM14_DMF", minsize = 4, ges.filter=T)

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

Make_gene_expr_boxplots <- function(RNAseq, Gene_to_focus, Gene_name, Comparison_A, Comparison_A_name, Comparison_B, Comparison_B_name, Phenotype){
  Phenotype_of_interest <- Phenotype
  
  Data_on_the_gene <- Focus_on_one_gene_not_TF(RNAseq, Gene_to_focus, Comparison_A, Comparison_A_name, Comparison_B, Comparison_B_name, Phenotype)
  ggplot(Data_on_the_gene, aes(x=Phenotype, y = Gene, fill=Phenotype))+
    geom_boxplot() +
    geom_jitter(Data_on_the_gene, inherit.aes = FALSE, mapping = aes(y = Gene, x = Phenotype), width = 0.25, alpha = 0.5, colour = "darkred")+
    ggtitle(paste("Expression", Gene_name, Comparison_A_name, "vs", Comparison_B_name))+
    ylab("Gene expression")
}

ID_RELA <- dplyr::filter(Annot, SYMBOL == "RELA") %>% rownames(.)
ID_MYC <- dplyr::filter(Annot, SYMBOL == "MYC") %>% rownames(.)

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

HL60_Mut_IDHi_vs_no_treat_network <- All_workflow(Differential_Cell_lines_analysis$`HL60.Mut.AGI5198-HL60.Mut.DMF`, c(10, 2, 5), HL60_Mut_IDHi_vs_no_treat$mrs_table, c(1,2,3), PPI_TF_target_Network)
HL60_Mut_vs_HL60_WT_network  <- All_workflow(Differential_Cell_lines_analysis$`HL60.Mut.None-HL60.WT.None`, c(10, 2, 5), HL60_Mut_vs_HL60_WT$mrs_table, c(1,2,3), PPI_TF_target_Network)
MOLM14_Mut_IDHi_vs_no_treat_network <- All_workflow(Differential_Cell_lines_analysis$`MOLM14.Mut.AGI5198-MOLM14.Mut.DMF`, c(10, 2, 5), MOLM14_Mut_IDHi_vs_no_treat$mrs_table, c(1,2,3), PPI_TF_target_Network)

write.csv(HL60_Mut_IDHi_vs_no_treat_network$features, "~/tmp/HL60_IDHi.tsv", quote = F)
HL60_Mut_IDHi_vs_no_treat_network$network %>% igraph::as_data_frame() %>% write.csv("~/tmp/HL60_IDHi_net.csv", quote = F)

write.csv(HL60_Mut_vs_HL60_WT_network$features, "~/tmp/HL60_Mut_vs_HL60_WT_features.csv", quote = F)
HL60_Mut_vs_HL60_WT_network$network %>% igraph::as_data_frame() %>% write.csv("~/tmp/HL60_Mut_vs_HL60_WT_network.csv", quote = F)


write.csv(MOLM14_Mut_IDHi_vs_no_treat_network$features, "~/tmp/MOLM14_Mut_IDHi_vs_no_treat_features.csv", quote = F)
MOLM14_Mut_IDHi_vs_no_treat_network$network %>% igraph::as_data_frame() %>% write.csv("~/tmp/MOLM14_Mut_IDHi_vs_no_treat_network.csv", quote = F)

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

HL60_m_vs_M14_m <- run_msviper(Transcriptomes_SYMBOL, 
                                   regulons, use_aracne = T, 
                                   Phenotype == "HL60.Mut.DMF", Phenotype == "MOLM14.Mut.DMF",  
                                   "HL60_Mut", "M14_Mut", 
                                   minsize = 4, ges.filter=T)

HL60_m_agi_vs_M14_m_agi <- run_msviper(Transcriptomes_SYMBOL, 
                                           regulons, use_aracne = T, 
                                           Phenotype == "HL60.Mut.AGI5198", Phenotype == "MOLM14.Mut.AGI5198",  
                                           "HL60_IDHi", "M14_IDHi", 
                                           minsize = 4, ges.filter=T)

HL60_m_vs_M14_m_network <- All_workflow(Differential_Cell_lines_analysis$`HL60.Mut.DMF-MOLM14.Mut.DMF`, c(10, 2, 5),
                                        HL60_m_vs_M14_m$mrs_table, c(1,2,3), 
                                                  PPI_TF_target_Network)

HL60_m_agi_vs_M14_m_agi_network  <- All_workflow(Differential_Cell_lines_analysis$`HL60.Mut.AGI5198-MOLM14.Mut.AGI5198`, c(10, 2, 5), 
                                                 HL60_m_agi_vs_M14_m_agi$mrs_table, c(1,2,3), 
                                             PPI_TF_target_Network)

write.csv(HL60_m_vs_M14_m_network$features, "~/tmp/HL60_M14_IDHm_features.csv", quote = F, row.names = F)
HL60_m_vs_M14_m_network$network %>% igraph::as_data_frame() %>% write.csv("~/tmp/HL60_M14_IDHm_net.csv", quote = F, row.names = F)

write.csv(HL60_m_agi_vs_M14_m_agi_network$features, "~/tmp/HL60_M14_IDHi_features.csv", quote = F, row.names = F)
HL60_m_agi_vs_M14_m_agi_network$network %>% igraph::as_data_frame() %>% write.csv("~/tmp/HL60_M14_IDHi_network.csv", quote = F)

Variability_in_cell_lines <- data.frame(Variability = sapply(Transcriptomes, function(x){var(x)}),
                                      Pheno = Phenotype)
Variability_in_cell_lines <- dplyr::filter(Variability_in_cell_lines, Pheno %in% c("HL60.Mut.AGI5198", "HL60.Mut.DMF", "MOLM14.Mut.DMF", "MOLM14.Mut.AGI5198"))

gc()
