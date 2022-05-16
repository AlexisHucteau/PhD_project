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
library(RCy3)

# Preparation of transcriptomes
Annotations <- data.frame(ACCNUM = sapply(contents(hugene20sttranscriptclusterACCNUM), paste, collapse = ", "), 
                          SYMBOL = sapply(contents(hugene20sttranscriptclusterSYMBOL), paste, collapse = ", "), 
                          DESC = sapply(contents(hugene20sttranscriptclusterGENENAME), paste, collapse = ", "))
SDRF <- read.csv("~/GitHub/Epigenomic_integration/DATA_RNAseq_Lucille/Samplesheet.sdrf.csv")
celFiles <- list.celfiles("~/GitHub/Epigenomic_integration/DATA_RNAseq_Lucille", full.names = T)
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
    dplyr::select(tf, target) %>%
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

# Transcriptomes_SYMBOL <- merge(Transcriptomes, Annot, by.x = 0, by.y = 0, all = T) %>%
#   dplyr::filter(., SYMBOL != "NA") %>%
#   split(., .$SYMBOL) %>% 
#   lapply(., function(x){
#     l <- length(x[1,])-3
#     cnames <- colnames(x)[c(2:l)]
#     df <- x[,c(2:l)] %>%
#       as.matrix(.) %>%
#       colMeans(.) %>% 
#       data.frame(.) %>%
#       t(.) %>%
#       data.frame(.)
#     colnames(df) <- cnames
#     df
#   }) %>%
#   rbindlist(.) %>% 
#   data.frame(.)
# 
# rownames(Transcriptomes_SYMBOL) <- unique(Annot$SYMBOL)[-1]
# colnames(Transcriptomes_SYMBOL) <- str_remove_all(colnames(Transcriptomes_SYMBOL), "X")

Transcriptomes_SYMBOL <- read.csv("GitHub/Koichi_gene_expression_analyses_git/Koichi_gene_expression_analyses/DATA/Cell_line_Transcriptomes_symbol.csv", check.names = F, row.names = 1)

HL60_Mut_IDHi_vs_no_treat <- run_msviper(Transcriptomes_SYMBOL, regulons, use_aracne = T, Phenotype == "HL60.Mut.DMF", Phenotype == "HL60.Mut.AGI5198", "HL60_DMF", "HL60_IDHi", minsize = 4, ges.filter=T)
HL60_Mut_vs_HL60_WT <- run_msviper(Transcriptomes_SYMBOL, regulons, use_aracne = T, Phenotype == "HL60.WT.None", Phenotype == "HL60.Mut.None", "HL60_WT", "HL60_Mut", minsize = 4, ges.filter=T)
MOLM14_Mut_IDHi_vs_no_treat <- run_msviper(Transcriptomes_SYMBOL, regulons, use_aracne = T, Phenotype == "MOLM14.Mut.DMF", Phenotype == "MOLM14.Mut.AGI5198", "MOLM14_DMF", "MOLM14_IDHi", minsize = 4, ges.filter=T)
HL60_Mut_IDHi_vs_Molm14_mut_IDHi <- run_msviper(Transcriptomes_SYMBOL, regulons, use_aracne = T, Phenotype == "MOLM14.Mut.AGI5198", Phenotype == "HL60.Mut.AGI5198", "MOLM14_IDHi", "HL60_IDHi", minsize = 4, ges.filter=T)
HL60_Mut_vs_Molm14_mut <- run_msviper(Transcriptomes_SYMBOL, regulons, use_aracne = T, Phenotype == "MOLM14.Mut.DMF", Phenotype == "HL60.Mut.DMF", "MOLM14_DMF", "HL60_DMF", minsize = 4, ges.filter=T)

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
    res <- res[res$Gene != "NA",]
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
  
  V_of_interest <- V(Big_Network) %>% .[which(names(.) %in% c(DEG_of_interest, TF_of_interest))]
  
  filtered_graph <- induced_subgraph(Big_Network, V_of_interest)
  
  eigen_centrality_result <- eigen_centrality(filtered_graph, directed = F)$vector
  
  page_rank_result <- igraph::page.rank(filtered_graph, directed = F)$vector
  
  features <- merge(DEG_analysis, TF_analysis, by = "Gene", all = T)
  features <- dplyr::filter(features, (abs(logFC) > logFC_treshold & P.Value < 0.1) | (pval < 0.1))
  features <- merge(features, eigen_centrality_result, by.x = "Gene", by.y = 0, all = T)
  colnames(features)[ncol(features)] <- "Eigen_centrality"
  features <- merge(features, page_rank_result, by.x = "Gene", by.y = 0, all = T)
  colnames(features)[ncol(features)] <- "Page_rank"
  
  set(features,which(is.na(features[["nes"]])),"nes",0)
  set(features,which(is.na(features[["pval"]])),"pval",1)
  set(features,which(is.na(features[["Eigen_centrality"]])),"Eigen_centrality",0)
  set(features,which(is.na(features[["Page_rank"]])),"Page_rank",0)
  
  features$TF <- ifelse(features$nes == 0, F, T)
  
  res <- list("features" = features,
              "network" = filtered_graph
  )
  res$most_important_network <- Find_most_importants_genes(res)
  return(res)
}

FIsInGene_020720_with_annotations <- read.csv("Documents/List_genes_from_transcriptome/FIsInGene_020720_with_annotations.tsv", sep = "\t")

build_network_from_aracn <- function(PPI_net = FIsInGene_020720_with_annotations, aracn_regulon){
  aracn_regulon$merging <- paste(aracn_regulon$tf, aracn_regulon$target, sep = ",")
  PPI_net$merging <- paste(PPI_net$Gene1, PPI_net$Gene2, sep = ",")
  merged <- merge(aracn_regulon, PPI_net, by = "merging", all.x = T, all.y = T)
  merged$mor <- ifelse(is.na(merged$mor), 0, merged$mor)
  merged$tf <- ifelse(is.na(merged$tf), "", merged$tf)
  merged$target <- ifelse(is.na(merged$target), "", merged$target)
  merged$state <- ifelse(is.na(merged$state), "", merged$state)
  merged$likelihood <- ifelse(is.na(merged$likelihood), "", merged$likelihood)
  merged
}

All_workflow <- function(feature_DEG_df, column_DEG, feature_tf_df, column_TF, Aracn_net, PPI_net = FIsInGene_020720_with_annotations, logFC_treshold = 0.75){
  DEG <- Prepare_features(feature_DEG_df, column_DEG, "DEG")
  TF <- Prepare_features(feature_tf_df, column_TF, "TF")
  NET <- build_network_from_aracn(PPI_net, Aracn_net)[,c(7,8,2:6,9:11)]
  igraph_net <- graph_from_data_frame(NET, directed = T)
  res <- Prepare_Cytoscape_network(igraph_net, DEG, TF, logFC_treshold)
  return(res)
}

Stylish_the_network <- function(PPI = T, TF = T, network_features, network_edges_features, title){
  defaults <- list(NODE_SHAPE="diamond",
                   NODE_SIZE=30,
                   EDGE_TRANSPARENCY=120,
                   NODE_LABEL_POSITION="c,c,c,0.00,0.00")
  nodeLabels <- mapVisualProperty(visual.prop = 'node label', table.column = 'name', mapping.type = 'p')
  createVisualStyle(title, defaults, list(nodeLabels))
  setVisualStyle(title)
  setNodeShapeMapping(table.column = 'TF', 
                      table.column.values = c(T, F), 
                      shapes = c('ELLIPSE', "RECTANGLE"), 
                      style.name = title)
  
  setNodeColorMapping(table.column = 'logFC', 
                      table.column.values = c(min(network_features$logFC), 0.0, max(network_features$logFC)), 
                      colors = c ('#0000FF', '#FFFFFF', '#FF0000'), 
                      style.name = title)
  
  setNodeFillOpacityMapping(table.column = 'P.Value', 
                            table.column.values = c(0, 0.1, 1), 
                            opacities = c(255, 200, 100), 
                            style.name = title)
  
  setNodeSizeMapping (table.column = 'Eigen_centrality', 
                      table.column.values = c(0, 1), 
                      sizes = c(20, 200), 
                      style.name = title)
  
  setNodeFontSizeMapping(table.column = 'Eigen_centrality', 
                         table.column.values = c(0, 1), 
                         sizes = c(10, 75), 
                         style.name = title)
  
  if (TF){
    setEdgeLineWidthMapping(table.column = 'mor',
                            table.column.values = c(min(network_edges_features$mor), 0.0, max(network_edges_features$mor)),
                            widths = c(10, 2, 10), 
                            style.name = title)
    
    setEdgeColorMapping(table.column = 'mor',
                        table.column.values = c(min(network_edges_features$mor), 0.0, max(network_edges_features$mor)),
                        colors = c('#0000FF', "#555555", '#FF0000'), 
                        style.name = title)
    
    setEdgeOpacityMapping(table.column = 'mor',
                          table.column.values = c(min(network_edges_features$mor), 0.0, max(network_edges_features$mor)),
                          opacities = c(255, 100, 255), 
                          style.name = title)
  }
  if(PPI){
    setEdgeTargetArrowShapeMapping(table.column = 'Direction',
                                   table.column.values = names(table(network_edges_features$Direction)),
                                   shapes = c('NONE', 'ARROW', 'T', 'NONE', 'ARROW', 'T', 'NONE', 'ARROW'), 
                                   style.name = title)
  }else{
    createColumnFilter(filter.name = "Activation", column = "mor", criterion = 0, type = "edges", predicate = "GREATER_THAN")
    setEdgeTargetArrowShapeBypass(getSelectedEdges(), "ARROW")
    
    createColumnFilter(filter.name = "Inhibition", column = "mor", criterion = 0, type = "edges", predicate = "LESS_THAN")
    setEdgeTargetArrowShapeBypass(getSelectedEdges(), "T")
  }
  
  createDegreeFilter(filter.name = "single", criterion = c(0,0))
  getSelectedNodes()
  deleteSelectedNodes()
  
  createColumnFilter(filter.name = "NES_DOWN", column = "nes", criterion = 0, type = "nodes", predicate = "LESS_THAN")
  setNodeColorBypass(getSelectedNodes(), new.colors = "#C411FF")
  createColumnFilter(filter.name = "NES_UP", column = "nes", criterion = 0, type = "nodes", predicate = "GREATER_THAN")
  setNodeColorBypass(getSelectedNodes(), new.colors = "#5AFF00")  
  createColumnFilter(filter.name = "TF_activity_pval", column = "pval", criterion = c(0.1,0.99), type = "nodes", predicate = "BETWEEN")
  setNodeColorBypass(getSelectedNodes(), new.colors = "#BBBBBB")
}

df2cytos_cell_lines <- function(Aracn_network = T, Aracn_net, DEGs, TFs, PPI_network = FIsInGene_020720_with_annotations, PPI = T, Network_title, Networks_collection){
  DEG <- Prepare_features(DEGs, c(10, 2, 5), "DEG")
  TF <- Prepare_features(TFs, c(1, 3, 4), "TF")
  if (PPI & Aracn_network){
    message("Both")
    NET <- build_network_from_aracn(PPI_network, Aracn_net)[,c(7,8,2:6,9:11)]
  }else if(PPI){
    message("Only PPI")
    NET <- PPI_network
  }else{
    message("only TFs")
    NET <- Aracn_net
  }
  iNet <- graph_from_data_frame(NET, directed = T)
  res <- Prepare_Cytoscape_network(iNet, DEG, TF)
  
  NET <- igraph::as_data_frame(res$network)
  colnames(NET)[c(1,2)] <- c("source","target")
  colnames(res$features)[1] <- "id"
  createNetworkFromDataFrames(edges = NET, nodes = res$features, title = Network_title, collection = Networks_collection)
  Stylish_the_network(PPI, Aracn_network, res$features, NET, title = Network_title)
  return(res)
}


HL60_Mut_IDHi_vs_no_treat_only_TFs <- df2cytos_cell_lines(Aracn_network = T, 
                                               HL60_Mut_IDHi_vs_no_treat$regulons, 
                                               Differential_Cell_lines_analysis$`HL60.Mut.AGI5198-HL60.Mut.None`, 
                                               HL60_Mut_IDHi_vs_no_treat$mrs_table, 
                                               FIsInGene_020720_with_annotations, 
                                               F, 
                                               "HL60_Mut_IDHi_vs_no_treat_only_TFs", 
                                               "Cell line collection")

HL60_Mut_IDHi_vs_no_treat_Only_PPI <- df2cytos_cell_lines(Aracn_network = F, 
                                     HL60_Mut_IDHi_vs_no_treat$regulons, 
                                     Differential_Cell_lines_analysis$`HL60.Mut.AGI5198-HL60.Mut.None`, 
                                     HL60_Mut_IDHi_vs_no_treat$mrs_table, 
                                     FIsInGene_020720_with_annotations, 
                                     T, 
                                     "HL60_Mut_IDHi_vs_no_treat_Only_PPI", 
                                     "Cell line collection")

HL60_Mut_IDHi_vs_no_treat_Both_net <- df2cytos_cell_lines(Aracn_network = T, 
                                      HL60_Mut_IDHi_vs_no_treat$regulons, 
                                      Differential_Cell_lines_analysis$`HL60.Mut.AGI5198-HL60.Mut.None`, 
                                      HL60_Mut_IDHi_vs_no_treat$mrs_table, 
                                      FIsInGene_020720_with_annotations, 
                                      T, 
                                      "HL60_Mut_IDHi_vs_no_treat_both_net", 
                                      "Cell line collection")

HL60_Mut_vs_HL60_WT_only_TFs <- df2cytos_cell_lines(Aracn_network = T, 
                                         HL60_Mut_vs_HL60_WT$regulons, 
                                         Differential_Cell_lines_analysis$`HL60.Mut.DMF-HL60.WT.None`, 
                                         HL60_Mut_vs_HL60_WT$mrs_table, 
                                         FIsInGene_020720_with_annotations, 
                                         F, 
                                         "HL60_Mut_vs_HL60_WT_only_TFs", 
                                         "Cell line collection")

HL60_Mut_vs_HL60_WT_Only_PPI <-df2cytos_cell_lines(Aracn_network = F, 
                               HL60_Mut_vs_HL60_WT$regulons, 
                               Differential_Cell_lines_analysis$`HL60.Mut.DMF-HL60.WT.None`, 
                               HL60_Mut_vs_HL60_WT$mrs_table, 
                               FIsInGene_020720_with_annotations, 
                               T, 
                               "HL60_Mut_vs_HL60_WT", 
                               "Cell line collection")

HL60_Mut_vs_HL60_WT_Both_net <-df2cytos_cell_lines(Aracn_network = T, 
                               HL60_Mut_vs_HL60_WT$regulons, 
                               Differential_Cell_lines_analysis$`HL60.Mut.DMF-HL60.WT.None`, 
                               HL60_Mut_vs_HL60_WT$mrs_table, 
                               FIsInGene_020720_with_annotations, 
                               T, 
                               "HL60_Mut_vs_HL60_WT_Both_net", 
                               "Cell line collection")

MOLM14_Mut_IDHi_vs_no_treat_only_TFs <- df2cytos_cell_lines(Aracn_network = T, 
                                                 MOLM14_Mut_IDHi_vs_no_treat$regulons, 
                                                 Differential_Cell_lines_analysis$`MOLM14.Mut.AGI5198-MOLM14.Mut.DMF`, 
                                                 MOLM14_Mut_IDHi_vs_no_treat$mrs_table, 
                                                 FIsInGene_020720_with_annotations, 
                                                 F, 
                                                 "MOLM14_Mut_IDHi_vs_no_treat_only_TFs", 
                                                 "Cell line collection")

MOLM14_Mut_IDHi_vs_no_treat_Only_PPI <-df2cytos_cell_lines(Aracn_network = F, 
                                       MOLM14_Mut_IDHi_vs_no_treat$regulons, 
                                       Differential_Cell_lines_analysis$`MOLM14.Mut.AGI5198-MOLM14.Mut.DMF`, 
                                       MOLM14_Mut_IDHi_vs_no_treat$mrs_table, 
                                       FIsInGene_020720_with_annotations, 
                                       T, 
                                       "MOLM14_Mut_IDHi_vs_no_treat_Only_PPI", 
                                       "Cell line collection")

MOLM14_Mut_IDHi_vs_no_treat_Both_net <- df2cytos_cell_lines(Aracn_network = T, 
                                       MOLM14_Mut_IDHi_vs_no_treat$regulons, 
                                       Differential_Cell_lines_analysis$`MOLM14.Mut.AGI5198-MOLM14.Mut.DMF`, 
                                       MOLM14_Mut_IDHi_vs_no_treat$mrs_table, 
                                       FIsInGene_020720_with_annotations, 
                                       T, 
                                       "MOLM14_Mut_IDHi_vs_no_treat_Both_net", 
                                       "Cell line collection")

HL60_Mut_IDHi_vs_Molm14_mut_IDHi_only_TFs <- df2cytos_cell_lines(Aracn_network = T, 
                                                      HL60_Mut_IDHi_vs_Molm14_mut_IDHi$regulons, 
                                                      Differential_Cell_lines_analysis$`HL60.Mut.AGI5198-MOLM14.Mut.AGI5198`, 
                                                      HL60_Mut_IDHi_vs_Molm14_mut_IDHi$mrs_table, 
                                                      FIsInGene_020720_with_annotations, 
                                                      F, 
                                                      "HL60_Mut_IDHi_vs_Molm14_mut_IDHi_only_TFs", 
                                                      "Cell line collection")

HL60_Mut_IDHi_vs_Molm14_mut_IDHi_Only_PPI <-df2cytos_cell_lines(Aracn_network = F, 
                                            HL60_Mut_IDHi_vs_Molm14_mut_IDHi$regulons, 
                                            Differential_Cell_lines_analysis$`HL60.Mut.AGI5198-MOLM14.Mut.AGI5198`, 
                                            HL60_Mut_IDHi_vs_Molm14_mut_IDHi$mrs_table, 
                                            FIsInGene_020720_with_annotations,
                                            T, 
                                            "HL60_Mut_IDHi_vs_Molm14_mut_IDHi_Only_PPI", 
                                            "Cell line collection")

HL60_Mut_IDHi_vs_Molm14_mut_IDHi_Both_net <-df2cytos_cell_lines(Aracn_network = T, 
                                            HL60_Mut_IDHi_vs_Molm14_mut_IDHi$regulons, 
                                            Differential_Cell_lines_analysis$`HL60.Mut.AGI5198-MOLM14.Mut.AGI5198`, 
                                            HL60_Mut_IDHi_vs_Molm14_mut_IDHi$mrs_table, 
                                            FIsInGene_020720_with_annotations, 
                                            T, 
                                            "HL60_Mut_IDHi_vs_Molm14_mut_IDHi_Both_net", 
                                            "Cell line collection")

HL60_Mut_vs_Molm14_mut_only_TFs <- df2cytos_cell_lines(Aracn_network = T, 
                                            HL60_Mut_vs_Molm14_mut$regulons, 
                                            Differential_Cell_lines_analysis$`HL60.Mut.DMF-MOLM14.Mut.DMF`, 
                                            HL60_Mut_vs_Molm14_mut$mrs_table, 
                                            FIsInGene_020720_with_annotations, 
                                            F, 
                                            "HL60_Mut_vs_Molm14_mut_only_TFs", 
                                            "Cell line collection")

HL60_Mut_vs_Molm14_mut_Only_PPI <- df2cytos_cell_lines(Aracn_network = F, 
                                  HL60_Mut_vs_Molm14_mut$regulons, 
                                  Differential_Cell_lines_analysis$`HL60.Mut.DMF-MOLM14.Mut.DMF`, 
                                  HL60_Mut_vs_Molm14_mut$mrs_table, 
                                  FIsInGene_020720_with_annotations, 
                                  T, 
                                  "HL60_Mut_vs_Molm14_mut_Only_PPI", 
                                  "Cell line collection")

HL60_Mut_vs_Molm14_mut_Both_net <- df2cytos_cell_lines(Aracn_network = T, 
                                   HL60_Mut_vs_Molm14_mut$regulons, 
                                   Differential_Cell_lines_analysis$`HL60.Mut.DMF-MOLM14.Mut.DMF`, 
                                   HL60_Mut_vs_Molm14_mut$mrs_table, 
                                   FIsInGene_020720_with_annotations, 
                                   T, 
                                   "HL60_Mut_vs_Molm14_mut_Both_net", 
                                   "Cell line collection")


Do_cool_scatterplot <- function(Feature, title){
  Feature <- dplyr::filter(Feature, Eigen_centrality > 0.0005 & Page_rank != 0 & ((P.Value < 0.05 & abs(logFC) > 1.5) | pval < 0.05))
  DEG <- ifelse(Feature$logFC < 0, "DOWN", "UP")
  ggplot(Feature, aes(x = log(Page_rank), y = log(Eigen_centrality), label = Gene, colour = DEG))+
    geom_text(check_overlap = T, size = 2, nudge_x = 0.05, hjust = 0, outlier.size = 0)+
    geom_point(size = 0.5)+
    labs(title = paste0("Network-based node prioritization ", title))+
    xlab("Page Rank (log)")+
    ylab("Eigen Centrality (log)")+
    scale_colour_manual(values=c("#0000FF", "#FF0000"))
}


gc()
