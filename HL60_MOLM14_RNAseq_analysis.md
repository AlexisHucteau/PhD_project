# Analysis of HL60 MOLM14 IDHm cell lines +/- IDHi

# 1. Preparation of transcriptomes

# 2. PCA analysis

    PCA_analysis <- PCA(as.matrix(t(Transcriptomes)), graph = T)

    ## Warning: ggrepel: 6 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/unnamed-chunk-1-1.png)![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/unnamed-chunk-1-2.png)

    fviz_eig(PCA_analysis, addlabels = T)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/unnamed-chunk-1-3.png)

    plot(HL60_Mut_IDHi_vs_no_treat$mrs)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/Do%20TF%20plot-1.png)

    plot(HL60_Mut_vs_HL60_WT$mrs)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/Do%20TF%20plot-2.png)

    plot(MOLM14_Mut_IDHi_vs_no_treat$mrs)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/Do%20TF%20plot-3.png)

# 3. RELA and MYC expression in cell lines

    ID_RELA <- dplyr::filter(Annot, SYMBOL == "RELA") %>% rownames(.)
    ID_MYC <- dplyr::filter(Annot, SYMBOL == "MYC") %>% rownames(.)


    Make_gene_expr_boxplots(Transcriptomes, ID_RELA, "RELA", "HL60.Mut.AGI5198", "HL60_IDHi", "HL60.Mut.DMF", "HL60_DMF", Phenotype)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/Making%20box%20plot%20of%20Expression-1.png)

    Make_gene_expr_boxplots(Transcriptomes, ID_MYC, "MYC", "HL60.Mut.AGI5198", "HL60_IDHi", "HL60.Mut.DMF", "HL60_DMF", Phenotype)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/Making%20box%20plot%20of%20Expression-2.png)

    Make_gene_expr_boxplots(Transcriptomes, ID_RELA, "RELA", "HL60.Mut.None", "HL60_Mut", "HL60.WT.None", "HL60_WT", Phenotype)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/Making%20box%20plot%20of%20Expression-3.png)

    Make_gene_expr_boxplots(Transcriptomes, ID_MYC, "MYC", "HL60.Mut.None", "HL60_Mut", "HL60.WT.None", "HL60_WT", Phenotype)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/Making%20box%20plot%20of%20Expression-4.png)

    Make_gene_expr_boxplots(Transcriptomes, ID_RELA, "RELA", "MOLM14.Mut.AGI5198", "MOLM14_IDHi", "MOLM14.Mut.DMF", "MOLM14_DMF", Phenotype)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/Making%20box%20plot%20of%20Expression-5.png)

    Make_gene_expr_boxplots(Transcriptomes, ID_MYC, "MYC", "MOLM14.Mut.AGI5198", "MOLM14_IDHi", "MOLM14.Mut.DMF", "MOLM14_DMF", Phenotype)

![](HL60_MOLM14_RNAseq_analysis_files/figure-markdown_strict/Making%20box%20plot%20of%20Expression-6.png)

# 4. Network signatures

### HL60 IDHm + IDHi vs HL60 IDHm + DMF
![HL60 IDHi](Pictures/HL60_IDHi.png)

### HL60 IDHm vs HL60 IDHwt
![HL60 m/wt](Pictures/HL60_Mut.png)

### MOLM14 IDHm + IDHi vs MOLM14 IDHm + DMF
![MOLM14 IDHi](Pictures/MOLM14_IDHi.png)
