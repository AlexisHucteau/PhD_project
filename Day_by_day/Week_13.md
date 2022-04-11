# 28.03.2022 - 01.04.2022 Week 13

## 28.03.2022

Looking at aracn effect on msviper TF activity prediction

Testing DeMAND-algorithm for Interrogating Drug Mechanism of
Action using Network Dysregulation analysis

Using DeMAND, Mechanism of Action dyregulated NR vs R:

| moaGene |     Pvalue |        FDR |
|----|----|----|
|**NOS2**   |   1.517842e-11|	4.551391e-08|
|**CYP17A1**	|  2.937329e-11	|4.551391e-08|
|**BIRC3**	  |  3.851648e-09	|3.978752e-06|
|**IL2**	   |   8.886857e-08	|6.885092e-05|
|**SELE**  |  2.475087e-07	|1.534059e-04|
|**BCL2L1**	  |9.864902e-07	|4.485631e-04|

Kullback-Leibler (KL) divergence:

gene1|gene2|KLD|KLD.p
----|----|----|----
**DLGAP1**|**SOX2**|9.08160044288928|0.001842360444131
**MYC**|**WDR36**|8.27935400426403|0.002225292942782
**MYC**|**SFMBT1**|8.27767843611966|0.002226212776101
**CYP17A1**|**NFIC**|8.0767423217575|0.00234076542414
**PRF1**|**STAT3**|8.04440062134892|0.002360021079851
**JUN**|**NOS2**|8.03554034015082|0.002365337512273
**NOS2**|**RELA**|7.59128947444795|0.002656597252669
**FOS**|**NOS2**|7.4782755495627|0.0027392172276

About NOS2: "Up-Regulation of iNOS in AML Blasts Creates an Immunosuppressive Microenvironment, Inhibits T-Cell Proliferation and Transforms T-Cells Towards a Tumor-Tolerating Phenotype"

In myeloid cells, the NF-κB complexes that bind to the nos2 promoter are p65/p65 and p50/p50 homodimers

## 29.03.2022

Downloading TCGA methylomes
Running the RNAseq pipeline on PTBP1 data

## 30.03.2022

Trying to make a pchic network of genes of interest.
Idea: Instead of 1 pchic w/ genes/tfs of interest + DMP/DMR
-> 1 pchic w/ genes/tf of interest
-> 1 pchic w/ DMR/DMP
