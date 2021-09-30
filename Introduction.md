# Introduction

The main topic of the thesis is in silico, epigenomic and functional investigations of resistance to IDH inhibitors in acute myeloid leukemia.

Acute myeloid leukemia (AML) is a deadly disease associated with poor outcomes. The median age of patients harbouring AML is around 70 years old and the prognoses decline with the age. For decades, intensive chemotherapies were the only treatments available. The combination of Cytarabine (AraC) and daunorubicin shows good results but not all the patients are eligible for intensive therapy. Hypomethylating agents are the current alternative therapy but only 16% of patients achieve complete remission with a median overall survival about 21 months.

Novel therapies have been accepted like the combinaison Venetoclax + Azacitidine to avoid intensive therapy but there is still a crucial need of novel therapies to overcome the relapses and refractory diseases. Although 60–80% of AML adult patients achieve complete remission after the first induction chemotherapy, roughly 20% will show primary refractory disease and more than 50% will relapse.

Specific mutations have been found to induce the relapses like FLT3-ITD, IDH1/2 or the presence of CD33 and inhibitors are already studied but there are still resistances to that therapies. IDH1 and IDH2 are metabolic genes that convert isocitrate into alpha ketoglutarate (aKG), a metabolite that is itself converted into 2 hydroxyglutarate (2HG) by the mutated IDH protein. The abnormal production of the oncometabolite 2-hydroxyglutarate (2-HG) induces epigenetic and transcriptional reprogramming, differentiation bias, and susceptibility to mitochondrial inhibitors in cancer cells.

In a recent study, it have been shown that patients with AML harboring an IDH mutation displayed an enhanced mitochondrial oxidative metabolism through the methylation-driven CEBPa induction but inhibition of the mutation, despite the reduction of CEBPa enhancement, failed to reverse the mitochondrial oxidative metabolism through an AKT/PPARGg pathway (ref: [Mitochondrial metabolism supports resistance to IDH mutant inhibitors in acute myeloid leukemia](https://pubmed.ncbi.nlm.nih.gov/33760042/)).

Another study from our group showed that the presence of the mutation induces Vitamin D receptor related programs priming the AML cells to differentiate with pharmacological doses of ATRA or/and VD comforting the relevance of metabolic pathway in the mutation(ref: [Activation of Vitamin D Receptor Pathway Enhances Differentiating Capacity in Acute Myeloid Leukemia with Isocitrate Dehydrogenase Mutations](https://www.preprints.org/manuscript/202108.0529/v1)).

Another study focusing on the inhibition of IDH suggests that stemness is associated with primary resistance and selection of mutations in RUNX1/CEBPA or RAS/RTX pathway genes are the driver of acquired resistance (ref: [Leukemia stemness and co-occurring mutations drive resistance to IDH inhibitors in acute myeloid leukemia](https://www.nature.com/articles/s41467-021-22874-x)). Despite those advances, the mechanism of the resistances remain unclear.

To investigate the resistance, we analyzed transcriptomes of different cohorts of patients. TCGA and Verhaak datasets have been used to compare the patients samples harbouring IDH mutation to patients samples without the mutation. Verhaak data have been used to look at primary resistance through the analysis of Overall Survival and the patient cohort from Koichi’s data have been used to investigate the acquired resistance by looking at relapsed or refractory AML behaviour after IDH inhibitor therapy.

By combining transcription factor activity, differential gene expression, functional protein-protein interaction network and network analysis, we found some transcription factor like **RELA**, **MYC**, **HIF1a**, **SMAD3** or **REL** that might be the main actor of the acquired resistance and **STAT3**/**STAT4** for the primary resistance.

1.  1. About **RELA**, also called **NF-kB**, it has been shown that it is involved in several cellular functions in hematological malignancies, i.e. inflammation, apoptosis, cell survival, proliferation, angiogenesis, and innate and acquired immunity (ref: [NF-κB pathways in hematological malignancies](https://pubmed.ncbi.nlm.nih.gov/24419302/)).
    2. But also linked to **MYC** as a **NF-kB** repressing factor (**NKRF**), shifted to the cytosplasm by a protein synthesis inhibitor **HTT**, attenuates the transactivation activity of **p65** on the **MYC** gene (ref: [Homoharringtonine deregulates MYC transcriptional expression by directly binding NF-κB repressing factor](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6369765/)).
    3. In AML, **RELA** targets a long noncoding RNA **uc002jit.1**. Its KO impaired the stability of **PARP1**, a prtein that help DNA damage repair. Its KO also inhibit AML cells proliferation and increased the sensitivity to drugs. (ref: [The role of the novel LincRNA uc002jit.1 in NF-kB-mediated DNA damage repair in acute myeloid leukemia cells](https://www-sciencedirect-com.proxy.insermbiblio.inist.fr/science/article/pii/S0014482720302007))
    4. The cleavage of **p65**/**RELA** by **RIP3** induces apoptosis. **RIP3** silencing in leukemia cells results in suppression of the complex regulation of the apoptosis/necroptosis switch and **NF-κB** activity. (ref: [RIP3 is downregulated in human myeloid leukemia cells and modulates apoptosis and caspase-mediated p65/RelA cleavage](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4454320/))
    5. Stabilization of NF-κB-inducing kinase (**NIK**) suppresses AML. Stabilization of NIK-induced activation of NF-κB non-canonical signaling upregulates **Dnmt3a** and downregulates **Mef2c**, which suppresses and promotes AML development, respectively. [Stabilization of NF-κB-Inducing Kinase Suppresses MLL-AF9-Induced Acute Myeloid Leukemia](https://pubmed.ncbi.nlm.nih.gov/29320732/)

2.  1. [A Myc enhancer cluster regulates normal and leukaemic haematopoietic stem cell hierarchies](https://drive.google.com/file/d/1Z0EUt4SEl1aKonIxNyraw3rbo27JNyu6/view?usp=sharing)


## Relevant analyses

To investigate deeper in those results, we need to look into the epigenomic features that might explain that dysregulation of RELA and MYC by in vivo experiments and in vitro analysis.

About in vivo analysis, it would be interesting to check in cell lines the protein abundance of RELA as NF-κB before and after the treatment to IDH inhibitors and before and after chemotherapy with AraC.
In a more large investigation, an overview of the epigenomics features induced by resistance to IDH inhibitors may highlight some mechanisms behind those dysregulation.

The idea of this project around TF analysis is to find the mechanisms that is behind the dysregulation from an epigenomic point of view to a metabolic shift or metabolism adaptability and resistance.

There is different points that can lead to a TF dysregulation. From the accessibility of the chromatin analysed by ATACseq, to the effective binding of a specific TF targeted by a CHiPseq data. The DNA methylation can also be analyzed as it can explain the effective or absence of TF binding to enhancer/promoter in a specific phenotype.

In this paper: [Identification of transcription factor binding sites using ATAC-seq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2), they developed a framework that uses open chromatin data to identify the active transcription factor binding sites. Their method is originally proposed to model the active binding sites by simultaneous analysis of DNase-seq and the ChIP-seq profiles of histone modifications on a genome-wide level.

In this one: [methyl-ATAC-seq measures DNA methylation at accessible chromatin](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6581052/) they present a methyl-ATAC-seq, which implements modifications to ATAC-seq, including subjecting the output to BS-seq. Merging these assays into a single protocol identifies the locations of open chromatin and reveals, unambiguously, the DNA methylation state of the underlying DNA. Such combinatorial methods eliminate the need to perform assays independently and infer where features are coincident.

1. In a huge paper, they did a network genomic integration of phenotypic, structural, and functional relationships. [Metabolic resilience is encoded in genome plasticity](https://www.biorxiv.org/content/10.1101/2021.06.25.449953v2) they showed that a lot of epigenome features are linked to the metabolism.
2. By combining Chip-seq of VDR, FAIRE-seq and RNAseq, they highlighted signatures of VDR pathway [A hierarchical regulatory network analysis of the vitamin D induced transcriptome reveals novel regulators and complete VDR dependency in monocytes](https://pubmed.ncbi.nlm.nih.gov/33753848/)
