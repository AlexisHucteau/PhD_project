# Introduction

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
    2.


3.  1. **HIF1a**
    2. test


4.  1. **SMAD3**
    2. test


5.  1. **REL**
    2. test


6.  1. **STAT3**
    2. test


7.  1. **STAT4**
    2. test


## Relevant analyses

1. In a huge paper, they did a network genomic integration of phenotypic, structural, and functional relationships. [Metabolic resilience is encoded in genome plasticity](https://www.biorxiv.org/content/10.1101/2021.06.25.449953v2) they showed that a lot of epigenome features are linked to the metabolism.
2. By combining Chip-seq of VDR, FAIRE-seq and RNAseq, they highlighted signatures of VDR pathway [A hierarchical regulatory network analysis of the vitamin D induced transcriptome reveals novel regulators and complete VDR dependency in monocytes](https://pubmed.ncbi.nlm.nih.gov/33753848/)
3. Splicing biblio & Rmats
4. Integrated Stress response biblio
5.

## Reflections on the topic

The main topic of the thesis is the resistance of refractory/relapsed AML patients cells to IDH inhibitor by metabolic plasticity and/or stress response through epigenomic changes.

In a bioinformatic point of view, it's important to know the tools and the data that permit to dig into their regulations. Transcriptions data are the first part of the analysis as it permits to find transcriptional regulation and by different tools, it's possible to infer the activity of proteins involved in the transcription (like Transcription factor activity), metabolic analysis, deconvolution. But transcriptions data don't exhibit regulation downstream regulations like translational regulation, splicing event neither upstream regulations like epigenetic features. In this project where there are many epigenetic changes, looking at upstream regulations can exhibit what trnascriptional data can't or can explain some underestimated gene expression deregulations and make them part of the main mechanism.

In a bioinformatic point of view, it's important to know the tools and the data that permit to investigate those epigenetic changes. Epigenetic in IDH is characterised, as we said, in a DNA hypermethylation phenotype. DNA methylation can be catched by Bisulfite Sequencing through CpGs methylation sequencing. It permits to find promoter hyper/hypomethylation that can interfer the binding of TFs on DNA and block the transcription of genes. We can notice that the hyper/hypomethylation of a promoter is directly linked to the transcription and if the gene is not up/down expressed, it may be not relevant to notice it. In other terms, DNA methylation in promoter and even in enhancer can't highlight different genes than transcription analysis. It only permits to highlight mechanisms behind the transcription. The mechanisms that can be found are the mechanisms that explain why some genes are dysregulated.

DNA methylation can appear in enhancer but the data on it is not very robust as sequencing is often focused on promoter and it's not possible for some data to suggest if an enhancer is hyper or hypomethylated. It's possible to bypass the problem by looking at overall methylation in a 3D context by using a 3D chromatin contact network but it need to be specific of the cell and, to me, AML cells that are close to HSC cells have their own chromatin contact network that change due to genomic mutations and alterations. Methylations alteractions can also be directly /indirectly linked to that chromatin modifications (need to read this paper : [DNA Methylation and Chromatin Remodeling: The Blueprint of Cancer Epigenetics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4826949/)). Methylation and chormatin remodeling looks like a feedback loop where we need the two data to understand their effect on transcription. It's important to add histone modifications in the loop too.
