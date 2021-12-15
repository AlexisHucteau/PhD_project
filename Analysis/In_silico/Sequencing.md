# Sequencing from scratch

## I - Select datasets

## II - Downloading sequences

> (base) <span style="color:green">alexis@PC-bioinfo</span>:<span style="color:blue">/media/alexis/DATA/KoichiRNAseq</span>$ prefetch --option-file ~/tmp/SRR_Acc_List.txt -O /media/alexis/DATA/KoichiRNAseq/

## III - Extract fastq from .SRA files

> (base) <span style="color:green">alexis@PC-bioinfo</span>:<span style="color:blue">/media/alexis/DATA/KoichiRNAseq</span>$ fastq-dump SRR13928140/SRR13928140.sra -O ./FASTQ_files/

Read 18581205 spots for SRR13928140/SRR13928140.sra  
Written 18581205 spots for SRR13928140/SRR13928140.sra

## IV - Analysis fastq

> (base) <span style="color:green">alexis@PC-bioinfo</span>:<span style="color:blue">/media/alexis/DATA/KoichiRNAseq</span>$ fastqc FASTQ_files/SRR13928140.fastq

Started analysis of SRR13928140.fastq  
Approx 5% complete for SRR13928140.fastq  
Approx 10% complete for SRR13928140.fastq  
Approx 15% complete for SRR13928140.fastq  
Approx 20% complete for SRR13928140.fastq  
Approx 25% complete for SRR13928140.fastq  
Approx 30% complete for SRR13928140.fastq  
Approx 35% complete for SRR13928140.fastq  
Approx 40% complete for SRR13928140.fastq  
Approx 45% complete for SRR13928140.fastq  
Approx 50% complete for SRR13928140.fastq  
Approx 55% complete for SRR13928140.fastq  
Approx 60% complete for SRR13928140.fastq  
Approx 65% complete for SRR13928140.fastq  
Approx 70% complete for SRR13928140.fastq  
Approx 75% complete for SRR13928140.fastq  
Approx 80% complete for SRR13928140.fastq  
Approx 85% complete for SRR13928140.fastq  
Approx 90% complete for SRR13928140.fastq  
Approx 95% complete for SRR13928140.fastq  
Analysis complete for SRR13928140.fastq  
