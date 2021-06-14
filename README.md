# BinAC_varcalling
Workflow for read mapping (bwa), variant calling (GATK4) and filtering from short-read Whole Genome Sequencing data, for medium and large datasets of non-model species.
System setup: Scripts are made for running on a linux-based cluster with PBS queueing system.

The workflow is meant for the analysis of paired-end short reads. Based on the GATK4 [best practices for germline short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-), is specifically customized for medium and large WGS datasets of non-model species. Hence, it can handle a large number of samples, large and highly fragmented genomes and the absence of previous variant sets.

BRIEF WORKFLOW DESCRIPTION:
0_multi_trim_BinAC.sh
    Adaptors trimming with [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
1_BWA_multi_align_BinAC.sh
    Single-sample read alignment with [bwa-mem](http://bio-bwa.sourceforge.net/bwa.shtml) and detection of duplicates with [MarkDuplicatesSpark](https://gatk.broadinstitute.org/hc/en-us/articles/360046221811-MarkDuplicatesSpark).
2_HaploCaller_multi_BinAC.sh
    Single-sample variant calling with local reassembly ([Haplotypecaller](https://gatk.broadinstitute.org/hc/en-us/articles/360036715891-HaplotypeCaller)) to obtain to obtain individual sample GVCF files.
3_4_GenDB_GenoGVCFs_BinAC.sh
    
5_GatherVcfs_BinAC.sh
6_filt_spls_BinAC.sh
7_filt_variants_BinAC.sh




