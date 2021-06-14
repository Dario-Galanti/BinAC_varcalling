#PBS -l nodes=1:ppn=4 #Nodes and cores
#PBS -l walltime=20:00:00
#PBS -l mem=30Gb 
#PBS -S /bin/bash
#PBS -N filter_spls #Job name
#PBS -j oe
#PBS -q short
#PBS -o /beegfs/work/bbmdg01/logs/6_filt_spls.out
#PBS -e /beegfs/work/bbmdg01/logs/6_filt_spls.err

### Aim: Subset samples from multisample vcf file using GATK 4.1
### Author: Dario Galanti Jan 2021
### Documentation (Best practices workflow): https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
### Documentation (SelectVariants): https://gatk.broadinstitute.org/hc/en-us/articles/360036365752-SelectVariants#--exclude-sample-name
### Run: qsub -q short 6_filt_spls_BinAC.sh
### Dependencies: gatk4 v4.1.8.1

## First step of the filtering process. First we select relevant samples, in order to have a more accurate positions filtering.
## This allows to run downsteam analysis with different sets of samples without repeating the variant calling

## Define input and output
work=/beegfs/work/bbmdg01
gatk=~/miniconda3/envs/bwa/bin/gatk

wDir=${work}/GATK_v3_unfilt_VCF
fin=${wDir}/Ta_v3_vrts_unfilt.vcf.gz
fout=${wDir}/Ta_v3_vrts_unfilt_GWAspls.vcf.gz
delSPLS_file=${wDir}/spls_to_delete.args
#tmp=${work}/tmp

## SAMPLES TO DELETE: Here I only delete samples for which no information for any GWAS analysis is available (no WGBS, no CG_GH_2019, no CG_2020).
delSPLS=(TA_AM_02_01_F3_CC0_M1_1 TA_AM_03_01_F3_CC0_M1_1 TA_AM_04_01_F3_CC0_M1_1 TA_AM_05_01_F3_CC0_M1_1 TA_AM_06_01_F3_CC0_M1_1 TA_AM_07_01_F3_CC0_M1_1 TA_DE_17_01_F3_CC0_M1_1 TA_DE_17_02_F3_CC0_M1_1 TA_DE_18_01_F3_CC0_M1_1 TA_DE_19_01_F3_CC0_M1_1 TA_FR_04_01_F3_CC0_M1_1)
# store samples to delete to file
for i in ${delSPLS[@]}; do echo $i >> $delSPLS_file ;done

# Index vcf and Filter samples
${gatk} IndexFeatureFile -I $fin
${gatk} --java-options "-Xmx30g" SelectVariants -V $fin -O $fout -xl-sn $delSPLS_file --QUIET

## Print summary
fin_spls=$(zcat $fin | head -3000 | grep \#CHROM | cut -f10- | grep -o TA_ | wc -l)		#The "head" is just making it faster
fout_spls=$(zcat $fout | head -3000 | grep \#CHROM | cut -f10- | grep -o TA_ | wc -l)		#The "head" is just making it faster

echo $fin_spls samples were present in the original file. $fout_spls were kept in the filtered file



