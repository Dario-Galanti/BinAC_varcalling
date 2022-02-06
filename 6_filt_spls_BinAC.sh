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
### Run: qsub -q short -F "/beegfs/work/bbmdg01/GATK_v3_unfilt_VCF/Samples_to_keep.args" 6_filt_spls_BinAC.sh
### Dependencies: gatk4 v4.1.8.1

## First step of the filtering process. First we select relevant samples, in order to have a more accurate positions filtering.
## This allows to run downsteam analysis with different sets of samples without repeating the variant calling

## ARGUMENT: Samples_to_keep.args --> File in simple text format with one sample name per line (eg. Bioinformatics/WGS_align_SNP_calling/Samples/GWAspls.args)

## Define tools input and output
work=/beegfs/work/bbmdg01
gatk=~/miniconda3/envs/bwa/bin/gatk

## Define input and output
Spls_to_keep=$1					          # File in simple text format with samples to keep (one per line) and .args extention
wDir=${work}/GATK_v3_unfilt_VCF
fin=${wDir}/Ta_v3_vrts_unfilt.vcf.gz
fout=${wDir}/Ta_v3_vrts_unfilt_$(basename $Spls_to_keep .args).vcf.gz

# Index vcf and Filter samples
${gatk} IndexFeatureFile -I $fin
${gatk} --java-options "-Xmx30g" SelectVariants -V $fin -O $fout -sn $Spls_to_keep --QUIET

## Print summary
fin_spls=$(zcat $fin | head -4000 | grep \#CHROM | cut -f10- | tr "\t" "\n" | wc -l)		#The "head" is just making it faster
fout_spls=$(zcat $fout | head -4000 | grep \#CHROM | cut -f10- | tr "\t" "\n" | wc -l)		#The "head" is just making it faster

echo $fin_spls samples were present in the original file. $fout_spls were kept in the filtered file


