#!/bin/sh

## Author: Dario Galanti
## Aim: Find scaffolds harbouring 1 or few position in vcf file and remove them from the vcf file
## Input: vcf.gz file
## Dependencies: vcftools
## Run: bash lonely_vcf_pos.sh <input.vcf.gz> <output.vcf.gz> 
## Run: bash lonely_vcf_pos.sh SNPs_v2_ewas_filt/new_Ta_v2_ewas_filt.vcf.gz new_Ta_v2_ewas_filt.vcf.gz 

# Define tool paths
vcftools=~/miniconda3/envs/vcftools/bin/vcftools
bgzip=~/miniconda3/envs/bwa/bin/bgzip

#Retrieve scaffolds harbouring <x positions and save them in a txt file
x=3										#Minimum number of pos per scaffold to keep the scaffold in the vcf
zcat $1 | grep -v \# | awk '{print $1"\t"}' | uniq -c | awk -v x=$x '$1 < x {print $2"\t"}'> lonely_pos1.txt
zcat $1 | grep -v \# | grep -f lonely_pos1.txt | awk '{print $1"\t"$2}' > lonely_pos.txt
rm lonely_pos1.txt

#Filter out the lonely_positions
#vcftools --gzvcf $1 --exclude-positions lonely_pos.txt --recode --recode-INFO-all --stdout | bgzip -c > $2
$vcftools --gzvcf $1 --exclude-positions lonely_pos.txt --recode --recode-INFO-all --stdout | $bgzip -c > $2