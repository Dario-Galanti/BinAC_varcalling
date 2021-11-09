#PBS -l nodes=1:ppn=4 #Nodes and cores
#PBS -l walltime=48:00:00
#PBS -l mem=30Gb 
#PBS -S /bin/bash
#PBS -N GenotypeGVCFs #Job name
#PBS -j oe
#PBS -q short
#PBS -o /beegfs/work/bbmdg01/logs/GatherVCFs.out
#PBS -e /beegfs/work/bbmdg01/logs/GatherVCFs.err

### Aim: Produces a set of joint-called SNP and indel calls ready for filtering (multisample vcf) from the GenomicsDB database produced by GenomicsDBImport in the previous step
### Author: Dario Galanti Dec 2020
### Documentation (Best practices workflow): https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
### Documentation (GatherVcfs): https://gatk.broadinstitute.org/hc/en-us/articles/360037422071-GatherVcfs-Picard-
### Run: qsub -q short 5_GatherVcfs_BinAC.sh
### Dependencies: gatk4 v4.1.8.1

## NB: PROBLEM: In GATK4, combining GVCFs is a two steps process requiring to run 1) GenomicsDBimport and 2) GenotypeGVCFs.
## The first one can be parallelized and is quick, but GenotypeGVCFs cannot be parallelized and loads the whole DB into the RAM (even running only 1 interval with -L).
## My DB is havier than RAM available in BinAC, so everything bogs down and it's super slow!
## The only possible option with my data is to parallelize by chromosome right from GenomicsDBimport.
## Per-chromosome vcf files produced by GenotypeGVCFs have to be combined with GatherVcfs!!!

## Define input and output
work=/beegfs/work/bbmdg01
java=~/miniconda3/envs/bwa/bin/java
inDir=${work}/GATK_v3_scaff_VCFs
outDir=${work}/GATK_v3_unfilt_VCF
fout=${outDir}/Ta_v3_vrts_unfilt.vcf.gz
#tmp=${work}/tmp

mkdir -p $outDir

## NB: PROBLEM: Scaffolds are not ordered correctly in the ref genome!!! So the input file list has to be retrieved from a vcf header and not directly from the file names
## PROBLEM: In addition I excluded scaffolds < 10kb long, so I also have to delete those from the list
# Create input string in format -I chr1.vcf -I chr2.vcf ...
#fin_str=$(ls ${inDir}/*_variants.vcf.gz | sort -V | sed ':a;N;$!ba;s/\n/ -I /g')
#fin_str=$(zcat ${inDir}/Scaffold_4569_v3_Ta_variants.vcf.gz | grep 'contig=<ID' | awk -v inDir="${inDir}" -F'[=,]' '{print inDir"/"$3"_variants.vcf.gz"}' | sed ':a;N;$!ba;s/\n/ -I /g')
used_scaff=( $(ls ${inDir}/*_variants.vcf.gz | sort -V ) )
ordered_scaff=( $(zcat $(ls ${inDir}/*.vcf.gz | head -1) | grep 'contig=<ID' | awk -v inDir="${inDir}" -F'[=,]' '{print inDir"/"$3"_variants.vcf.gz"}') )

# Intersect arrays to delete unused scaffolds but keep genome scaffold order
fin_arr=()
for i1 in "${ordered_scaff[@]}"; do
	for i2 in "${used_scaff[@]}"; do
		if [[ $i1 = $i2 ]]; then
		fin_arr+=("$i1")
        fi
    done
done
fin_str=$(printf "%s\n" "${fin_arr[@]}" | sed ':a;N;$!ba;s/\n/ -I /g')

# It's a picard tool, I should probably use picard syntax
${java} -Xmx30g -jar ~/miniconda3/envs/bwa/share/picard-2.23.6-0/picard.jar GatherVcfs -I ${fin_str} -O ${fout} --QUIET
#~/miniconda3/envs/bwa/bin/gatk --java-options \"-Xmx30g\" GatherVcfs  -I ${fin_str} -O ${fout} --tmp-dir ${tmp} --QUIET


