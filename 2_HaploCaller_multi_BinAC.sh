#!/bin/bash

### Aim: Run GATK Haplotypecaller to obtain individual samples GVCF files
### Author: Dario Galanti Nov 2020
### Documentation (Best practices workflow): https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
### Documentation (HaplotypeCaller): https://gatk.broadinstitute.org/hc/en-us/articles/360046787532-HaplotypeCaller
### Run: bash 2_HaploCaller_multi_BinAC.sh
### Dependencies: gatk4 v4.1.8.1

### PRE-STEPS:
### Make reference genome dictionary
###	Dependencies: gatk4
### ~/miniconda3/envs/bwa/bin/java -jar ~/miniconda3/envs/bwa/share/picard-2.23.6-0/picard.jar CreateSequenceDictionary -R final.fasta

## Define home dir (work), tools, input and output
work=/beegfs/work/bbmdg01
gatk=~/miniconda3/envs/bwa/bin/gatk
ref=${work}/Tarvense_genome/v3_thlaspi_clr/thlaspi.fa
inDir=${work}/BWA_alignments
outDir=${work}/GATK

## Make directories for individual job scripts, logs and output
mkdir -p ${outDir}
mkdir -p ${work}/work/haplocall
mkdir -p ${work}/logs/haplocall
#mkdir -p ${work}/tmp				# Temporary directory

# RUN GATK HaplotypeCaller FOR MULTIPLE SAMPLES
for file in ${inDir}/TA_DE_01*_1.bam;
do
	sample=$(basename $file .bam)
	jobName=${work}/work/haplocall/hapcall.${sample}.sh
	(
	echo "#PBS -l nodes=1:ppn=10 #Nodes and cores"
	echo "#PBS -l walltime=20:00:00"
	echo "#PBS -l mem=40Gb"
	echo "#PBS -S /bin/bash"
	echo "#PBS -N hapcall.${sample}"
	echo "#PBS -j oe"
	echo "#PBS -q short"
	echo "#PBS -o ${work}/logs/haplocall/hapcall.${sample}.out"
	echo "#PBS -e ${work}/logs/haplocall/hapcall.${sample}.err"
	echo ""
	
	echo "cd ${work}"
	## Define input and output files
	echo "fout=${outDir}/$(basename $file).vcf.gz"
	##echo "ref=${work}/Tarvense_genome/v3_thlaspi_clr/thlaspi.fa"
	echo ""
	
	## Run HaplotypeCaller (https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller#--emit-ref-confidence)
	## We run in GVCF mode: https://gatk.broadinstitute.org/hc/en-us/articles/360035531812?id=4017
	## Briefly GVCF is a vcf with some extra mapping quality info either per bp or per variant (as in our command)
	#echo "${gatk} --java-options \"-Xmx40g\" HaplotypeCaller -R ${ref} -I ${file} -O ${fout} -ERC GVCF -G Standard -G AS_Standard --QUIET"
	echo "${gatk} --java-options \"-Xmx40g\" HaplotypeCaller -R ${ref} -I ${file} -O \${fout} -ERC GVCF --QUIET"
	) > $jobName
		chmod +x $jobName
		echo "Bash file $jobName created"
		qsub -q short ${jobName}

done
exit

