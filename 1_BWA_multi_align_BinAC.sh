#!/bin/bash

### Aim: WGS read alignment with BWA (BWA manual: http://bio-bwa.sourceforge.net/bwa.shtml) and Mark Duplicates with gatk MarkDuplicatesSpark
### Author: Dario Galanti Sept 2020
### Run: bash 1_BWA_multi_align_BinAC.sh
### Dependencies: bwa and gatk4 v4.1.8.1

## NB: This script is made for samples sequenced on a single lane. For samples sequenced on multiple lanes the "Extract read group" should be double checked

### PRE-STEPS:
### 1) Prepare conda environment
## conda create -n bwa
## conda install -c bioconda bwa
## conda install -c bioconda gatk4 (if issues with sometools check https://github.com/bioconda/bioconda-recipes/issues/12100 and use user's "Leipzig" solution)
### 2) Index reference genome (-p prefix)
## ~/miniconda3/envs/bwa/bin/bwa index -p Ta_genome /beegfs/work/bbmdg01/Tarvense_genome/v4_thlaspi_final/final.fasta
### 3) Check number of reads per sample
## mkdir fastqc
## echo -e Sample"\t"Mate"\t"Num_reads > Total_reads.txt
## for f in ../WGS_trimmed/trim.TA*;do spl=$(echo $f | rev | cut -c12-34 | rev); mate=$(echo $f | rev | cut -c10); reads=$(zcat $f | grep -c @);echo -e $spl"\t"$mate"\t"$reads >> Total_reads.txt;done

## Define home dir (work), tools, input and output
work=/beegfs/work/bbmdg01
bwa=~/miniconda3/envs/bwa/bin/bwa
gatk=~/miniconda3/envs/bwa/bin/gatk
genome=${work}/Tarvense_genome/v4_thlaspi_final/Ta_genome	# Basename of bwa reference genome index (see PRE-STEPS)
inDir=${work}/WGS_trimmed
outDir=${work}/BWA_alignments

## Make directories for individual job scripts, logs and output
mkdir -p ${outDir}
mkdir -p ${work}/work/bwa
mkdir -p ${work}/logs/bwa
mkdir -p ${work}/BWA_align_metrics	# Dir for alignment metrics
mkdir -p ${work}/tmp				# Temporary directory for MarkDuplicatesSpark
#mkdir -p ${work}/java_tmp			# Temporary directory for SortSam

## SUBMIT MULTIPLE ALIGNMENTS
for file in ${inDir}/trim.TA_SE_06_*_1.fastq.gz;
do
	sample=$(basename $file _1.fastq.gz | cut -d"." -f2)
	jobName=${work}/work/bwa/bwa.${sample}.sh
	(
	echo "#PBS -l nodes=1:ppn=8 #Nodes and cores"
	echo "#PBS -l walltime=48:00:00"
	echo "#PBS -l mem=12Gb"
	echo "#PBS -S /bin/bash"
	echo "#PBS -N bwa.${sample}"
	echo "#PBS -j oe"
	echo "#PBS -q short"
	echo "#PBS -o ${work}/logs/bwa/bwa.${sample}.out"
	echo "#PBS -e ${work}/logs/bwa/bwa.${sample}.err"
	echo ""
	
	##echo "cd ${work}/Tarvense_genome/v4_thlaspi_final/"
	
	## Define input and output files
	echo "fin_1=${inDir}/$(basename $file)"
	echo "fin_2=${inDir}/$(basename $file 1.fastq.gz)2.fastq.gz"
	echo "samfile=${outDir}/${sample}.sam"
	#echo "sort_sam=${work}/BWA_alignments/sort_${sample}.sam"
	echo "bamfile=${outDir}/${sample}.bam"
	echo "metrics=${work}/BWA_align_metrics/dups_metrics_${sample}.txt"
	echo ""
	## Extract read group ID (this is necessary for GATK!!! But NB: It only works for my dataset where each sample was sequenced on a single lane!!!).
	echo "header=\$(zcat \${fin_1} | head -n 1)"
	echo "id=\$(echo \${header} | cut -f 1-4 -d: | sed 's/@//' | sed 's/:/_/g')"
	echo "sm=\$(echo \${header} | grep -Eo \"[ATGCN]+\$\")"
	echo "RG=\"@RG\tID:\${id}\tSM:${sample}\tLB:\${id}_\${sm}\tPL:ILLUMINA\""
	echo ""
	## Run bwa alignment
	echo "${bwa} mem -t 8 -R $(echo \${RG}) -c 1 ${genome} \${fin_1} \${fin_2} > \${samfile}"	# -R attaches read group, necessary for MarkDuplicatesSpark
	echo "echo -e \"${sample} bwa alignment finished\n\""
	## Sort by queryname (this step may not be necessary and it causes some samples to run out of disk space, neet to specify /tmp if willing to use it)
	#echo "${gatk} --java-options \"-Djava.io.tmpdir=${work}/java_tmp -Xmx8G\" SortSam -I \${samfile} -O \${sort_sam} -SO queryname --QUIET"
	#echo "echo -e \"${sample} SortSam finished\n\""
	## MarkDuplicatesSpark (will output a coordinate-sorted bamfile)
	echo "${gatk} MarkDuplicatesSpark -I \${samfile} -O \${bamfile} -M \${metrics} --allow-multiple-sort-orders-in-input --conf 'spark.executor.cores=8' --conf 'spark.local.dir=${work}/tmp' --QUIET"
	#echo "rm \${sort_sam}"
	echo "rm \${samfile}"
	echo "[ -f \${bamfile} ] && echo sample ${sample} bamfile is ready for HaplotypeCaller"
	) > $jobName
		chmod +x $jobName
		echo "Bash file $jobName created"
		qsub -q short ${jobName}

done
exit



