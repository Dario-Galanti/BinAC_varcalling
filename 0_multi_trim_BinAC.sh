#!/bin/bash

### Aim: WGS read trimming with cutadapt (cutadapt manual: https://cutadapt.readthedocs.io/en/stable/)
### Author: Dario Galanti 2019
### Run: bash 0_multi_trim_BinAC.sh
### Dependencies: cutadapt

## Define indir, outdir (work) and tools
work=/beegfs/work/bbmdg01
cutadapt=~/miniconda3/envs/cutadapt/bin/cutadapt
inDir=${work}/WGS_US
outDir=${work}/WGS_UStrim

## Make directories for individual job scripts and logs
mkdir -p ${outDir}
mkdir -p ${work}/work/cutadapt
mkdir -p ${work}/logs/cutadapt

## RUN MULTIPLE TRIMMINGS
for file in ${inDir}/*_1.fastq.gz;
do
	jobName=${work}/work/cutadapt/trim."$(basename $file _1.fastq.gz)".sh
	(
	echo "#PBS -l nodes=1:ppn=4 #Nodes and cores"
	echo "#PBS -l walltime=06:00:00"
	echo "#PBS -S /bin/bash"
	echo "#PBS -N trim."$(basename $file _1.fastq.gz)""
	echo "#PBS -j oe"
	echo "#PBS -q short"
	echo "#PBS -o ${work}/logs/cutadapt/trim."$(basename $file _1.fastq.gz)".out"
	echo "#PBS -e ${work}/logs/cutadapt/trim."$(basename $file _1.fastq.gz)".err"

	echo "fin_1=${inDir}/$(basename $file)"
	echo "fin_2=${inDir}/$(basename $file 1.fastq.gz)2.fastq.gz"
	echo "fout_1=${outDir}/trim.$(basename $file)"
	echo "fout_2=${outDir}/trim.$(basename $file 1.fastq.gz)2.fastq.gz"

	echo "${cutadapt} -j 4 -q 20 --minimum-length 25 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -o \"\$fout_1\" -p \"\$fout_2\" \"\$fin_1\" \"\$fin_2\""
	) > $jobName
		chmod +x $jobName
		echo "Bash file $jobName created"
		qsub -q short ${jobName}

done
exit

