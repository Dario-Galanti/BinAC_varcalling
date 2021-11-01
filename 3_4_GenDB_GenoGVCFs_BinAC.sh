### Aim: Combine GVCF individual files produced by HaplotypeCaller in a two step process: 1) GenomicsDBImport; 2) GenotypeGVCFs
### Author: Dario Galanti Nov 2020
### Documentation (Best practices workflow): https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
### Documentation (GenomicsDBImport): https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport
### Documentation (GenotypeGVCFs): https://gatk.broadinstitute.org/hc/en-us/articles/360036899732-GenotypeGVCFs
### Documentation (Speed-up): https://gatk.broadinstitute.org/hc/en-us/community/posts/360063088471-Speeding-up-GenotypeGVCFS-GATK4
### Run: bash 3_4_GenDB_GenoGVCFs_BinAC.sh
### Dependencies: gatk4 v4.1.8.1

## NB: PROBLEM: In GATK4, combining GVCFs is a two steps process requiring to run 1) GenomicsDBimport and 2) GenotypeGVCFs.
## The first one can be parallelized and is quick, but GenotypeGVCFs cannot be parallelized and loads the whole DB into the RAM (even running only 1 interval with -L).
## My DB is havier than RAM available in BinAC, so everything bogs down and it's super slow!
## The only possible option with my data is to parallelize by chromosome right from GenomicsDBimport.
## Per-chromosome vcf files produced by this script have to be combined with GatherVcfs!!!

## RESOURCE ALLOCATION: Chr 1,2,3 -> 32Gb, 4 cpus, 100h | Chr 4,5,6 -> 30Gb , 4 cpus, 48h | Scaffolds -> 20Gb, 2 cpus, 48h

## Define input and output
work=/beegfs/work/bbmdg01
gatk=~/miniconda3/envs/bwa/bin/gatk
genome=${work}/Tarvense_genome/v3_thlaspi_clr/thlaspi.fa
index=${work}/Tarvense_genome/v3_thlaspi_clr/thlaspi.fa.fai
intervals=${work}/Tarvense_genome/v3_thlaspi_clr/GATK_thlaspi.bed
samples=${work}/Samples_map.txt
inDir=${work}/GATK
DBDir=${work}/GenoDB	# The tool will creat the DB folders, but since we are running per chromosome we have to make a folder containing ad DBs
outDir=${work}/GATK_v3_scaff_VCFs
tmp=${work}/tmp

## Make directories for individual job scripts, logs and output
mkdir -p $DBDir
mkdir -p $outDir
mkdir -p ${work}/work/GenDB_GenoGVCF
mkdir -p ${work}/logs/GenDB_GenoGVCF

## Make intervals file (NB: Only calling SNPs for Scaffolds > 10kb !!!!!!!) and chromosome array
## Intervals: GenomicsDBImport requires intervals to run. See (https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists)
if [ ! -f $intervals ]; then awk 'OFS="\t"{if($2>10000){print $1,0,$2}}' $index > $intervals ;fi
chr_arr=($(cut -f1 $intervals))

## Make samples map file
ls -d $inDir/*vcf.gz | awk '{OFS="\t"; print substr($0,27,23), $0}' > $samples

for chr in ${chr_arr[*]};
do
	# NB: We run only a portion of the scaffolds as larger scaffolds require more memory allocation then the smaller ones
	chr_num=$(echo $chr | grep -Eo '[0-9]{1,6}')
	if [ $chr_num -ge 1 ] && [ $chr_num -le 7 ];		#Define range of scaffolds to run (includes range extremes). Larger scaffolds require more memory
	then
		jobName=${work}/work/GenDB_GenoGVCF/GATK4.${chr}
		(
		echo "#PBS -l nodes=1:ppn=4 #Nodes and cores"
		echo "#PBS -l walltime=48:00:00"
		echo "#PBS -l mem=32Gb"
		echo "#PBS -S /bin/bash"
		echo "#PBS -N GATK4.${chr}"
		echo "#PBS -j oe"
		echo "#PBS -q short"
		echo "#PBS -o ${work}/logs/GenDB_GenoGVCF/GATK4.${chr}.out"
		echo "#PBS -e ${work}/logs/GenDB_GenoGVCF/GATK4.${chr}.err"
		echo ""
		
		echo "cd ${work}"
		
		echo "chr_DBDir=${DBDir}/${chr}"	# GenomicsBDimport will make this directory
		## Run GATK4 GenomicsDBImport
		echo "~/miniconda3/envs/bwa/bin/gatk --java-options \"-Xmx20g\" GenomicsDBImport --genomicsdb-workspace-path \${chr_DBDir} -L ${chr} --sample-name-map ${samples} --tmp-dir ${tmp} --reader-threads 4 --QUIET"
		## Run GATK4 GenotypeGVCFs
		echo "chr_fout=${outDir}/${chr}_v3_Ta_variants.vcf.gz"
		echo "${gatk} --java-options \"-Xmx20g\" GenotypeGVCFs -R ${genome} -V gendb://\${chr_DBDir} -O \${chr_fout} --tmp-dir ${tmp} --QUIET"
		) > $jobName
			chmod +x $jobName
			echo "Bash file $jobName created"
			qsub -q short ${jobName}
	
	fi
done
exit

