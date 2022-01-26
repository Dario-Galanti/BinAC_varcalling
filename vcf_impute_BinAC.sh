#PBS -l nodes=1:ppn=8 #Nodes and cores
#PBS -l walltime=40:00:00
#PBS -l mem=40gb 
#PBS -S /bin/bash
#PBS -N imputation #Job name
#PBS -j oe
#PBS -q short
#PBS -o /beegfs/work/bbmdg01/logs/beagle_impute_v4.out
#PBS -e /beegfs/work/bbmdg01/logs/beagle_impute_v4.err

### Aim: Phase and impute missing GT calls from vcf file using beagle 5.1; Adding ref to the vcf and formatting for GWAS
### Author: Dario Galanti Jan 2021
### Documentation: https://faculty.washington.edu/browning/beagle/old.beagle.html
### Input: Already filtered vcf file containing a minority of missing GT calls (./.)
### Dependencies: conda install -c bioconda beagle
### Run: qsub -q short vcf_impute_BinAC.sh

### IMPORTANT: Strict MAF filtering (0.05) should be applied AFTER imputation with BEAGLE to increase imputation accuracy as found in the paper below!
### "Impact of pre-imputation SNP-filtering on genotype imputation results"

## Define tools
java=~/miniconda3/envs/vcftools/bin/java
beaglejar=~/miniconda3/envs/vcftools/share/beagle-5.1_24Aug19.3e8-0/beagle.jar
bgzip=~/miniconda3/envs/bwa/bin/bgzip
plink=~/miniconda3/envs/vcftools/bin/plink

## Define input and outputs
work="/beegfs/work/bbmdg01"
#fin=${work}/Ta_v3_10NA_1MAF/Ta_v3_vrts_10NA_1MAF_biall_GWAspls.vcf.gz
fin=${work}/GATK_v4_filt_VCF/Ta_v4_vrts_10NA_biall_1MAF_withCNspls.vcf.gz

outDir=${work}/Ta_v4_10NA_imp
fout_base=${outDir}/Ta_v4_vrts_10NA_biall_1MAF_imp_withCNspls
fout=${fout_base}.vcf.gz
base_withref=${outDir}/Ta_v4_vrts_10NA_biall_1MAF_imp_withCNspls_withref
vcf_withref=${base_withref}.vcf.gz

mkdir -p $outDir

## 1) Run BEAGLE
$java -Xss5m -Xmx20g -jar $beaglejar nthreads=8 gt=$fin out=$fout_base

## 2) ADD REFERENCE GENOTYPE TO THE VCF FILE
Ref=MN106A		#Define name of the reference
zcat $fout | awk -v Ref="$Ref" '{if(substr($1,1,2)=="##"){print}else if(substr($1,1,2)=="#C"){OFS="\t";print $0,Ref}else{OFS="\t";print $0,"0|0"}}' | $bgzip -c > $vcf_withref

## OPTIONAL: INDEX VCF FILES
~/miniconda3/envs/bwa/bin/tabix -p vcf $fout
~/miniconda3/envs/bwa/bin/tabix -p vcf $vcf_withref


