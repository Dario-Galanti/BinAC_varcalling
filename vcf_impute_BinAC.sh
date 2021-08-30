#PBS -l nodes=1:ppn=8 #Nodes and cores
#PBS -l walltime=20:00:00
#PBS -l mem=20gb 
#PBS -S /bin/bash
#PBS -N imputation #Job name
#PBS -j oe
#PBS -q short
#PBS -o /beegfs/work/bbmdg01/logs/beagle_impute.out
#PBS -e /beegfs/work/bbmdg01/logs/beagle_impute.err

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
fin=${work}/Ta_v3_10NA_1MAF/Ta_v3_vrts_10NA_1MAF_biall_GWAspls.vcf.gz

outDir=${work}/Ta_v3_10NA_1MAF
fout_base=${outDir}/Ta_v3_vrts_1MAF_imputed_GWAspls
fout=${fout_base}.vcf.gz
base_withref=${outDir}/Ta_v3_vrts_1MAF_imputed_GWAspls_withref
vcf_withref=${base_withref}.vcf.gz

mkdir -p $outDir

## 1) Run BEAGLE
$java -Xss5m -Xmx20g -jar $beaglejar nthreads=8 gt=$fin out=$fout_base

## 2) ADD REFERENCE GENOTYPE TO THE VCF FILE
Ref=MN106A
zcat $fout | awk -v Ref="$Ref" '{if(substr($1,1,2)=="##"){print}else if(substr($1,1,2)=="#C"){OFS="\t";print $0,Ref}else{OFS="\t";print $0,"0|0"}}' | $bgzip -c > $vcf_withref

## OPTIONAL: INDEX VCF FILES
~/miniconda3/envs/bwa/bin/tabix -p vcf $fout
~/miniconda3/envs/bwa/bin/tabix -p vcf $vcf_withref

## OPTIONAL: RECODE VCF IN PLINK FORMAT FOR GWAS
## IMPORTANT!!! NB: THIS REQUIRES SOME RELEVANT FORMATTING DESCRIBED BELOW.
## 1) Substitute "_" with "-" in sample names or plink crashes
## 2) Add variant name -> "chr_pos"
## 3) IMPORTANT!!! Plink is very very bad at handling non-model species with Scaffold-level assemblies.
## 3) We use --allow-extra-chr, but this causes eg. Chr1 --> 1 but Scaffold_1 --> Scaffold_1.
## 3) SOLUTION: AFTER running plink we change all scaffold names in the map file with the following rationale:
## 3) SOLUTION: Scaffold_n --> n+7 (eg. Scaffold_1 --> 8) in order not to re-use Chr numbers

zcat $vcf_withref | awk 'OFS="\t"{if(!/\#/){$3=$1"_"$2;print} else {gsub("_","-");print}}' > ${outDir}/temp_ready.vcf
$plink --vcf ${outDir}/temp_ready.vcf --allow-extra-chr --recode --out $base_withref
rm ${outDir}/temp_ready.vcf
#$plink --file $base_withref --allow-extra-chr --make-bed --out $base_withref
$plink --file $base_withref --allow-extra-chr --recode A --out $base_withref
awk 'OFS="\t"{if ($1 ~ /^Scaffold/){gsub(/Scaffold_/,"",$1);$1=($1+7)};print}' ${base_withref}.map > ${outDir}/map_ready.map
mv ${outDir}/map_ready.map ${base_withref}.map

## OPTIONAL: Calculate IBD matrix for later use (pop structure correction)
#$plink --file $base_withref --allow-extra-chr --distance square ibs flat-missing
