#PBS -l nodes=1:ppn=4 #Nodes and cores
#PBS -l walltime=40:00:00
#PBS -l mem=30Gb 
#PBS -S /bin/bash
#PBS -N filter_variants #Job name
#PBS -j oe
#PBS -q short
#PBS -o /beegfs/work/bbmdg01/logs/7_filt_var.out
#PBS -e /beegfs/work/bbmdg01/logs/7_filt_var.err

### Aim: Filter variants from vcf file with GATK 4.1
### Author: Dario Galanti Jan 2021
### Documentation (Best practices workflow): https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
### Documentation (How to Filter Variants): https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
### Documentation (Hard Filtering Germline Short Variants): https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
### Documentation (Unable to use VQRS): https://gatk.broadinstitute.org/hc/en-us/articles/360037499012-I-am-unable-to-use-VQSR-recalibration-to-filter-variants
### Documentation (VariantFiltration): https://gatk.broadinstitute.org/hc/en-us/articles/360046788072-VariantFiltration#--genotype-filter-expression
### Run: qsub -q short 7_filt_variants_BinAC.sh
### Dependencies: gatk4 v4.1.8.1, R, java, vcftools, samtools(bgzip), lonely_vcf_pos.sh(custom script)

### IMPORTANT: Strict MAF filtering (0.05) should be applied AFTER imputation with BEAGLE to increase imputation accuracy as found in the paper below!
### "Impact of pre-imputation SNP-filtering on genotype imputation results"

## Process
## 1) Split SNPs and INDELS/MIXED variants
## 2) Make diagnostic plots, these can be checked to confirm the filtering parameters we used are not too restrictive.
## 3) Filter SNPs and INDELS/MIXED variants separately (VariantFiltration)
## 4) Remove filtered sites from the files (SelectVariants)
## 5) Rejoin SNP and INDEL files and sort (SortVcf)
## 6) Filter for NAs<0.1 (6a) and for NAs<0.1, biallelic sites and MAF (6b).
## 7) Filter scaffolds harbouring less than 3 variants
## 8) Fix vcftools converting phased missing positions ".|." into "."

## Define input and outputs
work=/beegfs/work/bbmdg01
tmp=${work}/tmp

wDir=${work}/GATK_v3_unfilt_VCF
fin=${wDir}/Ta_v3_vrts_unfilt_GWAspls.vcf.gz
SNP_unfilt=${wDir}/Ta_v3_SNPs_unfilt_GWAspls.vcf.gz
INDEL_unfilt=${wDir}/Ta_v3_INDELs_unfilt_GWAspls.vcf.gz

outDir=${work}/GATK_v3_filt_VCF
mkdir -p ${outDir}
SNP_VarFilt=${outDir}/Ta_v3_SNPs_VarFilt_GWAspls.vcf.gz
INDEL_VarFilt=${outDir}/Ta_v3_INDELs_VarFilt_GWAspls.vcf.gz
NAfilt_vcf=${outDir}/Ta_v3_vrts_NAfilt_GWAspls.vcf.gz		# After NAs<0.1 filtering
filtvcf_1=${outDir}/Ta_v3_vrts_filt1_GWAspls.vcf.gz			# After NAs<0.1 biallelic sites and 0.1MAF filtering
filtvcf_2=${outDir}/Ta_v3_vrts_filt2_GWAspls.vcf.gz			# After scaffold filtering
final_vcf=${outDir}/Ta_v3_vrts_filt_GWAspls.vcf.gz			# After fixing "./." to "." vcftools issue

## Define tools
gatk=~/miniconda3/envs/bwa/bin/gatk
Rscript=~/miniconda3/envs/R/bin/Rscript
java=~/miniconda3/envs/bwa/bin/java
vcftools=~/miniconda3/envs/vcftools/bin/vcftools
bgzip=~/miniconda3/envs/bwa/bin/bgzip
lonely_vcf_pos=${work}/scripts/lonely_vcf_pos.sh	# Custom script to filter scaffolds harbouring less than x positions (x=3, can be defined in the script)

## 1a) Subset SNPs
$gatk --java-options "-Xmx30g" SelectVariants -V $fin -select-type SNP -O $SNP_unfilt --QUIET
## 1b) Subset INDELS and MIXED variants
$gatk --java-options "-Xmx30g" SelectVariants -V $fin -select-type INDEL -select-type MIXED -O $INDEL_unfilt --QUIET

## 2) Make diagnostic plots
$gatk --java-options "-Xmx30g" VariantsToTable -V $SNP_unfilt -O ${wDir}/GVCFall_SNPs.table -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR --QUIET
$gatk --java-options "-Xmx30g" VariantsToTable -V $INDEL_unfilt -O ${wDir}/GVCFall_INDELs.table -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR --QUIET
$Rscript --vanilla ${work}/Filtering_diagnostic_plots.R ${wDir}/GVCFall_SNPs.table ${wDir}/GVCFall_INDELs.table

## 3a) Filter genotypes from SNPs file
## Commands taken from: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering#2
## Logic of filtering: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
$gatk --java-options "-Xmx30g" VariantFiltration -V $SNP_unfilt -O ${wDir}/Ta_v3_SNPs_tempfilt_GWAspls.vcf.gz \
 -filter "QD < 2.0" --filter-name "QD2" \
 -filter "SOR > 4.0" --filter-name "SOR4" \
 -filter "FS > 60.0" --filter-name "FS60" \
 -filter "MQ < 20.0" --filter-name "MQ20" \
 -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
 -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" --QUIET

#INFO field: Quality by Depth: "Variant Confidence/Quality by Depth"
#INFO field: Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias"
#INFO field: Fisher Strand: "Phred-scaled p-value using Fisher's exact test to detect strand bias"
#Careful, my distribution seems strange!!! INFO field: Description="RMS Mapping Quality"
#INFO field: MappingQualityRankSumTest: "Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities"
#INFO field: Description="Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias"

## 3a) Filter genotypes from INDELs and MIXED variants
$gatk --java-options "-Xmx30g" VariantFiltration -V $INDEL_unfilt -O ${wDir}/Ta_v3_INDELs_tempfilt_GWAspls.vcf.gz \
 -filter "QD < 2.0" --filter-name "QD2" \
 -filter "QUAL < 30.0" --filter-name "QUAL30" \
 -filter "FS > 200.0" --filter-name "FS200" \
 -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20"	--QUIET

## 4) Remove filtered sites from the files
$gatk --java-options "-Xmx30g" SelectVariants -V ${wDir}/Ta_v3_SNPs_tempfilt_GWAspls.vcf.gz --exclude-filtered -O $SNP_VarFilt --QUIET
$gatk --java-options "-Xmx30g" SelectVariants -V ${wDir}/Ta_v3_INDELs_tempfilt_GWAspls.vcf.gz --exclude-filtered -O $INDEL_VarFilt --QUIET

#rm ${wDir}/Ta_v3_*s_tempfilt_GWAspls.vcf.gz

## 5) Re-join SNP and VCF files
####$java -Xmx30g -jar ~/miniconda3/envs/bwa/share/picard-2.23.6-0/picard.jar GatherVcfs -I $SNP_filt -I $INDEL_filt -O ${outDir}/Ta_v3_vrts_VarFilt_GWAspls.vcf.gz
$java -Xmx30g -jar ~/miniconda3/envs/bwa/share/picard-2.23.6-0/picard.jar SortVcf -I $SNP_VarFilt -I $INDEL_VarFilt -O ${outDir}/Ta_v3_vrts_VarFilt_GWAspls.vcf.gz --TMP_DIR ${tmp}

## 6a) Filter for NAs<0.1 ONLY
$vcftools --gzvcf ${outDir}/Ta_v3_vrts_VarFilt_GWAspls.vcf.gz --max-missing 0.9 --recode --recode-INFO-all --stdout | $bgzip -c > $NAfilt_vcf
## 6b) Filter for NAs<0.1, biallelic sites and MAF
$vcftools --gzvcf ${outDir}/Ta_v3_vrts_VarFilt_GWAspls.vcf.gz --max-alleles 2 --max-missing 0.9 --maf 0.01 --recode --recode-INFO-all --stdout | $bgzip -c > $filtvcf_1

## 7) SCAFFOLD FILTERING: Scaffolds with too few positions in the vcf file should be removed or beagle will get stuck.
## zcat $filtvcf_1 | grep -v \# | cut -f1 | uniq -c | awk '$1<100'	# Check if necessary
bash ${lonely_vcf_pos} $filtvcf_1 $filtvcf_2

## 8) Fix vcftools converting phased missing positions ".|." into "." (https://github.com/vcftools/vcftools/issues/131)
## zcat $NAfilt_vcf | grep -o $'\t\.:' | wc -l		#Check whether this is necessary
zcat $NAfilt_vcf | sed $'s/\t\.:/\t\.|\.:/g' | $bgzip -c > ${outDir}/Ta_v3_vrts_NAfilt_tmp_GWAspls.vcf.gz
mv ${outDir}/Ta_v3_vrts_NAfilt_tmp_GWAspls.vcf.gz $NAfilt_vcf
zcat $filtvcf_2 | sed $'s/\t\.:/\t\.|\.:/g' | $bgzip -c > $final_vcf
rm $filtvcf_1 $filtvcf_2


echo ""
echo $(zcat $SNP_unfilt | grep -v \# | wc -l) original SNPs
echo $(zcat $SNP_VarFilt | grep -v \# | wc -l) SNPs after VariantFiltration
echo $(zcat $INDEL_unfilt | grep -v \# | wc -l) original INDELs
echo $(zcat $INDEL_VarFilt | grep -v \# | wc -l) INDELs after VariantFiltration
echo ""
echo $(zcat ${outDir}/Ta_v3_vrts_NAfilt_GWAspls.vcf.gz | grep -v \# | wc -l) variants \(SNPs + INDELs\) after filtering for --max-missing 0.9
echo $(zcat $final_vcf | grep -v \# | wc -l) variants \(SNPs + INDELs\) after filtering for --max-missing 0.9 --maf 0.05 --max-maf 0.95 --max-alleles 2 and scaffolds with <3 variants
echo $(zcat $final_vcf | grep -v \# | awk '{if (length($4)==1 && length($5)==1) {print}}' | wc -l) SNPs after filtering for --max-missing 0.9 --maf 0.05 --max-maf 0.95 --max-alleles 2 and scaffolds with \<3 variants
echo $(zcat $final_vcf | grep -v \# | awk '{ if (length($4)!=1 || length($5)!=1) {print}}' | wc -l) INDELs after filtering for --max-missing 0.9 --maf 0.05 --max-maf 0.95 --max-alleles 2 and scaffolds with \<3 variants
echo ""
echo Variants per chromosome
zcat $final_vcf | grep -v \# | cut -f1 | uniq -c

# Count SNPs in VCF file
#echo $(zcat $final_vcf | grep -v \# | awk '{if (length($4)==1 && length($5)==1 && $5!=".") {c++}} END {print c}')
#echo $(zcat $final_vcf | grep -v \# | awk '{if ($4~/^[ACGT]$/ && $5~/^[ACGT]$/){c++}} END {print c}')
