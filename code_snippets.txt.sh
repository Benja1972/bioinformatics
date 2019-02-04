$ bgzip -c WES_TLX3_TAP_Tum7-8.vcf > WES_TLX3_TAP_Tum7-8.vcf.gz

$ tabix -p vcf WES_TLX3_TAP_Tum7-8.vcf.gz 

$ bcftools stats WES_TLX3_TAP_Tum7-8_TAP-01.vcf.gz WES_TLX3_TAP_Tum7-8_TLX3-01.vcf.gz > joined_stats.txt


$ plot-vcfstats joined_stats.txt -p WES_joined 


####### Split VCF file by samples

#!/bin/bash

[ ! -f "$1" ] && printf "\n\n>> ERROR : File Not Found : $1 \n" && exit 1; 
FILE="$1";

for sample in `bcftools query -l $FILE`; do
  bcftools view -c1 -Oz -s $sample -o ${FILE/.vcf*/.$sample.vcf.gz} $FILE
done

#######################


# Split VCF file by samples
$ bcftools view -c1 -Oz -s AC3812,AC3813 -o TAP_WGS.vcf.gz FERRIER_09_Germline.allchr.snpEff.p.SAL.SAL10_1.vcf.gz

# Intersect VCF
bcftools isec -p dir A.vcf.gz B.vcf.gz
bcftools isec -p dir A.vcf.gz B.vcf.gz C.vcf.gz D.vcf.gz E.vcf.gz F.vcf.gz -n=6


########################################
# Call variants by bcftools
# http://samtools.github.io/bcftools/howtos/variant-calling.html

bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf


bcftools view calls.bcf | vcfutils.pl varFilter - > calls_filter.vcf