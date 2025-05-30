#!/bin/bash

bcftools norm -f /path/to/hg38.fa \
  -m -any \
  -Oz \
  -o Chinese.norm.vcf.gz \
  Chinese.vcf.gz

plink --vcf Chinese.norm.vcf.gz --split-x hg38 --chr chrX,chrY --threads 80 --keep-allele-order --make-bed --allow-extra-chr 0 --vcf-half-call m --out Chinese.hg38.chrXY
plink --bfile Chinese.hg38.chrXY --impute-sex 0.2 0.8 --keep-allele-order --make-bed --out Chinese.hg38.chrXY.sex
awk 'BEGIN {print "#IID SEX"} {print $2, $5}' Chinese.hg38.chrXY.sex.fam > Chinese.hg38.sex.info.txt
plink2 --vcf Chinese.norm.vcf.gz --split-par hg38 --update-sex Chinese.hg38.sex.info.txt --make-bed --allow-extra-chr 0 \
--vcf-half-call m --out Chinese.hg38

plink2 --bfile ASA.HK1 --exclude ASAHK1.duplicated.txt --recode vcf --output-chr chrMT --out JN
plink2 --bfile ASA.HK2 --recode vcf --output-chr chrMT --out HK2

bcftools +fixref HK1.vcf -Oz -o HK1.fixed.vcf.gz \
-- -f /path/to/hg38.fasta -m top

bcftools +fixref HK2.vcf -Oz -o HK2.fixed.vcf.gz \
-- -f /path/to/hg38.fasta -m top

plink --vcf HK2.fixed.vcf.gz --make-bed --keep-allele-order --double-id --out HK2.fixed
plink --vcf HK1.fixed.vcf.gz --make-bed --keep-allele-order --double-id --out HK1.fixed

awk '{$2 = "chr"$1":"$4":"$6":"$5; print}' Chinese.hg38.bim > temp.txt && rm Chinese.hg38.bim && mv temp.txt Chinese.hg38.bim
awk '{$2 = $1":"$4":"$6":"$5; print}' Thai_chip_data.chr6.fixed.bim > temp.txt && rm Thai_chip_data.chr6.fixed.bim && mv temp.txt Thai_chip_data.chr6.fixed.bim

awk '{print $2}' Chinese.hg38.bim | sort > Chinese.hg38_id.txt
awk '{print $2}' HK1.fixed.bim | sort > HK1.fixed_id.txt
awk '{print $2}' HK2.fixed.bim | sort > HK2.fixed_id.txt

comm -12 HK1.fixed_id.txt HK2.fixed_id.txt > intersection.txt
comm -12 intersection.txt Chinese.hg38_id.txt > intersection2.txt

plink --bfile HK2.fixed --extract intersection2.txt --keep-allele-order --recode vcf --output-chr chrMT --out HK2.samesnps
plink --bfile HK1.fixed --extract intersection2.txt --keep-allele-order --recode vcf --output-chr chrMT --out HK1.samesnps
plink2 --bfile Chinese.hg38 --extract intersection2.txt --recode vcf --output-chr chrMT --out Chinese.samesnps

bgzip -@ 30 HK2.samesnps.vcf && tabix HK2.samesnps.vcf.gz
bgzip -@ 30 HK1.samesnps.vcf && tabix HK1.samesnps.vcf.gz
bgzip -@ 30 Chinese.samesnps.vcf && tabix Chinese.samesnps.vcf.gz

bcftools merge HK1.samesnps.vcf.gz HK2.samesnps.vcf.gz Chinese.samesnps.vcf.gz -Oz -o ASA.Chinese.vcf.gz

plink --vcf ASA.Chinese.vcf.gz --keep-allele-order --make-bed --out ASA.Chinese --double-id

awk '{$1=0; print $0}' ASA.Chinese.fam > temp.txt && mv temp.txt ASA.Chinese.fam
plink --bfile ASA.Chinese --het --out WGSASA_het
sed 's/ \+/ /g' WGSASA_het.het > WGSASA.het
plink --bfile ASA.Chinese --remove het.txt --keep-allele-order --make-bed --out ASA.Chinese.het

plink --bfile ASA.Chinese.het --genome --min 0.8 --out filtered_pihat_min0.8
sed 's/ \+/ /g' filtered_pihat_min0.8.genome | awk '{print $3, $4}' > filtered_pihat_min0.8.txt
plink --bfile ASA.Chinese.het --remove filtered_pihat_min0.8.txt --keep-allele-order --make-bed --out ASA.Chinese.het.genome

plink --bfile ASA.Chinese.het.genome --indep-pairwise 50 5 0.2 --out indepSNPThai.WGSASA
plink --bfile ASA.Chinese.het.genome --extract indepSNPThai.WGSASA.prune.in --pca --out pca
plink --bfile ASA.Chinese.het.genome --remove pcasamples.txt --keep-allele-order --make-bed --out ASA.Chinese.het.genome.pca

plink --bfile ASA.Chinese.het.genome.pca --indep-pairwise 50 5 0.2 --out indepSNPThai.WGSASA
plink --bfile ASA.Chinese.het.genome.pca --extract indepSNPThai.WGSASA.prune.in --pca --out pca.after

plink --bfile ASA.Chinese.het.genome.pca --genome --min 0.185 --out filtered_pihat_min0.185
awk '{$1=$1}1' filtered_pihat_min0.185.genome > temp.txt && mv temp.txt filtered_pihat_min0.185.genome
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' filtered_pihat_min0.185.genome > pihat.txt

plink --bfile ASA.Chinese.het.genome.pca --geno 0.05 --maf 0.01 --hwe 0.00001 --keep-allele-order --make-bed --out ASA.Chinese.het.genome.pca.QC

plink --bfile ASA.Chinese.het.genome.pca.QC --recode vcf --keep-allele-order --output-chr chrMT --out ASA.Chinese.QC
bcftools sort ASA.Chinese.QC.vcf -Ov -o ASA.Chinese.QC.Sorted.vcf 
bgzip -@ 30 ASA.Chinese.QC.Sorted.vcf && tabix -p vcf ASA.Chinese.QC.Sorted.vcf.gz



for chr in chr{1..22}
do
  bcftools view -r $chr --output-type z --output ASA.Chinese.before_phasing.${chr}.vcf.gz ASA.Chinese.QC.Sorted.vcf.gz
done


for chr in {1..22}
do
/path/to/eagle --vcf=ASA.Chinese.before_phasing.${chr}.vcf.gz \
--geneticMapFile=/path/to/genetic_map_hg38_withX.txt.gz --numThreads=5 --outPrefix ASA.Chinese.chr${chr}.phased \
--vcfOutFormat z &
done

