#!/bin/bash
#BSUB -J vcftools               # Job name
#BSUB -o vcftools.%J.out        # Standard output file
#BSUB -e vcftools.%J.err        # Standard error file
#BSUB -q q6240                  # Queue name
#BSUB -n 30                     # Number of CPU cores

# Normalize VCF file (left-align indels, split multiallelic sites)
bcftools norm -f /path/to/hg38.fa \
  -m -any \
  -Oz \
  -o Chinese.WGS1263.norm.vcf.gz \
  Chinese.WGS1263.vcf.gz

# Process each chromosome in parallel
for i in chr{1..22}
do
  # Extract chromosome-specific variants
  bcftools view -r $i Chinese.WGS1263.norm.vcf.gz > Chinese.WGS1263.${i}.vcf && 
  
  # Apply quality filters
  vcftools \
    --vcf Chinese.WGS1263.${i}.vcf \
    --min-meanDP 5 \
    --minGQ 20 \
    --max-missing 0.8 \
    --out Chinese.WGS1263.filter.${i} \
    --recode --recode-INFO-all && # Keep all INFO fields
  
  # Convert to PLINK binary format
  plink --vcf /Chinese.WGS1263.filter.${i}.recode.vcf \
    --vcf-half-call m \
    --make-bed \
    --out Chinese.WGS1263.filter.${i} && 
  
  # Calculate allele frequencies
  plink --bfile Chinese.WGS1263.filter.${i} \
    --freq \
    --out Chinese.WGS1263.filter.${i} &
done && wait  # Wait for all chromosomes to finish processing

# Extract genotype information in custom format
for i in chr{1..22}
do
  # Extract CHROM, POS, ID, REF, ALT and all genotype fields
  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' \
    Chinese.WGS1263.filter.${i}.recode.vcf > Chinese.WGS1263.filter.${i}.modify.txt &
done && wait

# Convert genotype coding to numerical format
for i in chr{1..22}
do
  gawk '{
    # Create new variant ID: CHROM:POS:REF:ALT
    $3 = $1":"$2":"$4":"$5; 
    print
  }' Chinese.WGS1263.filter.${i}.modify.txt | 
  gawk '{
    # Convert genotype coding:
    # 0/0 or 0|0 -> 0 (homozygous reference)
    # 1/1 or 1|1 -> 2 (homozygous alternate)
    # 0/1, 1/0, etc -> 1 (heterozygous)
    # Others -> -1 (missing)
    for(i=6; i<=NF; i++) { 
      if($i=="0/0" || $i=="0|0") $i=0; 
      else if($i=="1/1" || $i=="1|1") $i=2; 
      else if($i=="1/0" || $i=="0/1" || $i=="1|0" || $i=="0|1") $i=1; 
      else $i=-1; 
    } 
    # Reformat output keeping first 5 columns and all genotype columns
    output = $1; 
    for(j=2; j<=5; j++) output = output "\t" $j; 
    for(k=6; k<=NF; k++) output = output "\t" $k; 
    print output
  }' > Chinese.WGS1263.${i}.modify.txt &
done && wait
