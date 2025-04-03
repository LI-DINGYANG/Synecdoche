#!/bin/bash
#BSUB -J gl  # Job name
#BSUB -o gl.out  # Standard output file
#BSUB -e gl.error  # Standard error file
#BSUB -q q8358  # Queue/partition name
#BSUB -n 320  # Number of CPU cores


# Define tool paths (USER MUST UPDATE THESE BASED ON THEIR ENVIRONMENT)
singularity=~/path/to/singularity  # Path to Singularity executable
bcftools=~/path/to/bcftools  # Path to bcftools
bgzip=~/path/to/bgzip  # Path to bgzip
tabix=~/path/to/tabix  # Path to tabix
java=~/path/to/java  # Path to Java executable

# Define input/output directories (USER MUST UPDATE THESE)
INPUT_DIR=/path/to/input  # Directory containing input files
OUTPUT_DIR=/path/to/output  # Directory for results

# Create output directories if they don't exist
mkdir -p "${OUTPUT_DIR}/population"
mkdir -p "${OUTPUT_DIR}/population/log"

# Cohort name (modify as needed)
corhotname="WGS"  

# Run GLnexus for joint genotyping
$singularity exec \
     --bind "${OUTPUT_DIR}:/outputdir" \
    /path/to/glnexus.sif \
    /path/to/glnexus_cli \
    --dir /outputdir/population/temp \
    --config gatk \
    --list /outputdir/gvcf.txt \
    --trim-uncalled-alleles \
    --more-PL > "${OUTPUT_DIR}/population/GLnexus.${corhotname}.bcf" && echo "** GLnexus done **"  

# Convert BCF to VCF
$bcftools view "${OUTPUT_DIR}/population/GLnexus.${corhotname}.bcf" -o "${OUTPUT_DIR}/population/GLnexus.${corhotname}.vcf"  

# Compress and index the VCF
$bgzip -@ 100 "${OUTPUT_DIR}/population/GLnexus.${corhotname}.vcf"  # Multithreaded compression
$tabix -p vcf "${OUTPUT_DIR}/population/GLnexus.${corhotname}.vcf.gz"  # Create tabix index



