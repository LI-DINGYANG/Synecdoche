#!/bin/bash
#BSUB -J GPU_A40  # Job name
#BSUB -o GPU_A40.out  # Standard output file
#BSUB -e GPU_A40.error  # Standard error file
#BSUB -q qa40  # Queue/partition name
#BSUB -m "gpu05"  # Specific host to run on
#BSUB -n 64  # Number of CPU cores
#BSUB -gpu "num=8:aff=yes"  # GPU resources (8 GPUs with affinity)

# Load required modules and define paths (USER MUST UPDATE THESE)
SINGULARITY_IMAGE=/path/to/clara_parabricks-4.4.0-1.sif  # Path to Singularity image
singularity=~/path/to/singularity  # Path to Singularity executable
tmp_dir=/path/to/tmpdir  # Temporary directory
bcftools=~/path/to/bcftools  # Path to bcftools
bgzip=~/path/to/bgzip  # Path to bgzip
tabix=~/path/to/tabix  # Path to tabix
samtools=~/path/to/samtools  # Path to samtools
java=~/path/to/java  # Path to Java executable
refpath=/path/to/reference_files  # Directory containing reference files
INPUT_DIR=/path/to/input  # Input directory
OUTPUT_DIR=/path/to/output  # Output directory
output_samples_file=/path/to/sample_list.txt  # File containing sample information

# Function to process genomic data
processGenomicsData() {
  local SAMPLE_NAME=$1
  local FQ1=$2
  local FQ2=$3
  local rg_tag=$4
  local Corhot=$5

  # Run Parabricks germline pipeline using Singularity
  $singularity exec --nv \
    --bind ${INPUT_DIR}:/workdir \
    --bind ${OUTPUT_DIR}:/outputdir \
    --bind ${refpath}:/refdir \
    --bind ${tmp_dir}:/tmpdir \
    $SINGULARITY_IMAGE \
    pbrun germline \
    --ref /refdir/hg38.fa \
    --in-fq /workdir/${FQ1} /workdir/${FQ2} $rg_tag \
    --knownSites /refdir/Homo_sapiens_assembly38.known_indels.vcf.gz \
    --knownSites /refdir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --knownSites /refdir/Homo_sapiens_assembly38.dbsnp138.vcf.gz \
    --out-recal-file /outputdir/${Corhot}/$SAMPLE_NAME/bwa/${SAMPLE_NAME}.GPU.BQSR_report.txt \
    --out-bam /outputdir/${Corhot}/$SAMPLE_NAME/bwa/${SAMPLE_NAME}.GPU.BQSR.bam \
    --tmp-dir /tmpdir \
    --memory-limit 500 \
    --num-cpu-threads-per-stage 16 \
    --bwa-cpu-thread-pool 16 \
    --gvcf \
    --out-variants /outputdir/${Corhot}/${SAMPLE_NAME}/gatk/${SAMPLE_NAME}.GPU.g.vcf \
    --run-partition \
    --bwa-options="-M -Y -K 10000000" \
    --haplotypecaller-options '-A QualByDepth -A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A InbreedingCoeff' \
    --gpusort \
    --gpuwrite \
    --num-htvc-threads 8 \
    --gpu-num-per-partition 2 \
    --num-gpus 8 && echo "** ${SAMPLE_NAME} GPU germline done **"  # Total GPUs used

  # Compress the GVCF file
  $bgzip -f -@ 50 ${OUTPUT_DIR}/${Corhot}/${SAMPLE_NAME}/gatk/${SAMPLE_NAME}.GPU.g.vcf
}

# Process each sample from the sample list file
while read sample fq1 fq2 library; do
    outsampledir=${OUTPUT_DIR}/${library}/${sample}

    # Create output directories if they don't exist
    if [ ! -d $outsampledir/bwa ]; then
        mkdir -p $outsampledir/bwa
    fi

    if [ ! -d $outsampledir/gatk ]; then
        mkdir -p $outsampledir/gatk
    fi
    
    echo "Processing $sample"
    RGTAG="@RG\tID:$sample\tLB:lib\tPL:Illumina\tSM:$sample\tPU:$sample"  # Read group information
    processGenomicsData $sample $fq1 $fq2 $RGTAG $library && echo "** ${sample} GPU GATK GVCF done **"
done < $output_samples_file
wait

# Create population directory if it doesn't exist
if [ ! -d ${OUTPUT_DIR}/population ]; then
    mkdir -p ${OUTPUT_DIR}/population
fi




