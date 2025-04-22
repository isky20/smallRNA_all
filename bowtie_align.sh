#!/bin/bash

# -----------------------------
# Setup
# -----------------------------
GENOME="genome_sus_scrofa.fa"   # Reference genome (Sus scrofa)
GTF="sus_scrofa.gtf"            # Annotation file in GTF format
THREADS=8                       # Number of threads to use
INDEX_PREFIX="genome_index"     # Prefix for Bowtie1 index

# -----------------------------
# Index the genome (if not already)
# -----------------------------
echo "Checking for genome index..."
if [ ! -f "${INDEX_PREFIX}.1.ebwt" ]; then
    echo "Indexing genome with bowtie-build..."
    bowtie-build $GENOME $INDEX_PREFIX
fi

# -----------------------------
# Process each FASTQ file
# -----------------------------
for FASTQ in *.fastq; do
    BASENAME=$(basename "$FASTQ" .fastq)
    
    echo "Processing $FASTQ..."

    # Run Bowtie1 alignment
    bowtie -v 1 -a --best --strata -p $THREADS -q $INDEX_PREFIX "$FASTQ" -S "${BASENAME}.sam"

    # Convert SAM to BAM and sort
    samtools view -bS "${BASENAME}.sam" > "${BASENAME}.bam"
    samtools sort "${BASENAME}.bam" -o "${BASENAME}_sorted.bam"
    samtools index "${BASENAME}_sorted.bam"

    # Cleanup intermediate files
    rm "${BASENAME}.sam" "${BASENAME}.bam"

    # Run featureCounts for gene quantification
    featureCounts -T $THREADS -a $GTF -o "${BASENAME}_counts.txt" -t exon -g gene_id "${BASENAME}_sorted.bam"

    echo "Finished processing $FASTQ. Count matrix: ${BASENAME}_counts.txt"
done

echo "✅ All FASTQ files processed and quantified."
