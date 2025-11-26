# Snakemake workflow for reference preparation and analysis
# Run with: snakemake --cores 16 --configfile config.yaml

import os

# Load configuration
configfile: "config.yaml"

# Configuration variables
PICARDJAR = config["picard_jar"]
SAMPLELIST = config["sample_list"]
FASTQDIR = config["fastq_dir"]
FQENDING = config["fastq_ending"]
REF = config["reference"]["fasta"]
REFPREFIX = config["reference"]["prefix"]
REFDIR = config["reference"]["output_dir"]

# Reference files
MASKEDREF = f"{REFDIR}/{REFPREFIX}.fa"

# Read sample list
with open(SAMPLELIST, 'r') as f:
    SAMPLES = [line.strip() for line in f if line.strip()]

# Rule to run everything (reference prep + alignment + picard processing + samtools stats)
rule all:
    input:
        f"{MASKEDREF}.bwt.2bit.64",
        f"{MASKEDREF}.fai",
        f"{REFDIR}/{REFPREFIX}.dict",
        expand("bam-alignments/{sample}-RG.bam", sample=SAMPLES),
        expand("bam-alignments/{sample}-RG.bam.bai", sample=SAMPLES),
        expand("samtools-stats/{sample}-coverage.txt", sample=SAMPLES),
        expand("samtools-stats/{sample}-stats.txt", sample=SAMPLES),
        expand("samtools-stats/{sample}-idxstats.txt", sample=SAMPLES),
        expand("samtools-stats/{sample}-flagstat.txt", sample=SAMPLES)

# Rule to run just the reference preparation
rule prepare_reference:
    input:
        f"{MASKEDREF}.bwt.2bit.64",
        f"{MASKEDREF}.fai",
        f"{REFDIR}/{REFPREFIX}.dict"

# Step 1: Run nucmer to identify repeats
rule nucmer_self_alignment:
    input:
        ref = REF
    output:
        delta = f"{REFDIR}/{REFPREFIX}.delta"
    params:
        prefix = f"{REFDIR}/{REFPREFIX}"
    threads: config["resources"]["threads"]["nucmer"]
    resources:
        mem_mb = config["resources"]["memory"]["nucmer"]
    log:
        f"{REFDIR}/logs/nucmer.log"
    shell:
        """
        mkdir -p {REFDIR}/logs
        nucmer --maxmatch --nosimplify -p {params.prefix} -t {threads} {input.ref} {input.ref} > {log} 2>&1
        """

# Step 2: Generate repeat coordinates
rule show_coords:
    input:
        delta = f"{REFDIR}/{REFPREFIX}.delta"
    output:
        coords = f"{REFDIR}/masked_ref_BEFORE_ORDER.bed"
    threads: config["resources"]["threads"]["default"]
    resources:
        mem_mb = config["resources"]["memory"]["default"]
    log:
        f"{REFDIR}/logs/show_coords.log"
    shell:
        """
        show-coords -r -T -H {input.delta} > {output.coords} 2> {log}
        """

# Step 3: Filter and format repeat coordinates
rule format_repeat_bed:
    input:
        coords = f"{REFDIR}/masked_ref_BEFORE_ORDER.bed"
    output:
        bed = f"{REFDIR}/masked_ref.bed"
    threads: config["resources"]["threads"]["default"]
    resources:
        mem_mb = config["resources"]["memory"]["default"]
    log:
        f"{REFDIR}/logs/format_bed.log"
    shell:
        """
        awk '{{if ($1 != $3 && $2 != $4) print $0}}' {input.coords} | \
        awk '{{print $8 "\\t" $1 "\\t" $2}}' > {output.bed} 2> {log}
        """

# Step 4: Mask reference with bedtools
rule mask_reference:
    input:
        ref = REF,
        bed = f"{REFDIR}/masked_ref.bed"
    output:
        masked = MASKEDREF
    threads: config["resources"]["threads"]["default"]
    resources:
        mem_mb = config["resources"]["memory"]["maskfasta"]
    log:
        f"{REFDIR}/logs/maskfasta.log"
    shell:
        """
        bedtools maskfasta -fi {input.ref} -bed {input.bed} -fo {output.masked} > {log} 2>&1
        """

# Step 5: Index with BWA-MEM2
rule bwa_index:
    input:
        ref = MASKEDREF
    output:
        multiext(MASKEDREF, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    threads: config["resources"]["threads"]["bwa_index"]
    resources:
        mem_mb = config["resources"]["memory"]["bwa_index"]
    log:
        f"{REFDIR}/logs/bwa_index.log"
    shell:
        """
        bwa-mem2 index {input.ref} > {log} 2>&1
        """

# Step 6: Index with samtools
rule samtools_faidx:
    input:
        ref = MASKEDREF
    output:
        fai = f"{MASKEDREF}.fai"
    threads: config["resources"]["threads"]["default"]
    resources:
        mem_mb = config["resources"]["memory"]["samtools"]
    log:
        f"{REFDIR}/logs/samtools_faidx.log"
    shell:
        """
        samtools faidx {input.ref} > {log} 2>&1
        """

# Step 7: Create Picard sequence dictionary
rule picard_dict:
    input:
        ref = MASKEDREF
    output:
        dict = f"{REFDIR}/{REFPREFIX}.dict"
    threads: config["resources"]["threads"]["default"]
    resources:
        mem_mb = config["resources"]["memory"]["picard"]
    log:
        f"{REFDIR}/logs/picard_dict.log"
    shell:
        """
        java -Xmx16g -jar {PICARDJAR} CreateSequenceDictionary \
            -R {input.ref} \
            -O {output.dict} > {log} 2>&1
        """

# ==============================================================================
# READ ALIGNMENT SECTION
# ==============================================================================

# Rule to run alignment for all samples
rule align_all:
    input:
        expand("bam-alignments/{sample}.bam", sample=SAMPLES)

# Step 8: Align reads with BWA-MEM2
rule bwa_mem2_align:
    input:
        ref = MASKEDREF,
        idx = f"{MASKEDREF}.bwt.2bit.64",
        fq1 = f"{FASTQDIR}/{{sample}}{FQENDING}1.fastq.gz",
        fq2 = f"{FASTQDIR}/{{sample}}{FQENDING}2.fastq.gz"
    output:
        sam = temp("bam-alignments/{sample}.sam")
    params:
        rg = lambda wildcards: f"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA"
    threads: config["resources"]["threads"]["bwa_align"]
    resources:
        mem_mb = config["resources"]["memory"]["bwa_align"]
    log:
        "bam-alignments/logs/{sample}_bwa.log"
    shell:
        """
        mkdir -p bam-alignments/logs
        bwa-mem2 mem -R '{params.rg}' -t {threads} {input.ref} {input.fq1} {input.fq2} \
            > {output.sam} 2> {log}
        """

# Step 9: Sort SAM to BAM with samtools
rule samtools_sort:
    input:
        sam = "bam-alignments/{sample}.sam"
    output:
        bam = "bam-alignments/{sample}.bam"
    threads: config["resources"]["threads"]["samtools_sort"]
    resources:
        mem_mb = config["resources"]["memory"]["samtools_sort"]
    log:
        "bam-alignments/logs/{sample}_sort.log"
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.sam} > {log} 2>&1
        """

# ==============================================================================
# PICARD PROCESSING SECTION
# ==============================================================================

# Rule to run all Picard processing steps
rule process_bams:
    input:
        expand("bam-alignments/{sample}-RG.bam", sample=SAMPLES),
        expand("bam-alignments/{sample}-RG.bai", sample=SAMPLES)

# Step 10: Mark duplicates with Picard
rule mark_duplicates:
    input:
        bam = "bam-alignments/{sample}.bam"
    output:
        bam = temp("bam-alignments/{sample}-markdup.bam"),
        metrics = "bam-alignments/{sample}-markdup-metrics.txt"
    threads: config["resources"]["threads"]["default"]
    resources:
        mem_mb = config["resources"]["memory"]["picard_markdup"]
    log:
        "bam-alignments/logs/{sample}_markdup.log"
    shell:
        """
        java -Xmx4g -jar {PICARDJAR} MarkDuplicates \
            -I {input.bam} \
            -O {output.bam} \
            -M {output.metrics} \
            --REMOVE_DUPLICATES true \
            --ASSUME_SORT_ORDER coordinate \
            --VALIDATION_STRINGENCY LENIENT > {log} 2>&1
        """

# Step 11: Clean SAM with Picard
rule clean_sam:
    input:
        bam = "bam-alignments/{sample}-markdup.bam"
    output:
        bam = temp("bam-alignments/{sample}-cleansam.bam")
    threads: config["resources"]["threads"]["default"]
    resources:
        mem_mb = config["resources"]["memory"]["picard_cleansam"]
    log:
        "bam-alignments/logs/{sample}_cleansam.log"
    shell:
        """
        java -Xmx4g -jar {PICARDJAR} CleanSam \
            -I {input.bam} \
            -O {output.bam} > {log} 2>&1
        """

# Step 12: Fix mate information with Picard
rule fix_mate_information:
    input:
        bam = "bam-alignments/{sample}-cleansam.bam"
    output:
        bam = temp("bam-alignments/{sample}-fixmate.bam")
    threads: config["resources"]["threads"]["default"]
    resources:
        mem_mb = config["resources"]["memory"]["picard_fixmate"]
    log:
        "bam-alignments/logs/{sample}_fixmate.log"
    shell:
        """
        java -Xmx4g -jar {PICARDJAR} FixMateInformation \
            -I {input.bam} \
            -O {output.bam} > {log} 2>&1
        """

# Step 13: Add or replace read groups with Picard (final BAM)
rule add_replace_read_groups:
    input:
        bam = "bam-alignments/{sample}-fixmate.bam"
    output:
        bam = "bam-alignments/{sample}-RG.bam",
        bai = "bam-alignments/{sample}-RG.bam.bai"
    params:
        sample = "{sample}"
    threads: config["resources"]["threads"]["picard_rg"]
    resources:
        mem_mb = config["resources"]["memory"]["picard_rg"]
    log:
        "bam-alignments/logs/{sample}_RG.log"
    shell:
        """
        java -Xmx4g -XX:ParallelGCThreads={threads} -jar {PICARDJAR} AddOrReplaceReadGroups \
            --INPUT {input.bam} \
            --OUTPUT {output.bam} \
            -ID {params.sample} \
            -LB {params.sample} \
            -PL Illumina \
            -SM {params.sample} \
            -CREATE_INDEX true \
            -PU {params.sample} > {log} 2>&1
        mv bam-alignments/{params.sample}-RG.bai bam-alignments/{params.sample}-RG.bam.bai
        """

# ==============================================================================
# SAMTOOLS STATISTICS SECTION
# ==============================================================================

# Rule to run all samtools statistics
rule samtools_statistics:
    input:
        expand("samtools-stats/{sample}-coverage.txt", sample=SAMPLES),
        expand("samtools-stats/{sample}-stats.txt", sample=SAMPLES),
        expand("samtools-stats/{sample}-idxstats.txt", sample=SAMPLES),
        expand("samtools-stats/{sample}-flagstat.txt", sample=SAMPLES)

# Step 14: Generate coverage statistics
rule samtools_coverage:
    input:
        bam = "bam-alignments/{sample}-RG.bam",
        bai = "bam-alignments/{sample}-RG.bam.bai"
    output:
        coverage = "samtools-stats/{sample}-coverage.txt"
    threads: config["resources"]["threads"]["default"]
    resources:
        mem_mb = config["resources"]["memory"]["samtools"]
    log:
        "samtools-stats/logs/{sample}_coverage.log"
    shell:
        """
        mkdir -p samtools-stats/logs
        samtools coverage -o {output.coverage} {input.bam} > {log} 2>&1
        """

# Step 15: Generate alignment statistics
rule samtools_stats:
    input:
        bam = "bam-alignments/{sample}-RG.bam",
        bai = "bam-alignments/{sample}-RG.bam.bai",
        ref = MASKEDREF
    output:
        stats = "samtools-stats/{sample}-stats.txt"
    threads: config["resources"]["threads"]["samtools_sort"]
    resources:
        mem_mb = config["resources"]["memory"]["samtools"]
    log:
        "samtools-stats/logs/{sample}_stats.log"
    shell:
        """
        samtools stats -@ {threads} --reference {input.ref} {input.bam} > {output.stats} 2> {log}
        """

# Step 16: Generate index statistics
rule samtools_idxstats:
    input:
        bam = "bam-alignments/{sample}-RG.bam",
        bai = "bam-alignments/{sample}-RG.bam.bai"
    output:
        idxstats = "samtools-stats/{sample}-idxstats.txt"
    threads: config["resources"]["threads"]["samtools_sort"]
    resources:
        mem_mb = config["resources"]["memory"]["samtools"]
    log:
        "samtools-stats/logs/{sample}_idxstats.log"
    shell:
        """
        samtools idxstats -@ {threads} {input.bam} > {output.idxstats} 2> {log}
        """

# Step 17: Generate flag statistics
rule samtools_flagstat:
    input:
        bam = "bam-alignments/{sample}-RG.bam",
        bai = "bam-alignments/{sample}-RG.bam.bai"
    output:
        flagstat = "samtools-stats/{sample}-flagstat.txt"
    threads: config["resources"]["threads"]["samtools_sort"]
    resources:
        mem_mb = config["resources"]["memory"]["samtools"]
    log:
        "samtools-stats/logs/{sample}_flagstat.log"
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat} 2> {log}
        """

# ==============================================================================
# UTILITY RULES
# ==============================================================================

# Rule to clean up intermediate alignment BAM files
# Run with: snakemake clean_alignment_bams --cores 1
rule clean_alignment_bams:
    input:
        expand("bam-alignments/{sample}-RG.bam", sample=SAMPLES)
    output:
        touch("bam-alignments/.cleanup_complete")
    params:
        bams = expand("bam-alignments/{sample}.bam", sample=SAMPLES)
    shell:
        """
        echo "Removing intermediate alignment BAM files..."
        rm -f {params.bams}
        echo "Cleanup complete. Intermediate BAM files removed."
        """
