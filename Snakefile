import pandas as pd
from scripts.aux_func import parse_xml, get_targets_from_csv

# Configuration file specifying parameters and settings for the workflow
configfile: "config.yaml"

parse_xml(
    input="data/Sample_Information/Sample_Information.xml",
    output="data/Sample_Information/Extracted_Sample_Information.csv",
)

# Load the sample information CSV and create a dict to store the information
sample_info_df = pd.read_csv('data/Sample_Information/Extracted_Sample_Information.csv')
specs = list(sample_info_df['Species'].unique())
sample_map = sample_info_df.set_index('Accession')['Sample_ID'].to_dict()
species_accessions = sample_info_df.groupby("Species")["Accession"].apply(list).to_dict()

#accessions = list(sample_map.keys())


# Main rule specifying the final targets to be generated
rule all:
    input:
        reference='data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna',
        
        reference_index='data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna.fai',
        
        dict=f"data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.dict",
        
        qc=get_targets_from_csv(
        	'data/Sample_Information/Extracted_Sample_Information.csv',
        	'_1_fastqc.html'
        ),
        
        re_qc=get_targets_from_csv(
        	'data/Sample_Information/Extracted_Sample_Information.csv',
        	'_R1_trimmed_fastqc.html'	
        ),
        
        merged_rg_file=[
            f"data/Genomic_data/{Species}/merged_bam/{Species}_rg.txt" for Species in species_accessions.keys()
        ],
        
        final_vcf=[
            f"data/Genomic_data/{Species}/final_vcf/{Species}_gatk_filtered.vcf" for Species in species_accessions.keys()
        ]
        

# Rule to download and decompress a reference genome from NCBI
rule download_reference:
    output:
        reference_gz='data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna.gz',
        reference='data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna'
    params:
        link='https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/775/205/GCF_002775205.1_X_maculatus-5.0-male/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna.gz'
    conda:
        "environment.yml"
    shell:
        """
        # Download the genome file and decompress it
        wget -O {output.reference_gz} {params.link}
        gunzip -c {output.reference_gz} > {output.reference}
        """

# Rule to download paired-end FASTQ files based on sample information in the CSV file
rule download_fastq:
    input:
        csv="data/Sample_Information/Extracted_Sample_Information.csv"
    output:
        fastq_1="data/Genomic_data/{Species}/{Accession}/{Accession}_1.fastq",
        fastq_2="data/Genomic_data/{Species}/{Accession}/{Accession}_2.fastq"
    params:
        fastq_dump = config['fastq_dump_path']
    conda:
        "environment.yml"
    shell:
        """
        # Create directories and download FASTQ files using fastq-dump
        mkdir -p data/Genomic_data/{wildcards.Species}/{wildcards.Accession}
        {params.fastq_dump} --split-3 --outdir data/Genomic_data/{wildcards.Species}/{wildcards.Accession} {wildcards.Accession}
        """

# Rule to generate FastQC reports for quality control of the downloaded FASTQ files
rule check_quality_control:
    input:
        fastq_1="data/Genomic_data/{Species}/{Accession}/{Accession}_1.fastq",
        fastq_2="data/Genomic_data/{Species}/{Accession}/{Accession}_2.fastq"
    output:
        fastqc_1="data/Genomic_data/{Species}/{Accession}/{Accession}_1_fastqc.html",
        fastqc_2="data/Genomic_data/{Species}/{Accession}/{Accession}_2_fastqc.html"
    params:
        fastq_folder="data/Genomic_data/{Species}/{Accession}"
    shell:
        """
        # Run FastQC for the first FASTQ file and save the report in the specified folder
        fastqc {input.fastq_1} {input.fastq_2} -o {params.fastq_folder}
        """

# Rule to trim FASTQ files using Trimmomatic.
rule trimming_fastq:
    input:
        fastq_1="data/Genomic_data/{Species}/{Accession}/{Accession}_1.fastq",
        fastq_2="data/Genomic_data/{Species}/{Accession}/{Accession}_2.fastq"
    output:
        fastq_trim_1="data/Genomic_data/{Species}/{Accession}/{Accession}_R1_trimmed.fastq",
        fastq_trim_2="data/Genomic_data/{Species}/{Accession}/{Accession}_R2_trimmed.fastq",
        log="data/Genomic_data/{Species}/{Accession}/{Accession}_trimmomatic.log"
    params:
        trimmomatic_params="ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    conda:
        "environment.yml"
    shell:
        """
        # Execute Trimmomatic on the input FASTQ files with specified parameters.
        trimmomatic PE -phred33 {input.fastq_1} {input.fastq_2} \
        {output.fastq_trim_1} {output.fastq_trim_1}.unpaired \
        {output.fastq_trim_2} {output.fastq_trim_2}.unpaired \
        {params.trimmomatic_params} \
        > {output.log} 2>&1
        """

# Rule to generate FastQC reports from trimmed FASTQ files for quality reassessment
rule recheck_quality_control:
    input:
        fastqc_1="data/Genomic_data/{Species}/{Accession}/{Accession}_1_fastqc.html",
        fastq_trim_1="data/Genomic_data/{Species}/{Accession}/{Accession}_R1_trimmed.fastq",
        fastq_trim_2="data/Genomic_data/{Species}/{Accession}/{Accession}_R2_trimmed.fastq"
    output:
        fastqc_1="data/Genomic_data/{Species}/{Accession}/{Accession}_R1_trimmed_fastqc.html",
        fastqc_2="data/Genomic_data/{Species}/{Accession}/{Accession}_R2_trimmed_fastqc.html"
    params:
        fastq_folder="data/Genomic_data/{Species}/{Accession}"
    shell:
        """
        # Run FastQC for the first FASTQ file and save the report in the specified folder
        fastqc {input.fastq_trim_1} {input.fastq_trim_2} -o {params.fastq_folder}
        """

# Index reference genome and create auxiliary files for each weekly run
rule index_reference:
    input:
        # Reference genome
        reference='data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna',
    output:
        # Indexed reference files
        index=expand("data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna.{suffix}", suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
        fai=f"data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna.fai",
        dicti=f"data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.dict",
    conda:
        # Specify the Conda environment.
        "environment.yml"
    shell:
        """
        # Index the reference genome
        bwa index {input.reference}
        
        # Generate the FASTA index (.fai) file using samtools faidx
        samtools faidx {input.reference}
        
        # Generate the sequence dictionary (.dict) file using Picard
        picard CreateSequenceDictionary R={input.reference} O={output.dicti}
        """


# Align trimmed FASTQ files using BWA for each weekly run
rule align_fastqc:
    input:
        # Reference genome file and paired-end FASTQ files (R1 and R2)
        reference='data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna',
        # File paths for FASTQ files (R1, R2) for each sample
        fastq_trim_1="data/Genomic_data/{Species}/{Accession}/{Accession}_R1_trimmed.fastq",
        fastq_trim_2="data/Genomic_data/{Species}/{Accession}/{Accession}_R2_trimmed.fastq"
    output:
        # Output SAM, BAM, sorted BAM, and BAM index files
        sam="data/Genomic_data/{Species}/{Accession}/{Accession}_aligned_genome.sam",
        bam="data/Genomic_data/{Species}/{Accession}/{Accession}_aligned_genome.bam",
        sorted_bam="data/Genomic_data/{Species}/{Accession}/{Accession}_sorted_aligned_genome.bam",
        metrics="data/Genomic_data/{Species}/{Accession}/{Accession}_dedup_metrics.txt",
        dedup_sorted_bam="data/Genomic_data/{Species}/{Accession}/{Accession}_dedup_sorted_aligned_genome.bam",
        index_bam="data/Genomic_data/{Species}/{Accession}/{Accession}_dedup_sorted_aligned_genome.bam.bai"
    params:
        # Parameters for each sample's read group in BWA
        RGID=lambda wildcards: sample_map[wildcards.Accession],  # Read group identifier
        RGLB=lambda wildcards: sample_map[wildcards.Accession],  # Library identifier
        RGPL="ILLUMINA",       # Sequencing platform
        RGPM="HISEQ",          # Platform model
        RGSM=lambda wildcards: sample_map[wildcards.Accession],   # Sample name
    conda:
        # Specify the Conda environment.
        "environment.yml"
    threads: config["ncpus"]  # Number of threads for BWA
    shell:
        """
        # Align FASTQ reads to the reference genome, output SAM file
        bwa mem -M -R "@RG\\tID:{params.RGID}\\tLB:{params.RGLB}\\tPL:{params.RGPL}\\tPM:{params.RGPM}\\tSM:{params.RGSM}" {input.reference} {input.fastq_trim_1} {input.fastq_trim_2} > {output.sam}
               
        # Convert SAM to BAM
        samtools view -Sb {output.sam} > {output.bam}
        
        # Sort the BAM file
        samtools sort {output.bam} -o {output.sorted_bam}
        
        picard MarkDuplicates I={output.sorted_bam} O={output.dedup_sorted_bam} M={output.metrics} REMOVE_DUPLICATES=true
        
        # Index the sorted BAM file
        samtools index {output.dedup_sorted_bam}
        
    """

# Generated read group @RG entries
rule generate_rg:
    input:
        dedup_sorted_bam="data/Genomic_data/{Species}/{Accession}/{Accession}_dedup_sorted_aligned_genome.bam",
    output:
        rg_file="data/Genomic_data/{Species}/{Accession}/{Accession}_rg.txt",
    params:
        sample_id=lambda wildcards: sample_map[wildcards.Accession],
    shell:
        """
        # Generate @RG tags dynamically
        : > {output.rg_file}  # Clear existing file
        
        echo -e "@RG\tID:{params.sample_id}\tLB:{params.sample_id}\tPL:ILLUMINA\tPM:HISEQ\tSM:{params.sample_id}" >> {output.rg_file}
        """

# Combine all individual read group (@RG) files for each species
rule merge_rg:
    input:
        lambda wildcards: expand(
            "data/Genomic_data/{Species}/{Accession}/{Accession}_rg.txt", 
            Species=[wildcards.Species], 
            Accession=species_accessions[wildcards.Species]
        ),
    output:
        "data/Genomic_data/{Species}/merged_bam/{Species}_rg.txt",
    shell:
        """
        # Combine all individual @RG files into one
        mkdir -p $(dirname {output})
        cat {input} > {output}
        """

# Rule to merge individual BAM files into a single merged BAM file
rule merge_bam:
    input:
        lambda wildcards: expand(
            "data/Genomic_data/{Species}/{Accession}/{Accession}_dedup_sorted_aligned_genome.bam",
            Species=[wildcards.Species],
            Accession=species_accessions[wildcards.Species]
        )
    output:
        merged_bam="data/Genomic_data/{Species}/merged_bam/{Species}_merged.bam"
    conda:
        "environment.yml"
    threads: config["ncpus"]  # Number of threads for samtools merge
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.merged_bam})
        
        # Merge BAM files using samtools
        samtools merge -@ {threads} {output.merged_bam} {input}
        
        # Index the merged BAM file
        samtools index {output.merged_bam}
        """

# Calls variants, extracts SNPs, and applies quality filters using GATK for genotype analysis
rule call_genotypes:
    input:
        merged_bam=lambda wildcards: f"data/Genomic_data/{wildcards.Species}/merged_bam/{wildcards.Species}_merged.bam",
        reference="data/Genomic_data/Reference_Genome/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna"
    output:
        raw_vcf="data/Genomic_data/{Species}/merged_bam/{Species}_gatk_raw_variants.vcf",
        snps_vcf="data/Genomic_data/{Species}/merged_bam/{Species}_gatk_only_SNPs_variants.vcf",
        filtered_vcf="data/Genomic_data/{Species}/final_vcf/{Species}_gatk_filtered.vcf"
    conda:
        "environment.yml"
    threads: config["ncpus"]  # Number of threads for GATK
    shell:
        """
        # Create directory for final VCF files
        mkdir -p $(dirname {output.filtered_vcf})
        
        # Step 1: Call variants
        gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" HaplotypeCaller \
            -R {input.reference} \
            -I {input.merged_bam} \
            -O {output.raw_vcf}
        
        # Step 2: Extract only SNPs
        gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" SelectVariants \
            -R {input.reference} \
            -V {output.raw_vcf} \
            --select-type-to-include SNP \
            -O {output.snps_vcf}
        
        # Step 3: Filter SNPs
        gatk --java-options "-Xmx4g -XX:ParallelGCThreads={threads}" VariantFiltration \
            -R {input.reference} \
            -V {output.snps_vcf} \
            --filter-name QD_filter --filter-expression "QD < 2.0" \
            --filter-name FS_filter --filter-expression "FS > 60.0" \
            --filter-name MQ_filter --filter-expression "MQ < 40.0" \
            --filter-name SOR_filter --filter-expression "SOR > 10.0" \
            -O {output.filtered_vcf}
        """
