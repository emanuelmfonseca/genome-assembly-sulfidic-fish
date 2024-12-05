# Genomic Data Processing Pipeline for Sulfidic Fish Species

## Overview

This pipeline automates the genomic data processing workflow for sulfidic fish species, including downloading reference genomes, preparing paired-end FASTQ files, performing quality control, aligning sequences, and calling variants using GATK. It is implemented using [Snakemake](https://snakemake.readthedocs.io/en/stable/) for workflow management and utilizes Conda for dependency management.

---

## Key Features

1. **Download and Index Reference Genome**: Automates retrieval and preparation of reference genomes for downstream analysis.
2. **FASTQ Processing**: Downloads paired-end sequencing data, performs quality control, and trims reads.
3. **Alignment**: Maps reads to the reference genome using BWA and marks duplicates.
4. **Variant Calling**: Identifies SNPs using GATK and applies quality filtering.
5. **Scalability**: Configurable to handle multiple species and samples efficiently.

---

## Requirements

1. **Software**:
   - [Snakemake](https://snakemake.readthedocs.io/en/stable/)
   - [Conda](https://docs.conda.io/)
   - [BWA](http://bio-bwa.sourceforge.net/)
   - [GATK](https://gatk.broadinstitute.org/)
   - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
   - [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
   - [Samtools](http://www.htslib.org/)

2. **Input Data**:
   - Sample metadata in `Sample_Information.xml`.

3. **Configuration**:
   - A `config.yaml` file specifying parameters such as number of CPUs and file paths.

---

## Installation

1. Clone the repository:
   ```bash
   git clone genome-assembly-sulfidic-fish
   cd genome-assembly-sulfidic-fish
   ```

2. Set up the Conda environment:
   ```bash
   conda env create -f environment.yml
   conda activate genomic-pipeline
   ```

3. Ensure Snakemake is installed:
   ```bash
   conda install -c bioconda snakemake
   ```

---

## Usage

1. **Set Up Configuration**:
   Edit `config.yaml` to customize parameters such as the number of CPUs, input paths, and output directories.
   The only user configuration required is to specify the path to the fastq-dump executable from the [sratoolkit](https://github.com/ncbi/sra-tools/) in the config.yaml file.

2. **Run the Pipeline**:
   ```bash
   snakemake --use-conda --cores <number-of-cores>
   ```

3. **Visualize Workflow**:
   Generate a DAG of the workflow for visualization:
   ```bash
   snakemake --dag | dot -Tpng > workflow.png
   ```

---

## File Structure

```
data/
├── Genomic_data/
│   ├── Reference_Genome/
│   ├── <Species>/<Accession>/
│   ├── merged_bam/
│   └── final_vcf/
├── Sample_Information/
│   ├── Sample_Information.xml
│   ├── Extracted_Sample_Information.csv
scripts/
├── aux_func.py
config.yaml
environment.yml
Snakefile
```

---

## Rules Summary

1. **`download_reference`**: Downloads and decompresses the reference genome.
2. **`download_fastq`**: Retrieves paired-end FASTQ files.
3. **`check_quality_control`**: Generates quality reports using FastQC.
4. **`trimming_fastq`**: Trims FASTQ files with Trimmomatic.
5. **`index_reference`**: Indexes the reference genome.
6. **`align_fastqc`**: Aligns reads to the reference genome and marks duplicates.
7. **`call_genotypes`**: Calls SNPs using GATK and applies filtering.

---

## Outputs

1. **Quality Control Reports**:
   - Raw and trimmed FASTQ FastQC reports.
2. **Alignment Files**:
   - BAM and sorted BAM files for each sample.
3. **Variant Calls**:
   - Raw, SNP-only, and filtered VCF files.

---

## Notes

- Update the `config.yaml` file to match your dataset structure and computational resources.
- This pipeline is modular and can be extended to include additional analysis steps.
- Ensure that all required software is installed and accessible in your environment.

---

For questions or issues, please submit an issue on the repository or contact Emanuel M. Fonseca directly at [emanuelmfonseca@gmail.com](mailto:emanuelmfonseca@gmail.com).