# pseuPIRA

**pseuPIRA.py** is a pipeline for estimating plasmid copy number (PCN) from microbial genome sequencing data. It combines **themisto** for pseudoalignment and **minimap2** for multiread mapping. This project uses **uv** to manage dependencies and packaging.

---

## Overview

pseuPIRA.py processes a gzipped GenBank reference file to extract replicon sequences (e.g., chromosomes and plasmids) and then:

1. **Extracts Replicon Sequences:**  
   - Reads a GenBank file (`*.gbff.gz`) and creates individual FASTA files for each replicon with custom headers containing metadata.

2. **Builds an Index & Pseudoaligns Reads:**  
   - Uses **themisto** to build an index from the replicon FASTA files.
   - Pseudoaligns sequencing reads (FASTQ format) to count how many reads map to each replicon.

3. **Estimates Plasmid Copy Numbers (Naïve):**  
   - Computes initial PCN estimates by normalizing read counts by replicon lengths.

4. **Refines PCN Estimates (PIRA):**  
   - Filters reads that map to multiple replicons.
   - Realigns ambiguous reads using **minimap2**.
   - Applies a Probabilistic Iterative Read Assignment (PIRA) algorithm to update and refine PCN estimates.

---

## Requirements

- **External Tools:**  
  - `themisto` – must be installed and available in your `$PATH`
  - `minimap2` – must be installed and available in your `$PATH`

- **Python Libraries:**  
  - Biopython
  - Polars
  - HTSeq
  - NumPy
  - argparse (standard library)

---

## Dependency Management with uv

This project uses **uv** for dependency management and packaging, similar to Poetry. The configuration is defined in the `pyproject.toml` file.

### Getting Started with uv

1. **Install uv:**

   ```bash
   pip install uv
   ```

2. **Install Project Dependences:**
In the project directory, run: 

  ```bash
  uv install
  ```
This will install all the required dependencies as specified in the pyproject.toml.

3. **Activate the Virtual Environment:**
uv automatically creates a virtual environment (in .venv/) during installation. To use it:

```bash
source .venv/bin/activate  # Linux/MacOS
```
Note: If you use uv run for executing scripts, it will automatically handle the environment. For manual execution (python file.py), activation ensures you're using the correct dependencies.


## Usage 
Run the pipeline with: 
```bash
python src/pseuPIRA.py -r <reference.gbff.gz> <reads1.fastq> [<reads2.fastq.gz> ...]
```
## Command-Line Options

- `-r`, `--reference`  
  **(Required)** Path to the gzipped GenBank file (`*.gbff.gz`) containing the reference genome.

- `-j`, `--num-processors`  
  **(Optional)** Specifies the number of threads/processors to use.  
  *Note: This option is currently not functional; the code is hardcoded to use 4 threads.*

- `-o`, `--output`  
  **(Optional)** Path to the output directory.  
  *Default: `../results/test-run`*

- `-q`, `--quick`  
  **(Optional)** Quick mode. Runs the pseudoalignment and naive PCN estimation steps but skips the PIRA refinement.

- `reads`  
  **(Positional)** One or more FASTQ read files (supporting gzipped or plain FASTQ).

---

## Example Data

Example data for testing pseuPIRA can be found in the `./data/RHBSTW-00316` directory. This dataset originates from NCBI BioProject [PRJNA605147](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA605147).

### Test Data Details

- **Reference Genome**:
  - **RefSeq ID**: GCF_013742375.1
  - **Assembly Level**: Complete genome
  - **Replicon Composition**: 1 chromosome and 9 plasmids
  - **BioSample**: [SAMN15148572](https://www.ncbi.nlm.nih.gov/biosample/SAMN15148572)
  - **Strain**: RHBSTW-00316
  - **Species**: *Enterobacter hormaechei*

- **Reference Genome Files**:
  - **NCBI Genome Datasets**:  
    [View Dataset](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013742375.1/)
  - **FTP Link**:  
    [Download from FTP](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/742/375/GCF_013742375.1_ASM1374237v1/)

- **Sequencing Data**:
  - **Illumina Data Accession**: [SRX8493146](https://www.ncbi.nlm.nih.gov/sra/SRX8493146)
  - **Required FASTQ File**:  
    To test the pipeline, download the file `SRR11948691.fastq.gz` from the [NCBI Trace Archive for SRR11948691](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR11948691&display=metadata).

### Setup Instructions

1. **Download the Data**:
   - Save the reference genome file (e.g., `GCF_013742375.1_ASM1374237v1_genomic.gbff.gz`) to the `./data/RHBSTW-00316` directory.
   - Download and save the FASTQ file as:  
     `./data/RHBSTW-00316/SRR11948691.fastq.gz`

2. **Run pseuPIRA.py on the Example Data**:  
   Open a terminal in the project directory and execute:
   ```bash
   python pseuPIRA.py -o ../results/RHBSTW-00316 -r ../data/RHBSTW-00316/GCF_013742375.1_ASM1374237v1_genomic.gbff.gz ../data/RHBSTW-00316/SRR11948691.fastq.gz
   ```

