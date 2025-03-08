### pseuPIRA.py by Rohan Maddamsetti.

themisto, and minimap2 must be in the $PATH.

The command-line interface is modelled after breseq (https://github.com/barricklab/breseq)

Usage: python pseuPIRA.py -r reference.gbff.gz reads1.fastq [reads2.fastq.gz ...]

Run pseuPIRA for estimating plasmid copy number from microbial genome sequencing data.

FASTQ read files (which may be gzipped) are input as the last unnamed argument(s).

Allowed Options
  -h,--help                        Produce help message showing advanced options
  -r,--reference <arg>             File containing reference genome sequence in gzipped Genbank (*.gbff.gz) format. (REQUIRED)
  -n,--name <arg>                  Human-readable name of the analysis run for output (DEFAULT=<none>)
  -j,--num-processors <arg>        Number of processors to use in multithreaded steps (DEFAULT=1)
  -o,--output <arg>                Path to pseuPIRA output (DEFAULT=.)


## Example data are found in ./data/RHBSTW-00316

The test data comes from NCBI BioProject PRJNA605147
https://www.ncbi.nlm.nih.gov/bioproject/PRJNA605147/

RefSeq ID: GCF_013742375.1  
Assembly level: Complete genome  
1 chromosome and 9 plasmids  
BioSample: SAMN15148572  
strain: RHBSTW-00316  
Species: Enterobacter hormaechei  

The reference genome is found here:
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_013742375.1/
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/742/375/GCF_013742375.1_ASM1374237v1/

The Illumina sequencing data is found here:
https://www.ncbi.nlm.nih.gov/sra/SRX8493146[accn]


To test pseuPIRA.py on the sample data, you have to download the SRR1194861.fastq.gz from the link on this page:
https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR11948691&display=download  

Save these data in data/RHBSTW-00316/SRR11948691.fastq.gz.

## Running pseuPIRA.py on the example data:
python3 pseuPIRA.py -o ../results/RHBSTW-00316 -r ../data/RHBSTW-00316/GCF_013742375.1_ASM1374237v1_genomic.gbff.gz ../data/RHBSTW-00316/SRR11948691.fastq.gz

