#!/bin/zsh

## for a fair comparison, running on uncompressed fastq data

echo "pseuPIRA.py in quick mode\n"
## quick mode (no PIRA, just pseudoalignment)
/usr/bin/time -l -h -p python pseuPIRA.py -q -o ../results/quick_RHBSTW-00316 -r ../data/RHBSTW-00316/GCF_013742375.1_ASM1374237v1_genomic.gbff.gz ../data/RHBSTW-00316/SRR11948691.fastq

echo "\npseuPIRA in full mode\n"
## full mode (pseudoalignment and PIRA
/usr/bin/time -l -h -p python pseuPIRA.py -o ../results/full_RHBSTW-00316 -r ../data/RHBSTW-00316/GCF_013742375.1_ASM1374237v1_genomic.gbff.gz ../data/RHBSTW-00316/SRR11948691.fastq

echo "\ncoverm for comparison\n"
## run coverm for comparison
/usr/bin/time -l -h -p coverm genome --single ../data/RHBSTW-00316/SRR11948691.fastq --genome-fasta-files ../data/RHBSTW-00316/fasta/chromosome1.fna ../data/RHBSTW-00316/fasta/plasmid2.fna ../data/RHBSTW-00316/fasta/plasmid3.fna ../data/RHBSTW-00316/fasta/plasmid4.fna ../data/RHBSTW-00316/fasta/plasmid5.fna ../data/RHBSTW-00316/fasta/plasmid6.fna ../data/RHBSTW-00316/fasta/plasmid7.fna ../data/RHBSTW-00316/fasta/plasmid8.fna ../data/RHBSTW-00316/fasta/plasmid9.fna ../data/RHBSTW-00316/fasta/plasmid10.fna -t 4 -m mean relative_abundance covered_fraction -o ../results/RHBSTW-00316-output_coverm.tsv

echo "\n"

## try for a really big Illumina dataset.

echo "\ntrying a really big dataset"
echo "pseuPIRA.py in quick mode\n"
## quick mode (no PIRA, just pseudoalignment)
/usr/bin/time -l -h -p python pseuPIRA.py -q -o ../results/quick_PRJNA280982 -r ../data/PRJNA280982/GCF_002285515.1_ASM228551v1_genomic.gbff.gz ../data/PRJNA280982/SRR1974308_1.fastq ../data/PRJNA280982/SRR1974308_2.fastq

echo "\npseuPIRA in full mode\n"
## full mode (pseudoalignment and PIRA
/usr/bin/time -l -h -p python pseuPIRA.py -o ../results/full_PRJNA280982 -r ../data/PRJNA280982/GCF_002285515.1_ASM228551v1_genomic.gbff.gz ../data/PRJNA280982/SRR1974308_1.fastq ../data/PRJNA280982/SRR1974308_2.fastq

echo "\ncoverm for comparison\n"
## run coverm for comparison
/usr/bin/time -l -h -p coverm genome --coupled ../data/PRJNA280982/SRR1974308_1.fastq ../data/PRJNA280982/SRR1974308_2.fastq --genome-fasta-files ../data/PRJNA280982/chromosome01.fna ../data/PRJNA280982/plasmid02.fna -t 4 -m mean relative_abundance covered_fraction -o ../results/PRJNA280982-output_coverm.tsv
