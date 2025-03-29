#!/bin/zsh

## run coverm for comparison
/usr/bin/time -l -h -p coverm genome --single ../data/RHBSTW-00316/SRR11948691.fastq.gz --genome-fasta-files ../data/RHBSTW-00316/fasta/chromosome1.fna ../data/RHBSTW-00316/fasta/plasmid2.fna ../data/RHBSTW-00316/fasta/plasmid3.fna ../data/RHBSTW-00316/fasta/plasmid4.fna ../data/RHBSTW-00316/fasta/plasmid5.fna ../data/RHBSTW-00316/fasta/plasmid6.fna ../data/RHBSTW-00316/fasta/plasmid7.fna ../data/RHBSTW-00316/fasta/plasmid8.fna ../data/RHBSTW-00316/fasta/plasmid9.fna ../data/RHBSTW-00316/fasta/plasmid10.fna -t 4 -m mean relative_abundance covered_fraction -o ../results/output_coverm.tsv

echo "\n"

## quick mode (no PIRA, just pseudoalignment)
/usr/bin/time -l -h -p python pseuPIRA.py -q -o ../results/quick_RHBSTW-00316 -r ../data/RHBSTW-00316/GCF_013742375.1_ASM1374237v1_genomic.gbff.gz ../data/RHBSTW-00316/SRR11948691.fastq.gz

echo "\n"

## full mode (pseudoalignment and PIRA
/usr/bin/time -l -h -p python pseuPIRA.py -o ../results/full_RHBSTW-00316 -r ../data/RHBSTW-00316/GCF_013742375.1_ASM1374237v1_genomic.gbff.gz ../data/RHBSTW-00316/SRR11948691.fastq.gz

