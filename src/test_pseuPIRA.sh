#!/bin/bash

## test_pseuPIRA.sh by Rohan Maddamsetti

echo "running test suite for pseuPIRA.py"
echo "##########################################################################"
DATAFILE="../data/RHBSTW-00316/SRR11948691.fastq.gz"

if [ ! -f "$DATAFILE" ]; then
    echo "Error: File $DATAFILE not found."
    exit 1
fi

echo "Read file $DATAFILE exists. Continuing with testing script..."

echo "##########################################################################"
echo "testing quick failure when data is not present:"
python pseuPIRA.py -o ../results/RHBSTW-00316 -r ../data/RHBSTW-00316/GCF_013742375.1_ASM1374237v1_genomic.gbff.gz ../data/RHBSTW-00316/nonexistent123.fastq.gz

echo "##########################################################################"
echo "testing quick failure when data is wrong format:"
python pseuPIRA.py -o ../results/RHBSTW-00316 -r ../data/RHBSTW-00316/GCF_013742375.1_ASM1374237v1_genomic.gbff.gz ../data/RHBSTW-00316/fasta/GCF_013742375.1_ASM1374237v1_genomic.fna.gz

echo "##########################################################################"
echo "testing proper operation:"
python pseuPIRA.py -o ../results/RHBSTW-00316 -r ../data/RHBSTW-00316/GCF_013742375.1_ASM1374237v1_genomic.gbff.gz ../data/RHBSTW-00316/SRR11948691.fastq.gz

echo "test script complete."
