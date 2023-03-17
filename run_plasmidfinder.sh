#!/bin/bash

for fasta in FASTA/*.fasta; do
    ID=$(echo $fasta | cut -d '.' -f 1 | cut -d '/' -f 2)
    # echo $ID
    mkdir "PFOut/${ID}"
    python3 plasmidfinder.py -i $fasta -o "PFOut/${ID}" -p plasmidfinder_db
done