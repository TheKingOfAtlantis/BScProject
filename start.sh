#!/bin/sh

# Make the output directories if they do not exist
mkdir -p data/gc/cds data/gc/ncRNA data/genomes
mkdir -p plot/gc/cds

# Check if we have unzipped the genomes
if [ ! -d "./data/genomes/archaea" ]; then
    unzip ./data/Archaea_filtered_genomes.zip -d ./data
    mv data/Archaea_filtered_genomes/ data/genomes/archaea
fi
if [ ! -d "./data/genomes/bacteria" ]; then
    unzip ./data/Bacterial_filtered_genomes.zip -d ./data
    mv ./data/filtered_genomes/ ./data/genomes/bacteria
fi

python scripts/genomic-GC.py   # Process files to generate GC
python scripts/gene-GC.py      # Process files to generate GC of CDS
python scripts/gene-GC-plot.py # Produce the plots for cds
