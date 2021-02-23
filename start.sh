#!/bin/sh

# Check if we have unzipped the genomes
if [ ! -d "data/archaea_filtered_genomes"]
    unzip data/Archaea_filtered_genomes.zip
    mv data/Archaea_filtered_genomes/ data/archaea_filtered_genomes/
if [ ! -d "data/bacterial_filtered_genomes"]
    unzip data/Bacterial_filtered_genomes.zip
    mv data/filtered_genomes/ data/bacterial_filtered_genomes/

# Make the output directories if they do not exist
mkdir -p data/gc/cds data/gc/ncRNA

py scripts/genomic-GC.py   # Process files to generate GC
py scripts/gene-GC.py      # Process files to generate GC of CDS
py scripts/gene-GC-plot.py # Produce the plots for cds
