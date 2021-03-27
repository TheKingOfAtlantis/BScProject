#!/bin/sh

# Make the output directories if they do not exist
mkdir -p data/genomes

# Check if we have unzipped the genomes
if [ ! -d "./data/genomes/archaea" ]; then
    unzip ./data/Archaea_filtered_genomes.zip -d ./data
    mv data/Archaea_filtered_genomes/ data/genomes/archaea
fi
if [ ! -d "./data/genomes/bacteria" ]; then
    unzip ./data/Bacterial_filtered_genomes.zip -d ./data
    mv ./data/filtered_genomes/ ./data/genomes/bacteria
fi

python scripts/qc/avaliableIDs.py        # Identifies the minimal identifier set
python scripts/qc/Protein.py             # Analysis genomes for non-conforming CDSs
python scripts/qc/Protein-Isolate.py     # Produce list of non-conforming CDSs (Psuedogene check)
python scripts/qc/GenomeVSGenes.py       # Plots no. of genes vs genome size
python scripts/qc/TranslationTable.py    # Identify translation table usage

