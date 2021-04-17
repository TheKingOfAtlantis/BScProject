#!/bin/bash

# If something goes wrong we want to exit immediately
# Otherwise might cause downstream problems and error messeage would get lost
set -e

header=$(tput bold; tput smul)
headerend=$(tput sgr0)

# Check if we have the genomes
if [ ! -d "./data/genomes" ]; then
    bash buildDataset.sh -y
fi

echo
echo "${header}Running QC Scripts${headerend}"
echo

python scripts/qc/avaliableIDs.py                # Identifies the minimal identifier set
python scripts/qc/proteinCheck.py                # Analysis genomes for non-conforming CDSs
python scripts/qc/proteinIsolate.py              # Produce list of non-conforming CDSs (Psuedogene check)
python scripts/qc/genomeVsGenes.py               # Plots no. of genes vs genome size
python scripts/qc/translationTable.py            # Identify translation table usage

echo
echo "${header}Running Analysis Scripts${headerend}"
echo

python scripts/analysis/genomicGC.py             # Process files to generate GC of genomes
python scripts/analysis/stopUsageProkarytoes.py  # Process files to generate GC3 & Stop usage data
python scripts/analysis/stopUsageHuman.py        # Process human CDSs for GC3 & Stop usage
python scripts/analysis/stopUsagePlot.py         # Produce the plots for GC3 & stop usage
python scripts/analysis/stopUsagePlotLogistic.py # Logistic regression on human data
