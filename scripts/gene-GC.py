import os
from zipfile import ZipFile

from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation

from common import loadGlob, concat

import pandas as pd

def LoadRecord(file, toFind = "CDS"):
    record = next(SeqIO.parse(file, "embl")) # Load the record from the file
    data = [] # We will store the data we generate here

    # Filtering out the features which != CDS and iterate over remainder
    for feature in filter(lambda x: x.type in toFind, record.features):
        # We want to explore value across frameshifts
        # So shift the frame over 2 codons worth of nucleotides
        for shift in range(0, 6): # Check across 2 codons
            loc   = feature.location + shift # Add the shift the "in-frame" ORF location
            seq = loc.extract(record.seq)    # Extract the DNA sequence using the location

            # Need to ensure we have a stop codon, otherwise we don't care
            if seq[-3:] in ["TAA", "TGA", "TAG"]:
                data.append({                  # Record the following:
                    "shift": shift,            # Shift we applied to find codon
                    "gc": SeqUtils.GC123(seq), # GC (incl. GC123) of sequence (given the shift)
                    "stop": str(seq[-3:])      # What stop codon we found
                })
    return (record.id, data) # Pair the data with the record ID

def concatPreprocess(data): return { k:pd.DataFrame(v) for k,v in data }

import pathlib
if __name__ == "__main__":

    toFind = ["CDS", "tRNA"]

    for geneType in toFind:
        # Make sure we have the path to export to (including parents)
        pathlib.Path(f"data/gc/{geneType.lower()}").mkdir(parents=True, exist_ok=True)

        print(f"{geneType}:")
        archaea_gene_gc = loadGlob("data/genomes/archaea/*", LoadRecord, extra=toFind)
        concat(archaea_gene_gc, concatPreprocess).to_json(f"data/gc/{geneType.lower()}/archaea.json", orient="table")

        bacteria_gene_gc = loadGlob("data/genomes/bacteria/*", LoadRecord, extra=toFind)
        concat(bacteria_gene_gc, concatPreprocess).to_json(f"data/gc/{geneType.lower()}/bacteria.json", orient="table")
