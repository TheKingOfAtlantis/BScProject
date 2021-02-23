import os
from zipfile import ZipFile

from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation

from common import loadGlob

def LoadRecord(file):
    record = next(SeqIO.parse(file, "embl")) # Load the record from the file
    data = [] # We will store the data we generate here

    # Filtering out the features which != CDS and iterate over remainder
    for feature in filter(lambda x: x.type in "CDS", record.features):
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

import json
if __name__ == "__main__":
    archaea_gene_gc = dict(loadGlob("data/genomes/archaea/*", LoadRecord))
    with open("data/gc/cds/archaea.json", 'w') as file: json.dump(archaea_gene_gc, file)
    del archaea_gene_gc

    bacteria_gene_gc = dict(loadGlob("data/genomes/bacteria/*", LoadRecord))
    with open("data/gc/cds/bacteria.json", 'w') as file: json.dump(bacteria_gene_gc, file)
    del bacteria_gene_gc
