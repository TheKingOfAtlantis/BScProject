
import pandas as pd

from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation

def LoadRecord(feature, toFind = "CDS"):
    data = []  # We will store the data we generate here

    for shift in range(0, 1): # FIXME: With the fasta file we don't have access to the sequence around the genome
        seq = feature.seq # Extract the DNA sequence using the location

        # Need to ensure we have a stop codon, otherwise we don't care
        if seq[-3:] in ["TAA", "TGA", "TAG"]:
            data.append({                  # Record the following:
                "shift": shift,            # Shift we applied to find codon
                "gc": SeqUtils.GC123(seq), # GC (incl. GC123) of sequence (given the shift)
                "stop": str(seq[-3:])      # What stop codon we found
            })

    return (feature.id, data) # Pair the data with the record ID

gc = []
data = list(SeqIO.parse("data/genomes/human.fna", "fasta"))

import common, os

from multiprocessing import Pool
import numpy as np

def concatPreprocess(data): return { k:pd.DataFrame(v) for k,v in data }

if __name__ == "__main__":
    result = common.__loadParallel(LoadRecord, data, len(data), True, False)
    result = common.concat(result, concatPreprocess)
    result.to_json("data/gc/cds/human.json", orient="table")

# Humans have high GC due to high recombintation rates + GC Conversion
# Check for premature stops
