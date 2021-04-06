from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation

import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from common import Filesystem, Parallel, getID

import pandas as pd

invalidCDS = pd.read_csv("data/qc/proteins/failed.csv").drop(columns=["Unnamed: 0", "index"])
def isInvalidFeature(genome, feature):
    rows = invalidCDS[invalidCDS["genome"] == genome]
    rows = rows[rows["id"] == getID(feature)]
    return rows.empty

def LoadRecord(file, toFind):
    record = next(SeqIO.parse(file, "embl")) # Load the record from the file
    data = [] # We will store the data we generate here

    # Filtering out the features which != CDS and iterate over remainder
    for feature in filter(lambda x: x.type in toFind, record.features):

        # Check feature against the list of invalid features
        # Found from QCing
        # FIXME: I seem to be causing issues with the results just being random noise effectively
        # if(feature.type == "CDS"):
        #     if(isInvalidFeature(record.id, feature)):
        #         continue # Skip this feature

        # We want to explore value across frameshifts
        # So shift the frame over 2 codons worth of nucleotides
        for shift in range(0, 6): # Check across 2 codons
            loc   = feature.location + shift # Add the shift the "in-frame" ORF location
            seq = loc.extract(record.seq)    # Extract the DNA sequence using the location

            # Need to ensure we have a stop codon, otherwise we don't care
            if seq[-3:] in ["TAA", "TGA", "TAG", "TAC"]:
                data.append({                  # Record the following:
                    "shift": shift,            # Shift we applied to find codon
                    "gc": SeqUtils.GC123(seq), # GC (incl. GC123) of sequence (given the shift)
                    "stop": str(seq[-3:])      # What stop codon we found
                })
    return (record.id, data) # Pair the data with the record ID

def concatPreprocess(data): return { k:pd.DataFrame(v) for k,v in data }
if __name__ == "__main__":
    toFind = ["CDS", "tRNA"]

    # Make sure we have the path to export to (including parents)
    Filesystem.mkdir(f"data/gc/{geneType.lower()}")

    archaea_gene_gc = Filesystem.loadGlob("data/genomes/archaea/*", LoadRecord, toFind=toFind, desc = f"Processing Archaea {geneType}")
    Parallel.concat(archaea_gene_gc, concatPreprocess).to_json(f"data/gc/{geneType.lower()}/archaea.json", orient="table")

    bacteria_gene_gc = Filesystem.loadGlob("data/genomes/bacteria/*", LoadRecord, toFind=toFind, desc = f"Processing Bacteria {geneType}")
    Parallel.concat(bacteria_gene_gc, concatPreprocess).to_json(f"data/gc/{geneType.lower()}/bacteria.json", orient="table")


# def LoadRecord(file, toFind = "CDS"):
#     record = next(SeqIO.parse(file, "embl")) # Load the record from the file
#     record.
#  if __name__ == "__main__":
#     archaea_gene_gc = loadGlob("data/genomes/*/*", LoadRecord)
