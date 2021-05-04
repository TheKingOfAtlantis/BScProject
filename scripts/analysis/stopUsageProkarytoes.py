from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation

import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from common import Filesystem, Parallel, getID

import pandas as pd

invalidCDS = pd.read_csv("data/qc/proteins/failed.csv")
def isInvalidFeature(genome, feature):
    rows = invalidCDS[invalidCDS.genome == genome]
    rows = rows[rows.id == getID(feature)]
    return not rows.empty

def processCDS(record, feature):
    if(isInvalidFeature(record.id, feature)):
        return [] # Skip this feature

    data = []

    seq = feature.location.extract(record.seq)
    gc  = SeqUtils.GC123(seq)

    # We want to explore value across frameshifts
    # So shift the frame over 2 codons worth of nucleotides
    for shift in range(0, 6): # Check across 2 codons
        loc = feature.location + shift # Add the shift "in-frame" ORF location
        seq = loc.extract(record.seq)  # Extract the DNA sequence using the location

        # Need to ensure we have a stop codon, otherwise we don't care
        if seq[-3:] in ["TAA", "TGA", "TAG", "TAC"]:
            data.append({              # Record the following:
                "shift": shift,        # Shift we applied to find codon
                "gc":    gc,           # GC (incl. GC123) of sequence
                "stop":  str(seq[-3:]) # What stop codon we found
            })
    return data

def processRNA(record, feature):
    data = []

    seq = feature.location.extract(record.seq)
    gc  = SeqUtils.GC123(seq)

    for stop in ["TAA", "TGA", "TAG", "TAC"]:
        count = seq.count(stop) # Count number of non-overlapping instances of each codon
        data.extend([{          # Record the following:
            "shift": "all",     # Shift we applied to find codon
            "gc": gc,           # GC (incl. GC123) of sequence
            "stop": str(stop)   # What stop codon we found
        }] * count)

    # We want to explore value across frameshifts
    # So shift the frame over 2 codons worth of nucleotides
    for shift in range(0, 3): # Check across 2 codons
        loc = feature.location + shift # Add the shift to "in-frame" location
        seq = loc.extract(record.seq)  # Extract the DNA sequence using the location

        for i in range(0, len(seq) - 3, 3):
            codon = seq[i:i+3]
            if codon in ["TAA", "TGA", "TAG", "TAC"]:
                data.append({          # Record the following:
                    "shift": shift,    # Shift we applied to find codon
                    "gc": gc,          # GC (incl. GC123) of sequence
                    "stop": str(codon) # What stop codon we found
                })
    return data

def LoadRecord(file, toFind):
    record = next(SeqIO.parse(file, "embl")) # Load the record from the file
    data = [] # We will store the data we generate here

    # Filtering out the features which != CDS and iterate over remainder
    for feature in filter(lambda x: x.type in toFind, record.features):

        # Check feature against the list of invalid features
        # Found from QCing
        if(feature.type == "CDS"):    data.extend(processCDS(record, feature))
        elif(feature.type == "tRNA"): data.extend(processRNA(record, feature))

    return (record.id, data) # Pair the data with the record ID

def concatPreprocess(data): return { k:pd.DataFrame(v) for k,v in data }

for geneType in ["CDS", "tRNA"]:
    # Make sure we have the path to export to (including parents)
    Filesystem.mkdir(f"data/gc/{geneType.lower()}")

    archaea_gene_gc = Filesystem.loadGlob(
        "data/genomes/archaea/*", LoadRecord,
        toFind=geneType,
        desc = f"Archaea {geneType} Stop Usage"
    )
    Parallel.concat(archaea_gene_gc, concatPreprocess).to_json(f"data/gc/{geneType.lower()}/archaea.json", orient="table")

    bacteria_gene_gc = Filesystem.loadGlob(
        "data/genomes/bacteria/*", LoadRecord,
        toFind=geneType,
        desc = f"Bacteria {geneType} Stop Usage"
    )
    Parallel.concat(bacteria_gene_gc, concatPreprocess).to_json(f"data/gc/{geneType.lower()}/bacteria.json", orient="table")
