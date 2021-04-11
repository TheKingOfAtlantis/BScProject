import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))

from IPython import get_ipython
if(get_ipython() != None):
    import os
    os.chdir("../..")

from common import Parallel
from functools import partial
from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation

import gffutils
from gffutils.biopython_integration import to_seqfeature

import pandas as pd
from tqdm.auto import tqdm

# Load annotations and sequences

annotation = gffutils.FeatureDB("data/genomes/GRCh38_latest_genomic.gff3.sqlite")
sequence = list(tqdm(
    SeqIO.parse("data/genomes/GRCh38_latest_genomic.fna", "fasta"),
    desc = "Loading FASTA Sequences"
))
# Select for CDSs and process
CDSs = list(tqdm(
    annotation.features_of_type("CDS"),
    desc  = "Loading CDS Annotations",
    total = annotation.count_features_of_type("CDS")
))

def processFeature(cds):
    chromosome = next(filter(lambda x: x.id == cds.chrom, sequence))
    feature    = to_seqfeature(cds)

    data = []
    for shift in range(0, 6):
        loc = feature.location + shift # Add the shift the "in-frame" ORF location
        seq = loc.extract(chromosome.seq)  # Extract the DNA sequence using the location

        # Need to ensure we have a stop codon, otherwise we don't care
        if seq[-3:] in ["TAA", "TGA", "TAG"]:
            data.append({                  # Record the following:
                "shift": shift,            # Shift we applied to find codon
                "gc": SeqUtils.GC123(seq), # GC (incl. GC123) of sequence (given the shift)
                "stop": str(seq[-3:])      # What stop codon we found
            })
    return (feature.id, pd.DataFrame.from_dict(data)) # Pair the data with the record ID

results = Parallel.loadParallel(
    processFeature, CDSs,
    len(CDSs),
    desc = "Stop Usage in Human Genome"
)

# Humans have high GC due to high recombintation rates + GC Conversion
# Check for premature stops
