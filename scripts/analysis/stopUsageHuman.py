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

gffutils.constants.always_return_list = False

annotation = gffutils.FeatureDB("data/genomes/GRCh38_latest_genomic.gff3.sqlite")
sequence = {seq.id:seq for seq in tqdm(
    SeqIO.parse("data/genomes/GRCh38_latest_genomic.fna", "fasta"),
    desc = "Loading FASTA Sequences"
)}

# Select for CDSs and process
representativeCDSs = pd.read_csv("data/qc/proteins/human.csv")

CDSs = list(annotation[id] for id in tqdm(
    representativeCDSs.cds,
    desc  = "Loading CDS Annotations",
    total = len(representativeCDSs.cds)
))

def processFeature(cds):
    chromosome = sequence[cds.chrom]
    feature    = to_seqfeature(cds)
    location   = feature.location

    data = []
    for shift in range(0, 6):
        loc = location + shift            # Add the shift to the "in-frame" ORF location
        seq = loc.extract(chromosome.seq) # Extract the DNA sequence using the location

        # Need to ensure we have a stop codon, otherwise we don't care
        if seq[-3:] in ["TAA", "TGA", "TAG", "TAC"]:
            data.append({                  # Record the following:
                "shift": shift,            # Shift we applied to find codon
                "gc": SeqUtils.GC123(seq), # GC (incl. GC123) of sequence (given the shift)
                "stop": str(seq[-3:])      # What stop codon we found
            })
    return (feature.id, data) # Pair the data with the record ID

results = Parallel.loadParallel(
    processFeature, CDSs,
    len(CDSs),
    chunkSize = 500,
    desc = "Stop Usage in Human Genome"
)
def concat(x): return { k:pd.DataFrame(v) for k,v in x }
result = Parallel.concat(results, concat)
result.to_json("data/gc/cds/human.json", orient="table")
# Humans have high GC due to high recombintation rates + GC Conversion
# Check for premature stops
