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
CDSs = list(tqdm(
    annotation.features_of_type("CDS"),
    desc  = "Loading CDS Annotations",
    total = annotation.count_features_of_type("CDS")
))

mitochondrialGenomes = [id for id,seq in sequence.items() if "mitochondrion" in seq.description]
def filterCDS(cds):
    # Presence of these attributes means that they are true (checked that this is true)
    toRemove = [
        "pseudo", # We don't want pseudogenes => They don't code for anything
        "partial" # Missing start or stop codons
    ]

    for attrib in toRemove:
        if(attrib in cds.attributes):
            return None

    # Strange occurance but some CDSs have the same start and end
    # We'll drop these
    if(cds.start == cds.end):
        return None

    # Removing mitochondrial genomes from our list
    if(cds.seqid in mitochondrialGenomes):
        return None

    return cds

def produceCDSPosition(cds):
    # CDS position span the region which contains the CDS
    # However this may contain multiple introns/exons
    # So we need to find the position of all introns it contains
    # and join them
    db = gffutils.FeatureDB("data/genomes/GRCh38_latest_genomic.gff3.sqlite")
    exons = db.region(cds, featuretype="exon")

    positions = [{
        "start":  exon.start,
        "end":    exon.end,
        "strand": exon.strand
    } for exon in exons]

    for pos in positions:
        if(pos["strand"] == "+"):   pos["strand"] = +1
        elif(pos["strand"] == "-"): pos["strand"] = -1
        elif(pos["strand"] == "?"): pos["strand"] = 0
        else: pos["strand"] = None

    from operator import itemgetter
    positions = sorted(positions, key=itemgetter("start"))

    # The start of the CDS will be contained w/i a exon
    # So we limit the start of the first exon to start w/ the CDS
    # Similarly end will be contained w/i the last exon
    positions[0]["start"] = cds.start
    positions[-1]["end"]  = cds.end

    from Bio.SeqFeature import FeatureLocation, CompoundLocation
    return FeatureLocation(**positions[0]) if len(positions) == 1 else CompoundLocation([
        FeatureLocation(**pos) for pos in positions
    ])

def processFeature(cds):
    chromosome = sequence[cds.chrom]
    feature    = to_seqfeature(cds)
    location   = feature.location

    data = []
    correction = 0
    for shift in range(0, 6):
        # Some inconsistency w/ CDSs
        # Stop codon is contained in the sequences of some but not all
        #
        # So we check, if in the in-frame shift there is no stop codon then we
        # hopefully correct for this by shifting 1 codon along
        if(shift == 0):
            seq = location.extract(chromosome.seq)
            if seq[-3:] not in ["TAA", "TGA", "TAG"]:
                correction += 3

        loc = location + shift + correction # Add the shift to the "in-frame" ORF location
        seq = loc.extract(chromosome.seq)   # Extract the DNA sequence using the location

        # Need to ensure we have a stop codon, otherwise we don't care
        if seq[-3:] in ["TAA", "TGA", "TAG", "TAC"]:
            data.append({                  # Record the following:
                "shift": shift,            # Shift we applied to find codon
                "gc": SeqUtils.GC123(seq), # GC (incl. GC123) of sequence (given the shift)
                "stop": str(seq[-3:])      # What stop codon we found
            })
    return (feature.id, data) # Pair the data with the record ID

CDSs = Parallel.loadParallel(
    filterCDS, CDSs,
    len(CDSs),
    chunkSize = 1000
    desc = "Checking CDSs in Human Genome"
)


results = Parallel.loadParallel(
    processFeature, CDSs,
    chunkSize = 500,
    len(CDSs),
    desc = "Stop Usage in Human Genome"
)
def concat(x): return { k:pd.DataFrame(v) for k,v in x }
result = Parallel.concat(results, concat)
result.to_json("data/gc/cds/human.json", orient="table")
# Humans have high GC due to high recombintation rates + GC Conversion
# Check for premature stops
