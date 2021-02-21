from Bio import SeqIO, SeqUtils
from Bio.SeqFeature import FeatureLocation

from common import loadZip

def LoadRecord(file):
    record = next(SeqIO.parse(file, "embl"))
    data = []
    for feature in filter(lambda x: x.type in "CDS", record.features):
        for shift in range(0, 6): # Check across 2 codons
            loc   = feature.location + shift
            seq = loc.extract(record.seq)

            # Need to ensure we have a stop codon
            if seq[-3:] in ["TAA", "TGA", "TAG"]:
                data.append({
                    "shift": shift,
                    "gc": SeqUtils.GC123(seq),
                    "stop": str(seq[-3:])
                })
    return (record.id, data)

import json
if __name__ == "__main__":
    archaea_gene_gc = dict(loadZip("data/Archaea_filtered_genomes.zip", LoadRecord))
    with open("data/archaea_gene_gc.json", 'w') as file: json.dump(archaea_gene_gc, file)
    del archaea_gene_gc

    bacteria_gene_gc = dict(loadZip("data/Bacteria_filtered_genomes.zip", LoadRecord))
    with open("data/bacteria_gene_gc.json", 'w') as file: json.dump(bacteria_gene_gc, file)
    del bacteria_gene_gc
