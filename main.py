import io, itertools
from zipfile import ZipFile

import pandas as pd
from Bio import SeqIO, SeqUtils


def calculateGC(record):
    print(record.id)
    data = {
        "gc": SeqUtils.GC(record.seq),
        "features": []
    }
    # Filter out everything which isn't a coding sequence
    # Then perform calculations on those features
    for feature in filter(lambda x: "CDS" in x.type, record.features):
        # print("\t",feature.qualifiers["locus_tag"])
        # Now we use the position to extract the DNA substring
        # Then Calculate GC content on that segment
        location = feature.location
        seq      = record.seq[location.start:location.end]
        gc       = SeqUtils.GC123(seq)

        idTypes = feature.qualifiers.keys()
        if("locus_tag" in idTypes): id = feature.qualifiers["locus_tag"]
        elif("db_xref" in idTypes): id = list(filter(lambda x: "EnsemblGenomes-Gn" in x, feature.qualifiers["db_xref"]))

        data["features"].append({
            "gene_id": itertools.chain(id),
            "gc": {
                "total": gc[0],
                "third": gc[-1]
            }
        })
    return (record.id, data)

def GCFromZip(zipfile):
    out = {}
    with ZipFile(zipfile) as zip:
        i = 0;
        for path in zip.namelist():
            i += 1
            print("{}/{}".format(i, len(zip.namelist())))

            file    = io.TextIOWrapper(zip.open(path))
            records = SeqIO.parse(file, "embl")

            for record in records:
                (id, data) = calculateGC(record)
                out[id] = data
    return out

import json

archaea_gc  = GCFromZip("data/Archaea_filtered_genomes.zip")
with open("data/archaea_gc.json", "w")  as js: json.dump(archaea_gc,  js)
del archaea_gc

bacteria_gc = GCFromZip("data/bacterial_filtered_genomes.zip")
with open("data/bacteria_gc.json", "w") as js: json.dump(bacteria_gc, js)
del bacteria_gc
