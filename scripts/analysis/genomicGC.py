import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))

from common import Filesystem
from Bio import SeqIO, SeqUtils
import pandas as pd

def CalculateGC(file):
    record = next(SeqIO.parse(file, "embl"))
    return (record.id, SeqUtils.GC(record.seq))

import json, pathlib
if __name__ == "__main__":
    Filesystem.mkdir("data/gc/")

    archaea_gc = dict(Filesystem.loadGlob("data/genomes/archaea/*", CalculateGC, desc = "Archaea Genome GC"))
    archaea_gc = pd.DataFrame.from_dict(
        archaea_gc, orient="index", columns=["gc"]
    )

    bacteria_gc = dict(Filesystem.loadGlob("data/genomes/bacteria/*", CalculateGC, desc = "Bacteria Genome GC"))
    bacteria_gc = pd.DataFrame.from_dict(
        bacteria_gc, orient="index", columns=["gc"]
    )

    pd.concat([
        archaea_gc,
        bacteria_gc
    ]).to_csv("data/gc/genomic.csv", index_label="genome")
