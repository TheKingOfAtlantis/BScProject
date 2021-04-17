import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))

from common import Filesystem
from Bio import SeqIO, SeqUtils

def CalculateGC(file):
    record = next(SeqIO.parse(file, "embl"))
    return (record.id, SeqUtils.GC(record.seq))

import json, pathlib
if __name__ == "__main__":
    Filesystem.mkdir("data/gc/")

    archaea_gc = Filesystem.loadGlob("data/genomes/archaea/*", CalculateGC, desc = "Archaea Genome GC")
    with open("data/gc/archaea.json", "w") as js: json.dump(archaea_gc, js)

    bacteria_gc = Filesystem.loadGlob("data/genomes/bacteria/*", CalculateGC, desc = "Bacteria Genome GC")
    with open("data/gc/bacteria.json", "w") as js: json.dump(bacteria_gc, js)
