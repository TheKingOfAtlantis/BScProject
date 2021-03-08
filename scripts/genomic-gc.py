from common import loadGlob
from Bio import SeqIO, SeqUtils

def CalculateGC(file):
    record = next(SeqIO.parse(file, "embl"))
    return (record.id, SeqUtils.GC(record.seq))

import json, pathlib
if __name__ == "__main__":
    pathlib.Path(f"data/gc/").mkdir(parents=True, exist_ok=True)

    archaea_gc = loadGlob("data/genomes/archaea/*", CalculateGC)
    with open("data/gc/archaea.json", "w") as js: json.dump(archaea_gc, js)
    del archaea_gc

    bacteria_gc = loadGlob("data/genomes/bacteria/*", CalculateGC)
    with open("data/gc/bacteria.json", "w") as js: json.dump(bacteria_gc, js)
    del bacteria_gc
