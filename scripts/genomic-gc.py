import os
from zipfile import ZipFile
from common import loadZip
from Bio import SeqIO, SeqUtils

def CalculateGC(file):
    record = next(SeqIO.parse(file, "embl"))
    return (record.id, SeqUtils.GC(record.seq))

import json
if __name__ == "__main__":
    archaea_gc = loadZip("data/Archaea_filtered_genomes.zip", CalculateGC)
    with open("data/archaea_gc.json", "w") as js: json.dump(archaea_gc, js)
    del archaea_gc

    bacteria_gc = loadZip("data/bacterial_filtered_genomes.zip", CalculateGC)
    with open("data/bacteria_gc.json", "w") as js: json.dump(bacteria_gc, js)
    del bacteria_gc

    os.chdir("data/") # To avoid the zip file containing a "data/" subdirectory
    with ZipFile("gc.zip", "w") as out:
        out.write("archaea_gc.json")
        out.write("bacteria_gc.json")
        os.remove("archaea_gc.json")
        os.remove("bacteria_gc.json")
