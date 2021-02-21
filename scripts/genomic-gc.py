import os, io, itertools
from zipfile import ZipFile

import pandas as pd
from Bio import SeqIO, SeqUtils
from multiprocessing import Pool

def CalculateGC(x):
    zipPath, pos, total = x
    print("{}/{}".format(pos, total))
    with ZipFile(zipPath) as zipFile:
        path = zipFile.namelist()[pos]
        file = io.TextIOWrapper(zipFile.open(path))
        embl = list(SeqIO.parse(file, "embl"))
        assert len(embl) == 1
        record = embl[0]
        return (record.id, SeqUtils.GC(record.seq))

def GCFromZip(path):
    with ZipFile(path) as zipFile:
        count = len(zipFile.namelist())
    with Pool(os.cpu_count()) as pool:
        return pool.map(CalculateGC, zip(
            itertools.repeat(path),
            range(1, count),
            itertools.repeat(count)
        ))

import json
if __name__ == "__main__":
    archaea_gc  = GCFromZip("data/Archaea_filtered_genomes.zip")
    with open("data/archaea_gc.json", "w")  as js: json.dump(archaea_gc,  js)
    del archaea_gc

    start_time = time.time()
    bacteria_gc = GCFromZip("data/bacterial_filtered_genomes.zip")
    with open("data/bacteria_gc.json", "w") as js: json.dump(bacteria_gc, js)
    del bacteria_gc

os.chdir("data")
with ZipFile("gc.zip", "w") as out:
    out.write("archaea_gc.json")
    out.write("bacteria_gc.json")
    os.remove("archaea_gc.json")
    os.remove("bacteria_gc.json")
os.chdir("..")
