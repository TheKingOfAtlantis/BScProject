
import itertools, pathlib

from Bio import SeqIO
from scipy.spatial import distance
import pandas as pd
import numpy as np

def removePseudo(features):
    return list(filter(lambda x: not "pseudo" in x.qualifiers, features))

def calculate(file):
    record = next(SeqIO.parse(file, "embl"))
    return pd.DataFrame({
        "size": len(record.seq),
        "count": len(removePseudo([x for x in record.features if x.type == "CDS"]))
    }, index=[record.id])

from common import loadGlob
if __name__ == "__main__":
    pathlib.Path("data/qc/count/").mkdir(parents=True, exist_ok=True)
    pathlib.Path("plot/qc/count/").mkdir(parents=True, exist_ok=True)

    archaea = pd.concat(loadGlob("data/genomes/archaea/*", calculate))
    archaea.to_csv("data/qc/count/archaea.csv")

    bacteria = pd.concat(loadGlob("data/genomes/bacteria/*", calculate))
    bacteria.to_csv("data/qc/count/bacteria.csv")

    import os, subprocess
    subprocess.call("Rscript --vanilla --slave scripts/QC-GenomeVsGenes.R", shell=True)
    if(os.path.exists("Rplots.pdf")):
        os.remove("Rplots.pdf")
