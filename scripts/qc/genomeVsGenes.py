
import itertools, pathlib

from Bio import SeqIO
from scipy.spatial import distance
import pandas as pd
import numpy as np

def calculate(file):
    record = next(SeqIO.parse(file, "embl"))
    return pd.DataFrame({
        "size": len(record.seq),
        "count": len([x for x in record.features if x.type == "gene"])
    }, index=[record.id])

import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from common import loadGlob
if __name__ == "__main__":
    pathlib.Path("data/qc/count/").mkdir(parents=True, exist_ok=True)
    pathlib.Path("plot/qc/count/").mkdir(parents=True, exist_ok=True)

    archaea = pd.concat(loadGlob("data/genomes/archaea/*", calculate))
    archaea.to_csv("data/qc/count/archaea.csv")

    bacteria = pd.concat(loadGlob("data/genomes/bacteria/*", calculate))
    bacteria.to_csv("data/qc/count/bacteria.csv")

    # Use Mahalanobis Distance to find those that deviate significantly from the mean (3 sd)
    # Mahalanobis Distance generatlisation of Z-score (measure of the distance from mean)

    print("Running the R script via subprocess.call(...)")

    import subprocess
    subprocess.call("Rscript scripts/QC-GenomeVsGenes.R")
