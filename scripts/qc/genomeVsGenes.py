
import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from common import Filesystem

import pandas as pd

def calculate(file):
    from Bio import SeqIO
    record = next(SeqIO.parse(file, "embl"))
    return pd.DataFrame({
        "size": len(record.seq),
        "count": len([x for x in record.features if x.type == "gene"])
    }, index=[record.id])

if __name__ == "__main__":
    Filesystem.mkdir("data/qc/count/")
    Filesystem.mkdir("plot/qc/count/")

    archaea = pd.concat(Filesystem.loadGlob("data/genomes/archaea/*", calculate, desc = "Processing Archaea"))
    archaea.to_csv("data/qc/count/archaea.csv")

    bacteria = pd.concat(Filesystem.loadGlob("data/genomes/bacteria/*", calculate, desc = "Processing Bacteria"))
    bacteria.to_csv("data/qc/count/bacteria.csv")

    # Use Mahalanobis Distance to find those that deviate significantly from the mean (3 sd)
    # Mahalanobis Distance generatlisation of Z-score (measure of the distance from mean)

    print("Running the R script via subprocess.call(...)")

    import subprocess
    subprocess.call("Rscript scripts/qc/genomeVsGenes.R", shell=True)
