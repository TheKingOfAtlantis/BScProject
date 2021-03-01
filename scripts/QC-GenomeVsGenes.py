
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

    from rpy2 import robjects
    from rpy2.robjects import pandas2ri
    robjects.pandas2ri.activate()
    robjects.r.source('scripts/QC-GenomeVsGenes.R')
    archaea  = robjects.conversion.rpy2py(robjects.r["archaea"])
    bacteria = robjects.conversion.rpy2py(robjects.r["bacteria"])

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(16, 9))
    archaea.plot(kind="scatter", x="size", y="count", ax=ax)
    archaea[archaea["dist"] > 3].plot(kind="scatter", x="size", y="count", ax=ax, color="orange", label="outliers")
    archaea[archaea["dist"] > 3].apply(lambda x: ax.annotate(
        text   = x.name,
        xy     = x[["size", "count"]],
        xytext = (5, -1),
        textcoords="offset points"
    ), axis=1)
    fig.savefig("plot/qc/count/archaea.png")

    fig, ax = plt.subplots(figsize=(16, 9))
    bacteria.plot(kind="scatter", x="size", y="count", ax=ax)
    bacteria[bacteria["dist"] > 3].plot(kind="scatter", x="size", y="count", ax=ax, color="orange", label="outliers")
    bacteria[bacteria["dist"] > 3].apply(lambda x: ax.annotate(
        text   = x.name,
        xy     = x[["size", "count"]],
        xytext = (5, -4),
        textcoords="offset points"
    ), axis=1)
    fig.savefig("plot/qc/count/bacteria.png")
