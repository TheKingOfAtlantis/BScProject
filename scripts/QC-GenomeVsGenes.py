
import itertools, pathlib

from Bio import SeqIO
import pandas as pd

def calculate(file):
    record = next(SeqIO.parse(file, "embl"))
    return pd.DataFrame({
        "size": len(record.seq),
        "count": len(record.features)
    }, index=[record.id])

from common import loadGlob
if __name__ == "__main__":
    pathlib.Path("data/qc/count/").mkdir(parents=True, exist_ok=True)
    pathlib.Path("plot/qc/count/").mkdir(parents=True, exist_ok=True)

    archaea = pd.concat(loadGlob("data/genomes/archaea/*", calculate))
    archaea.plot(kind="scatter", x="size", y="count", figsize=(19,6)).get_figure().savefig("plot/qc/count/archaea.png")
    archaea.to_csv("data/qc/count/archaea.csv")

    bacteria = pd.concat(loadGlob("data/genomes/bacteria/*", calculate))
    bacteria.plot(kind="scatter", x="size", y="count", figsize=(19,6)).get_figure().savefig("plot/qc/count/bacteria.png")
    bacteria.to_csv("data/qc/count/bacteria.csv")
