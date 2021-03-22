
# QC

# Check translation tables: Some only use table 4 - Plot what is going on
# Check gene lengths are multiple of three

# Plot genomes size vs no of genes
# Should fit nicely to a straight line

import itertools, pathlib
from Bio import SeqIO
from collections import Counter

def getTranslationTable(file):
    record = next(SeqIO.parse(file, "embl"))

    # CDS contain translation table
    # To be completionist will pull out table for all CDS features
    return (
        record.id,
        Counter(itertools.chain.from_iterable(
            x.qualifiers["transl_table"] for x in record.features if x.type == "CDS"
        ))
    )

from common import loadGlob

import pandas as pd
import matplotlib.pyplot as plt

def plotDistribution(data):
    fig, ax = plt.subplots(1, len(data.columns), sharey=True)
    for i, table in enumerate(data.columns):
        values = data[table].sort_values()
        if(len(data.columns) == 1):
            ax.bar(values.index, values)
            ax.set_xticklabels([])
            ax.set_xlabel("Genomes")
            ax.set_ylabel("No. of CDS features")
        else:
            ax[i].bar(values.index, values)
            ax[i].set_xticklabels([])
            ax[i].set_xlabel("Genomes")
            ax[i].set_ylabel("No. of CDS features")
    return fig

if __name__ == "__main__":
    pathlib.Path("data/qc/trans_table/").mkdir(parents=True, exist_ok=True)
    pathlib.Path("plot/qc/trans_table/").mkdir(parents=True, exist_ok=True)

    archaea_tables = pd.DataFrame(dict(loadGlob("data/genomes/archaea/*", getTranslationTable))).T
    archaea_tables.to_csv("data/qc/trans_table/archaea.csv")
    archaea_tables.count().plot(kind='bar', figsize=(6,6)).get_figure().savefig("plot/qc/trans_table/archaea.png")
    plotDistribution(archaea_tables).savefig("plot/qc/trans_table/archaea-distribution.png")

    bacteria_tables = pd.DataFrame(dict(loadGlob("data/genomes/bacteria/*", getTranslationTable))).T
    bacteria_tables.to_csv("data/qc/trans_table/bacteria.csv")
    bacteria_tables.count().plot(kind='bar', figsize=(19,6)).get_figure().savefig("plot/qc/trans_table/bacteria.png")
    plotDistribution(bacteria_tables).savefig("plot/qc/trans_table/bacteria-distribution.png")

    pd.concat({
        "Bacteria": bacteria_tables.count() / bacteria_tables.count().sum(),
        "Archaea":  archaea_tables.count() / archaea_tables.count().sum()
    }, axis = 1).plot(
        kind="bar",
        xlabel="Translation Table", ylabel="Frequency",
        rot=0, logy=True
    ).get_figure().savefig("plot/qc/trans_table/both.png", dpi=600)
