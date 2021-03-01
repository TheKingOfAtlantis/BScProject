
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

if __name__ == "__main__":
    pathlib.Path("data/qc/trans_table/").mkdir(parents=True, exist_ok=True)
    pathlib.Path("plot/qc/trans_table/").mkdir(parents=True, exist_ok=True)

    archaea_tables = pd.DataFrame(dict(loadGlob("data/genomes/archaea/*", getTranslationTable))).T
    archaea_tables.to_csv("data/qc/trans_table/archaea.csv")
    archaea_tables.count().plot(kind='bar', figsize=(19,6)).get_figure().savefig("plot/qc/trans_table/archaea.png")
    archaea_tables.plot(kind='bar', figsize=(19, 6))

    fig, ax = plt.subplots(1, len(archaea_tables.columns), sharey=True)
    for i, table in enumerate(archaea_tables.columns):
        print(table)
        values = archaea_tables[table].sort_values()
        if(len(archaea_tables.columns) == 1): ax.bar(values.index, values)
        else: ax[i].bar(values.index, values)
        ax.set_xticklabels([])
    fig.savefig("plot/qc/trans_table/archaea-distribution.png")

    bacteria_tables = pd.DataFrame(dict(loadGlob("data/genomes/bacteria/*", getTranslationTable))).T
    bacteria_tables.to_csv("data/qc/trans_table/bacteria.csv")
    bacteria_tables.count().plot(kind='bar', figsize=(19,6)).get_figure().savefig("plot/qc/trans_table/bacteria.png")

    # Plots the number of CDS which use each translation table
    fig, ax = plt.subplots(1, len(bacteria_tables.columns), sharey=True)
    for i, table in enumerate(bacteria_tables.columns):
        print(table)
        values = bacteria_tables[table].sort_values()
        if(len(bacteria_tables.columns) == 1): ax.bar(values.index, values)
        else: ax[i].bar(values.index, values)
        ax.set_xticklabels([])
    fig.savefig("plot/qc/trans_table/bacteria_distribution.png")
