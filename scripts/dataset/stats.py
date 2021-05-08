
import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))

from common import Filesystem, Parallel

import os
import pandas as pd
import numpy as np

import statsmodels.api as sm

import matplotlib.ticker as mticker
import matplotlib.pyplot as plt

import seaborn as sns
sns.set(palette="colorblind")

def getDB():
    from taxadb.taxid import TaxID
    return TaxID(dbtype='sqlite', dbname='data/genomes/build/taxadb.sqlite')

def getLineage(taxid):
    taxaDB = getDB()
    out = taxaDB.lineage_id(taxid, ranks = True)
    if(out is None): return { "missing": taxid }
    else: return dict(out)

def lineageConcat(block):
    return {x[0]:pd.Series(x[1], dtype='int') for x in block}

Filesystem.mkdir("plot/stats/")

prefiltered = pd.read_csv("data/genomes/build/ncbi.csv")

# Reusing the dataset/Filter.py code to get genus

with Parallel.getPool() as pool:
    # Some taxa ids referred to by multiple accession no.
    # So we filter our list to a unique set of taxa IDs

    # Now we need to group by genus and select just one species to use
    # First lets get the lineage of each bacteria taxa ID
    taxa_result = pd.DataFrame({ "tax_id": pd.unique(prefiltered.tax_id) })
    taxa_result["lineage"] = Parallel.loadParallel(
        getLineage,
        taxa_result.tax_id,
        count = len(taxa_result),
        desc="Retrieving lineages"
    )

taxa = Parallel.concat(
    taxa_result.values,
    lineageConcat, axis = 1
).T.rename_axis("tax_id").reset_index()

prefiltered = pd.merge(
    prefiltered, taxa[taxa.genus.isnull() == False],
    on       = "tax_id",
    validate = "m:1" # Ensure its a many-to-one relation
)

# To get number of species per genus
genusStats = prefiltered.groupby("genus").size().sort_values(ascending=False)

figsize = (5 * 1.61803399, 5)

fig, ax = plt.subplots(figsize=figsize)
ax.scatter(genusStats.index.astype("str"), genusStats.values)
ax.set_xticklabels([])
ax.set_xticks([])
ax.set_xlabel("Genus")
ax.set_ylabel("No. of Genomes in Genus")
ax.get_figure().savefig("plot/stats/prefiltered_genus_distrib.png", dpi=400)

fig, ax = plt.subplots(figsize=figsize)
ax = genusStats.plot.hist(ax=ax, logy=True, bins=20)
ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
ax.set_xlabel("Frequency of Genomes in Each Genus")
ax.set_ylabel("Number of Genera")
ax.get_figure().savefig("plot/stats/prefiltered_genus_histogram.png")
