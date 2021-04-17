import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))

from IPython import get_ipython
from IPython.display import display
if(get_ipython() != None):
    import os
    os.chdir("../..")

import pandas as pd
from Bio import SeqIO
import itertools, json

from common import Parallel, getID

def isolate(df):
    ''' Isolates those which failed '''
    inv = ~df # Invert true <=> false
    return df[inv.any(axis=1)]

def getData(name):
    print(f"Isolating non-conforming {name} CDSs")
    data = pd.read_json(f"data/qc/proteins/{name}.json", orient = "table")
    return isolate(data)

archaea  = getData("archaea")
bacteria = getData("bacteria")

isolated = pd.concat({
    "bacteria": bacteria,
    "archaea": archaea
}, names = ["domain"])
isolated.to_csv("data/qc/proteins/failed.csv")

def __findPseudo(x):
    (domain, genome), df = x
    ids = df.index.get_level_values("id")

    isPseudo = pd.DataFrame()
    with open(f"data/genomes/{domain}/{genome}.embl") as file:
        record   = next(SeqIO.parse(file, "embl"))
        CDSs     = filter(lambda feature: feature.type == "CDS", record.features) # Filter down to CDSs (as our IDs work here)
        isolated = filter(lambda cds: getID(cds) in ids, CDSs)                    # Isolate thoses whose IDs we have as isolated
        psuedo   = filter(lambda cds: "pseudo" in cds.qualifiers, isolated)       # Filter for thoses with the pseudo annotation#

        psuedoIds = [getID(feature) for feature in psuedo]
        return df.loc[(domain, genome, psuedoIds)]

def findPseudo(df):
    # Check which genes that we encountered are psuedo genes
    return pd.concat(
        Parallel.loadParallel(
            __findPseudo, df.groupby(["domain", "genome"]),
            count = len(df.index.get_level_values("genome").unique()),
            desc = "Isolating pseudogenes"
        )
    )

psuedo = findPseudo(isolated)
psuedo = isolated.loc[psuedo.index.values]

isolated["psuedo"] = isolated.isin(psuedo).any(axis=1)

summary = pd.concat({
    "w/ Psuedo": isolated.groupby("domain").sum(),
    "w/o Psuedo": isolated[isolated.psuedo == False].groupby("domain").sum()
}, names=["group"])
summary["total genes"] = pd.concat({
     "w/ Psuedo": isolated.groupby("domain").size(),
     "w/o Psuedo": isolated[isolated.psuedo == False].groupby("domain").size()
})
summary = summary.drop(columns="psuedo").T
summary.to_csv("data/qc/proteins/failed-summary.csv")

print("Non-conforming Summary")
display(summary)
