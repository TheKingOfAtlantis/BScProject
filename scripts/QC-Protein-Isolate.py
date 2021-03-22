import pandas as pd
from Bio import SeqIO
import itertools, json

def isolate(df):
    ''' Isolates those which failed '''
    inv = ~df[["length", "start", "end", "internal"]] # Invert true <=> false
    return df[inv.any(axis=1)]

def getData(name):
    with open(f"data/qc/proteins/{name}.json") as file:
        data = json.load(file)
        data = pd.concat({
            k:pd.DataFrame(v)
            for k,v in data.items()
        })
    return isolate(data).rename_axis(index=["genome", "index"])

archaea  = getData("archaea")
bacteria = getData("bacteria")

pd.concat({
    "bacteria": bacteria,
    "arachaea": archaea
}).to_csv("data/qc/proteins/failed.csv")

def __findPseudo(x):
    genome, df, path = x
    isPseudo = pd.DataFrame()
    with open(f"{path}/{genome.split('.')[0]}.embl") as file:
        record = next(SeqIO.parse(file, "embl"))

        for feature in filter(lambda x: x.type == "CDS", record.features):
            # First we check that we have that ID
            row = None

            if("protein_id" in feature.qualifiers): row = df.loc[genome][df.loc[genome]["id"] == feature.qualifiers["protein_id"][0]]
            if("locus_tag" in feature.qualifiers):  row = df.loc[genome][df.loc[genome]["id"] == feature.qualifiers["locus_tag"][0]]

            if(row.empty): continue

            # This gene/feature is in our list so lets check it
            if("pseudo" in feature.qualifiers):
                row["genome"] = genome
                isPseudo = pd.concat([isPseudo, row])
    return isPseudo

def findPseudo(path, df):
    # Check which genes that we encountered are psuedo genes
    from multiprocessing import Pool
    import os, tqdm

    with Pool(os.cpu_count()) as pool:
        results = list(tqdm.tqdm(pool.imap(
            __findPseudo, zip(
                set(df.index.get_level_values(0)),
                itertools.repeat(df),
                itertools.repeat(path)
            )
        ), total=len(set(df.index.get_level_values(0)))))
    return pd.concat(results)

postPsuedo_archaea  = archaea[~archaea["id"].isin(findPseudo("data/genomes/archaea", archaea)["id"])]
postPsuedo_bacteria = bacteria[~bacteria["id"].isin(findPseudo("data/genomes/bacteria", bacteria)["id"])]

archaea.drop(["id_type", "id"], axis = 1, inplace = True)
bacteria.drop(["id_type", "id"], axis = 1, inplace = True)
postPsuedo_archaea.drop(["id_type", "id"], axis = 1, inplace = True)
postPsuedo_bacteria.drop(["id_type", "id"], axis = 1, inplace = True)

print("Fails w/ Psuedo")
print(pd.concat([pd.concat({
        "Archaea":  (~archaea).sum(),
        "Bacteria": (~bacteria).sum()
    }, axis = 1), pd.DataFrame({
        "Archaea":  [len(archaea)],
        "Bacteria": [len(bacteria)]
    }, index = ["total genes"])]
))

print("\nFails w/o Psuedo")
print(pd.concat([pd.concat({
        "Archaea":  (~postPsuedo_archaea).sum(),
        "Bacteria": (~postPsuedo_bacteria).sum()
    }, axis = 1), pd.DataFrame({
        "Archaea":  [len(postPsuedo_archaea)],
        "Bacteria": [len(postPsuedo_bacteria)]
    }, index = ["total genes"])]
))
