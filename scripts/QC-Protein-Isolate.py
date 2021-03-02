import pandas as pd
from Bio import SeqIO
import itertools, json

def isolate(df):
    ''' Isolates those which failed '''
    inv = ~df # Invert true <=> false
    return df[inv.any(axis=1)]

def getData(name):
    with open(f"data/qc/proteins/{name}.json") as file:
        data = json.load(file)
        data = pd.concat({
            k:pd.DataFrame(v).drop("position", axis = 1)
            for k,v in data.items()
        })
        return isolate(data)

archaea  = getData("archaea")
bacteria = getData("bacteria")

pd.concat({
    "bacteria": bacteria,
    "arachaea": archaea
}).to_csv("data/qc/proteins/failed.csv")

def findPseudo(path, df):
    # Check which genes that we encountered are psuedo genes
    isPseudo = []
    for genome in set(df.index.get_level_values(0)):
        with open(f"{path}/{genome.split('.')[0]}.embl") as file:
            record   = next(SeqIO.parse(file, "embl"))
            features = list(filter(lambda x: x.type in "CDS", record.features))
            for pos in df.loc[genome].index:
                if("pseudo" in features[pos].qualifiers):
                    isPseudo.append((genome, pos))
    return isPseudo

print("Fails w/ Psuedo")
print(pd.concat({
    "Archaea":  (~archaea).sum(),
    "Bacteria": (~bacteria).sum()
}, axis = 1))
print("Fails w/o Psuedo")
print(pd.concat({
    "Archaea":  (~archaea.drop(findPseudo("data/genomes/archaea", archaea))).sum(),
    "Bacteria": (~bacteria.drop(findPseudo("data/genomes/bacteria", bacteria))).sum()
}, axis = 1))
