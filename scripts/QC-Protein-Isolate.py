import pandas as pd
from Bio import SeqIO

def isolate(df):
    ''' Isolates those which failed '''
    inv = ~df # Invert true <=> false
    return df[inv.any(axis=1)]

if __name__ == "__main__":
    import itertools, json

    with open("data/qc/proteins/archaea.json") as file:
        archaea = json.load(file)
        archaea = pd.concat({
            k:pd.DataFrame(v).drop("position", axis = 1)
            for k,v in archaea.items()
        })
        archaea = isolate(archaea)

    with open("data/qc/proteins/bacteria.json") as file:
        bacteria = json.load(file)
        bacteria = pd.concat({
            k:pd.DataFrame(v).drop("position", axis = 1)
            for k,v in bacteria.items()
        })
        bacteria = isolate(bacteria)

    pd.concat({
        "bacteria": bacteria,
        "arachaea": archaea
    }).to_csv("data/qc/proteins/failed.csv")
