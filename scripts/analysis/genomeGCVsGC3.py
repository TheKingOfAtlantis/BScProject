
import pandas as pd

genomeGC = pd.concat({
    "bacteria": pd.read_json("data/gc/bacteria.json", orient="table"),
    "archaea": pd.read_json("data/gc/archaea.json", orient="table"),
})
