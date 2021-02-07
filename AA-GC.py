
# Calculate the GC content of each Amino Acid

import itertools

import pandas as pd
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.SeqUtils import GC

nt     = ["A", "T", "G", "C"]
codons = ["".join(x) for x in itertools.product(nt, repeat=3)]
codons = [Seq(x) for x in codons]

codons = pd.DataFrame({
    "aa": [x.translate() for x in codons],
    "codon": codons,
    "gc.total": [GC(x) for x in codons],
    "gc.third": [GC(x[2]) for x in codons]
}).set_index(["aa", "codon"]).sort_index()

AA_Names = pd.read_csv("aa_name.csv", index_col="single")
avg = codons.groupby(by="aa")[["gc.total", "gc.third"]].mean()
avg["count"] = codons.groupby(by="aa")["gc.total"].count()
avg["name"]  = AA_Names.loc[avg.index.map(lambda x: str(x))]["name"]

(fig, ax) = plt.subplots(figsize=(16, 9))
avg.plot(kind="bar", ax=ax, x="name", y=["gc.total", "gc.third"], label=["GC Total", "GC Third"])
ax.set_xlabel("Amino Acid")
ax.set_ylabel("%GC")

fig.savefig("plot/aa_gc.png")
fig.savefig("plot/aa_gc.pdf")


with pd.ExcelWriter('data/aa-gc.xlsx') as writer:
    codons.to_excel(writer, sheet_name="codon")
    avg.to_excel(writer, sheet_name="aa")
