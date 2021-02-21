
import json, itertools
import pandas as pd
import numpy as np
from scipy import stats

import matplotlib.pyplot as plt

with open("data/archaea_gene_gc.json", 'r') as file: archaea_gene_gc = json.load(file)

df = pd.concat({k:pd.DataFrame(v) for k,v in archaea_gene_gc.items()})

df[["gc", "gc1", "gc2", "gc3"]] = pd.DataFrame(df.gc.tolist(), index=df.index)
df = df[["shift", "stop", "gc", "gc3"]]

df = df.reset_index().rename(columns={"level_0": "id"}).drop("level_1", axis=1)
freq = df
freq["freq"] = 1
freq = pd.pivot_table(
    df, values=["freq", "gc", "gc3"], index="id", columns=["stop", "shift"],
    aggfunc={"freq": np.sum, "gc": np.mean, "gc3": np.mean }
)
freq["freq"] = freq["freq"].divide(freq["freq"].groupby("shift").sum(axis=1), axis=0)
freq = freq.reorder_levels((2, 1, 0), axis=1).sort_index(axis=1)


for shift in freq.columns.get_level_values("shift").unique():
    fig, ax = plt.subplots(figsize=(16,9))
    for codon in freq[shift].columns.get_level_values("stop").unique():
        xData = freq[shift][codon]["gc"].dropna()
        yData = freq[shift][codon]["freq"].dropna()

        colour = next(ax._get_lines.prop_cycler)['color']

        linreg = stats.linregress(xData, yData)
        sc = ax.scatter(xData, yData, label=codon, color=colour)

        x = np.linspace(xData.min(), xData.max())
        ax.plot(x, x * linreg.slope + linreg.intercept, color=colour)
        ax.annotate(
            f"y = {linreg.slope:.2g}x + {linreg.intercept:.2g}\n$r^2=${linreg.rvalue:.2g},p-value={linreg.pvalue:.2g}\nError: m={linreg.stderr:.2g}, c={linreg.intercept_stderr:.2g}",
            xy=(xData.max(), xData.max() * linreg.slope + linreg.intercept)
        )
    ax.legend()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xlabel("%GC3 Content")
    ax.set_ylabel("Relative Stop Frequency")

    fig.savefig(f"plot/archaea-stop-gc-shift{shift}.png")
    fig.savefig(f"plot/archaea-stop-gc-shift{shift}.pdf")
