import os, json, itertools, pathlib
import pandas as pd
import numpy as np

import statsmodels.api as sm

import matplotlib.pyplot as plt

def loadDataFrame(file):
    df = pd.read_json(file, orient="table")

    df[["gc", "gc1", "gc2", "gc3"]] = pd.DataFrame(df.gc.tolist(), index=df.index)
    df = df[["shift", "stop", "gc", "gc3"]]

    # We want to sort and bin the human genome
    if(pathlib.Path(file.name).stem == "human"):
        # Sort and bin
        df = df.sort_values(by=["gc3"])
        df = pd.concat(
            np.array_split(df, 100),
            keys=[f"Human-${x}" for x in range(100)]
        )


    return df.reset_index().rename(columns={"level_0": "id"}).drop("level_1", axis=1)

def calculateFreq(df):
    freq = df.copy()
    freq["freq"] = 1 # As of right now we can see the frequence of each is 1 (as this is true)

    # We then use pivot to perform the bulk of the magic
    # We use the IDs as the new indices (which is why we called reset_index earlier)
    # We create columns for each stop codon and frame-shift value
    # The the value of these columns correspond to the:
    #  - Total frequency of each codon in each frame
    #  - Average GC of each sequence in each frame
    #  - Average GC3 of each sequence in each frame
    freq = pd.pivot_table(
        freq, values=["freq","gc3"], index="id", columns=["stop", "shift"],
        aggfunc={"freq": np.sum, "gc3": np.mean }
    )
    # Normalise the frequence by dividing by the number of codons in each frame found in each genome
    freq["freq"] = freq["freq"].divide(freq["freq"].groupby("shift", axis=1).sum(), axis=0)
    return freq.reorder_levels((2, 1, 0), axis=1).sort_index(axis=1) # Swap around order of columns

def plotGCvsFreq(name, freq, group, pdf = False):
    # Create new plot for each frameshift
    for shift in freq.columns.get_level_values("shift").unique():
        fig, ax = plt.subplots(figsize=(16,9))
        for codon in freq[shift].columns.get_level_values("stop").unique():
            # Get the frequence and GC(3)
            data = pd.DataFrame({
                "x": freq[shift][codon]["gc3"],
                "y": freq[shift][codon]["freq"] * 100
            })

            data.dropna(inplace = True)

            # Use this so we can give both line and scatter points the same colour
            colour = next(ax._get_lines.prop_cycler)['color']

            # Plot the raw data values
            ax.scatter(data.x, data.y, label=codon, color=colour)

            fit = sm.OLS(data.y, sm.add_constant(data.x)).fit()
            ax.plot(data.x, fit.predict(sm.add_constant(data.x)), color=colour)

        ax.legend()

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.set_xlabel("GC3 Content (%)")
        ax.set_ylabel("Stop Codon Frequency (%)")

        ax.set_ylim((-5, 105))

        pathlib.Path(f"plot/gc/{group}").mkdir(parents=True, exist_ok=True)
        fig.savefig(f"plot/gc/{group}/{name}-stop-shift{int(shift)}.png")
        if(pdf): fig.savefig(f"plot/gc/{group}/{name}-stop-shift{int(shift)}.pdf")


import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from common import loadGlob

def load(file, group):
    name  = os.path.basename(file.name).split(".")[0]
    df    = loadDataFrame(file)
    noTAG = df[df["stop"] != "TAG"].copy() # Drop TAG => TAA, TGA & TAC
    noTAC = df[df["stop"] != "TAC"].copy() # Drop TAG => TAA, TGA & TAG

    plotGCvsFreq(name, calculateFreq(noTAG), group + "+TAC")
    plotGCvsFreq(name, calculateFreq(noTAC), group + "+TAG")

if __name__ in "__main__":
    loadGlob("data/gc/cds/*.json", load, extra="cds")
    loadGlob("data/gc/trna/*.json", load, extra="trna")
