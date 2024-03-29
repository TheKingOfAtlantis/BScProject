
import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from common import Filesystem

import os
import pandas as pd
import numpy as np

import statsmodels.api as sm
import matplotlib.pyplot as plt

import seaborn as sns
sns.set(palette="colorblind")

genomicGC = pd.read_csv("data/gc/genomic.csv")

def loadDataFrame(file):
    df = pd.read_json(file, orient="table")

    df[["gc", "gc1", "gc2", "gc3"]] = pd.DataFrame(df.gc.tolist(), index=df.index)
    # For prokaryotes we want to pull the genomic GC
    # We don't use it for humans so no need to bother
    if(pathlib.Path(file.name).stem != "human"):
        df["id"] = df.index.get_level_values(0).values
        df = df.drop(columns="gc").merge(
            genomicGC,
            how="left",
            left_on="id",
            right_on="genome"
        ).set_index("id")
    df = df[["shift", "stop", "gc", "gc3"]]

    # We want to sort and bin the human genome
    if(pathlib.Path(file.name).stem == "human"):
        # Sort and bin
        bins = 100
        df = df.sort_values(by=["gc3"])
        df = pd.concat(
            np.array_split(df, bins),
            keys=[f"Human-${x}" for x in range(bins)]
        )

    return df.reset_index().rename(columns={"level_0": "id"}).drop("level_1", axis=1, errors='ignore')

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
        freq, values=["freq","gc","gc3"], index="id", columns=["stop", "shift"],
        aggfunc={"freq": np.sum, "gc": np.mean, "gc3": np.mean }
    )
    # Normalise the frequence by dividing by the number of codons in each frame found in each genome
    freq["freq"] = freq["freq"].divide(freq["freq"].groupby("shift", axis=1).sum(), axis=0)
    return freq.reorder_levels((2, 1, 0), axis=1).sort_index(axis=1) # Swap around order of columns

def plotGCvsFreq(name, freq, group, pdf = False):
    # Create new plot for each frameshift
    for gcType in ["gc", "gc3"]:
        for shift in freq.columns.get_level_values("shift").unique():
            figsize = (5 * 1.61803399, 5)
            fig, ax = plt.subplots(figsize=figsize)
            for codon in freq[shift].columns.get_level_values("stop").unique():
                # Get the frequence and GC(3)
                data = pd.DataFrame({
                    "x": freq[shift][codon][gcType],
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

            ax.set_xlabel(f"{gcType.upper()} Content (%)")
            ax.set_ylabel("Stop Codon Frequency (%)")

            ax.set_ylim((-5, 105))

            Filesystem.mkdir(f"plot/{gcType}/{group}")
            fig.savefig(f"plot/{gcType}/{group}/{name}-stop-shift{shift}.png")
            if(pdf): fig.savefig(f"plot/{gcType}/{group}/{name}-stop-shift{shift}.pdf")
            plt.close(fig)

def load(file, group):
    name  = os.path.basename(file.name).split(".")[0]
    df    = loadDataFrame(file)
    noTAG = df[df["stop"] != "TAG"].copy() # Drop TAG => TAA, TGA & TAC
    noTAC = df[df["stop"] != "TAC"].copy() # Drop TAG => TAA, TGA & TAG

    plotGCvsFreq(name, calculateFreq(df),    group)
    plotGCvsFreq(name, calculateFreq(noTAG), group + "+TAC")
    plotGCvsFreq(name, calculateFreq(noTAC), group + "+TAG")

if __name__ in "__main__":
    Filesystem.loadGlob("data/gc/cds/*.json",  load, group="cds",  desc="Plotting CDS Data")
    Filesystem.loadGlob("data/gc/trna/*.json", load, group="trna", desc="Plotting tRNA Data")
