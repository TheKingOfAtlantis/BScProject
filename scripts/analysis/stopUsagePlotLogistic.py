import pandas as pd
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import statsmodels.api as sm

human = pd.read_json("data/gc/cds/human.json", orient = "table")
human[["gc", "gc1", "gc2", "gc3"]] = pd.DataFrame(human.gc.tolist(), index=human.index)
human = human[["shift", "stop", "gc3"]]

human = human.reset_index().rename(columns={"level_0": "id"}).drop("level_1", axis=1)

codons = ['TAA', 'TGA', 'TAG']
fig, axes = plt.subplots(nrows = len(codons), figsize = (16, 9), sharex=True)
for ax, stop in zip(axes, codons):
    data = pd.DataFrame({
        "x": human.gc3,
        "y": human.stop == stop # Whether each CDS has the stop we are checking for
    })

    print(f"Stop {stop}:")

    # Log logisitic regression curve
    fit = sm.Logit(data.y, sm.add_constant(data.x)).fit()
    predict_line = np.linspace(data.x.min(), data.x.max())
    ax.plot(predict_line, fit.predict(sm.add_constant(predict_line)), "-")

    ax.set_ylabel(f"{stop} Used")

    ax.set_yticks([0, 1])
    ax.set_yticklabels(["False", "True"])

    print(fit.summary2())

axes[-1].set_xlabel("Gene GC3 Content (%)")
plt.savefig(f"plot/gc/cds-human-logit.png")
