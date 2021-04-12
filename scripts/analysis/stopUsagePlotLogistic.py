import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from common import Filesystem

import pandas as pd
import numpy as np

import statsmodels.api as sm
import matplotlib.pyplot as plt

human = pd.read_json("./data/gc/cds/human.json", orient = "table")
human[["gc", "gc1", "gc2", "gc3"]] = pd.DataFrame(human.gc.tolist(), index=human.index)
human = human[["shift", "stop", "gc", "gc3"]]

human = human.reset_index().rename(columns={"level_0": "id"}).drop("level_1", axis=1)

for stop in ['TAA', 'TAG', 'TGA']:
    human["freq"] = human["stop"] == stop
    result = sm.Logit(human["freq"], human["gc3"]).fit()

    print(f"Stop {stop}:")
    print(result.summary2())

    fig, ax = plt.subplots(figsize = (16, 9))
    ax.scatter(human['gc3'],human["freq"])
    predict_x = np.linspace(human['gc3'].min(), human['gc3'].max())
    ax.plot(predict_x,result.predict(predict_x), "-")

    Filesystem.mkdir("plot/gc/cds/")
    plt.savefig(f"plot/gc/cds/human-logit-{stop}.png")
    plt.close()
