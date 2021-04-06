# This file fetches the genome EMBL files

import pathlib
import pandas as pd

from common import Download

from IPython import get_ipython
if(get_ipython() != None):
    import os
    os.chdir("..")

filtered = pd.read_csv("data/genomes/build/filtered.csv")

fileSuffix = "_genomic.gbff.gz"
files = filtered.link.str.split("/").str[-1] + fileSuffix

outputDir   = "data/genomes/build/download/"
outputFiles = outputDir + files
pathlib.Path(outputDir).mkdir(parents=True, exist_ok=True)

urls = "https://ftp.ncbi.nlm.nih.gov/genomes" + filtered.link + "/"
urls = urls + files

Download.getFile(urls, outputFiles)
