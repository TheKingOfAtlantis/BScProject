# This file fetches the genome EMBL files

import pathlib
import pandas as pd

from common import Download, Filesystem, Parallel

from IPython import get_ipython
if(get_ipython() != None):
    import os
    os.chdir("..")

filtered = pd.read_csv("data/genomes/build/filtered.csv")

buildDir    = "data/genomes/build/"
downloadDir = buildDir + "download/"
intermDir   = buildDir + "interm/"
outDir      = buildDir + "out/"

fileSuffix = "_genomic.gbff.gz"

domainMap = { 2: "bacteria",  2157: "archaea"}

files = pd.DataFrame()
files["accession"] = filtered.accession_refseq
files["domain"]    = filtered.superkingdom.map(domainMap)
files["url"]       = (
    "https://ftp.ncbi.nlm.nih.gov/genomes" + filtered.link + "/" + # Url to the assembly directory
    filtered.link.str.split("/").str[-1] + fileSuffix              # File to retrieve in the directory
)
files["gzFile"]    = downloadDir + files.accession +".gbff.gz"
files["gbFile"]    = intermDir   + files.accession + ".gbff"
files["emblFile"]  = outDir      + files.domain    + "/" + files.accession + ".embl"

Filesystem.mkdir([
    downloadDir,
    intermDir
] + [outDir + domain for domain in domainMap.values()])

# Now that we have all the files download as gz files
# We need to decompress them

def decompress(inputPath, isRetry = False):
    import gzip, shutil
    file = files[files.gzFile == inputPath].squeeze()
    try:
        with open(file.gbFile, "wb") as outputFile:
            with gzip.open(file.gzFile, "rb") as inputFile:
                shutil.copyfileobj(inputFile, outputFile)
    except Exception as e:
        if(not isRetry):
            # File probably failed to download properly
            # We'll give it another try
            Download.getFile(file.url, file.gzFile)
            decompress(inputFile, isRetry = True)
        else:
            # We've already tried to download the file again and still failed
            # One choice left
            print(accession)
            raise e

# To avoid having to completely rewrite our code for genbank files
# We should try to convert them to EMBL files

def convert(inputPath):
    import re
    from Bio import SeqIO

    file = files[files.gbFile == inputPath].squeeze()

    with open(file.emblFile, "w") as outputFile:
        with open(file.gbFile, "r") as inputFile:
            data = SeqIO.parse(inputFile, "genbank")
            SeqIO.write(data, outputFile, "embl")

# We can now download all the files
Download.getFile(
    files.url, files.gzFile,
    description = "Downloading Assemblies",
    disable=True, # Disables progress bar for individual files
)

Parallel.loadParallel(
    decompress, files.gzFile,
    len(files),
    desc = "Decompressing *.gbff.gz => *.gbff"
)

Parallel.loadParallel(
    convert, files.gbFile,
    len(files),
    desc = "Converting *.gbff => *.embl"
)
