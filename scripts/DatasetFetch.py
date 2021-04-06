# This file fetches the genome EMBL files

import pathlib
import pandas as pd

from common import Download, Filesystem, Parallel

from IPython import get_ipython
if(get_ipython() != None):
    import os
    os.chdir("..")

filtered = pd.read_csv("data/genomes/build/filtered.csv")

fileSuffix = "_genomic.gbff.gz"
files      = filtered.link.str.split("/").str[-1] + fileSuffix

outputDirRoot = "data/genomes/build/"
outputDir     = outputDirRoot + "download/"
outputFiles   = outputDir + files
Filesystem.mkdir(outputDir)

urls = "https://ftp.ncbi.nlm.nih.gov/genomes" + filtered.link + "/"
urls = urls + files


# We can now download all the files
# FIXME: Progress bar output is a mess

Download.getFile(urls, outputFiles)

# Now that we have all the files download as gz files
# We need to decompress them

filtered.superkingdom = filtered.superkingdom.map({ 2: "bacteria",  2157: "archaea"})

Filesystem.mkdir(outputDirRoot + "interm/bacteria")
Filesystem.mkdir(outputDirRoot + "interm/archaea")

def decompress(inputPath):
    import gzip, shutil, re
    # Remove the .gz extension
    # Change output directory from download/ to out/
    accession  = re.findall("(GCF_\d+\.\d+)", inputPath)[0]
    info       = filtered[filtered.accession_refseq == accession].squeeze()
    outputPath = inputPath[:-len(".gz")].replace("download", f"interm/{info.superkingdom}")

    with open(outputPath, "wb") as outputFile:
        with gzip.open(inputPath, "rb") as inputFile:
            shutil.copyfileobj(inputFile, outputFile)

Parallel.loadParallel(
    decompress,
    outputFiles,
    len(outputFiles),
    desc = "Decompressing *.gbff.gz => *.gbff"
)

# To avoid having to completely rewrite our code for genbank files
# We should try to convert them to EMBL files

def convert(x):
    import re
    from Bio import SeqIO

    inputPath, domain = x

    accession  = re.findall("(GCF_\d+\.\d+)", inputPath)[0]
    outputPath = outputDirRoot + f"out/{domain}/" + accession + ".embl"

    with open(outputPath, "w") as outputFile:
        with open(inputPath, "r") as inputFile:
            data = SeqIO.parse(inputFile, "genbank")
            SeqIO.write(data, outputFile, "embl")

import itertools, glob
for domain in ["archaea", "bacteria"]:
    Filesystem.mkdir(outputDirRoot + f"out/{domain}/")

    files = glob.glob(outputDirRoot + f"interm/{domain}/*")
    Parallel.loadParallel(
        convert, zip(files, itertools.repeat(domain)),
        len(files),
        desc = f"Converting {domain}/*.gbff => *.embl"
    )
