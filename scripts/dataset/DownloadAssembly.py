import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))

import pathlib
import pandas as pd

from common import Download, Filesystem, Parallel

from IPython import get_ipython
if(get_ipython() != None):
    import os
    os.chdir("..")

filtered = pd.read_csv("data/genomes/build/filtered.csv")

buildDir    = "data/genomes/build/"
downloadDir = buildDir  + "download/"
intermDir   = buildDir  + "interm/"
GbDir       = intermDir + "gb/"
EmblDir     = intermDir + "embl/"
outDir      = buildDir  + "out/"

# Maps NCBI taxanomic ID of each domain to a name
domainMap = { 2: "bacteria",  2157: "archaea"}

files = pd.DataFrame()
files["accession"] = filtered.accession_refseq
files["domain"]    = filtered.superkingdom.map(domainMap)
files["url"]       = (
    "https://ftp.ncbi.nlm.nih.gov/genomes" + filtered.link + "/" + # Url to the assembly directory
    filtered.link.str.split("/").str[-1] + "_genomic.gbff.gz"      # File to retrieve in the directory
)
files["gzFile"]    = downloadDir + files.accession + ".gbff.gz"
files["gbFile"]    = GbDir       + files.accession + ".gbff"
files["emblFile"]  = EmblDir     + files.accession + ".embl"

Filesystem.mkdir([
    downloadDir,
    GbDir,
    EmblDir
] + [outDir + domain for domain in domainMap.values()])

# Now that we have all the files download as gz files
# We need to decompress them

def decompress(inputPath, isRetry = False):
    import gzip, shutil
    file = files[files.gzFile == inputPath].squeeze()
    try:
        with gzip.open(file.gzFile, "rb") as srcFile,\
                open(file.gbFile, "wb") as dstFile:
                shutil.copyfileobj(srcFile, dstFile)
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
    with open(file.gbFile, "r") as srcFile,\
         open(file.emblFile, "w") as dstFile:
            SeqIO.convert(
                srcFile, "genbank",
                dstFile, "embl"
            )

# Finally to ensure we can still work with our preexisting code
# We need to rename the files to be the same as record.id

def rename(inputPath):
    from Bio import SeqIO
    import shutil
    file = files[files.emblFile == inputPath].squeeze()
    with open(file.emblFile, "r") as srcFile:
        seq     = next(SeqIO.parse(srcFile, "embl"))
        dstPath = outDir + f"{file.domain}/{seq.id}.embl"
        # Parsing with BioPython means the file is at EOF
        # so we need to reset the file to the start
        srcFile.seek(0, 0)
        with open(dstPath, "w") as dstFile:
            shutil.copyfileobj(srcFile, dstFile)
    return {
        "domain": file.domain,
        "accession_original": files.accession,
        "accession_renamed":  seq.id
    }

# We can now download all the files
Download.getFile(
    files.url, files.gzFile,
    description = "Downloading Assemblies",
    disable=True, # Disables progress bar for individual files
)

# Decompress the files
Parallel.loadParallel(
    decompress, files.gzFile,
    len(files),
    desc = "Decompressing *.gbff.gz => *.gbff"
)

# Convert to EMBL format
Parallel.loadParallel(
    convert, files.gbFile,
    len(files),
    desc = "Converting *.gbff => *.embl"
)

# Rename to sequence accession ID
map = Parallel.loadParallel(
    rename, files.emblFile,
    len(files),
    desc = "Renaming files to use Accession ID"
)

# This file is serves no functional purpose within the codebase (not even used to rebuild the same dataset)
# Instead provides a list, for those interested, of the genomes used and where they can be found
import pandas as pd
pd.DataFrame.from_dict(map).to_csv(outDir + "file_mapping.csv", index=False)
