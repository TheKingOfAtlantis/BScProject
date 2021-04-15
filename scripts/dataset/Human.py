
import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from IPython import get_ipython
if(get_ipython() != None):
    import os
    os.chdir("../..")

import gffutils
from common import Download, Parallel

# NCBI RefSeq has a lot of extra CDS annotated sequences which have very strange behaviours
# Switching over to GENCODE dataset, tried the comprehensive annotations first for a bit,
# but after a little bit of reading found that the basic annotation was more inline with what I
# was going to do with the comprehensive annotations anyway which was to filter for 1 CDS to

# url_human_seq = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"
url_human_seq = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz"
gz_human_seq  = "data/genomes/build/download/GRCh38_latest_genomic.fna.gz"
fna_human_seq = "data/genomes/build/out/GRCh38_latest_genomic.fna"

# url_human_annotation  = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz"
# url_human_annotation  = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gff3.gz"
url_human_annotation  = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.basic.annotation.gff3.gz"
gz_human_annotation   = "data/genomes/build/download/GRCh38_latest_genomic.gff3.gz"
gff3_human_annotation = "data/genomes/build/out/GRCh38_latest_genomic.gff3"
db_human_annotation   = "data/genomes/build/out/GRCh38_latest_genomic.gff3.sqlite"

Download.getFile(url_human_seq, gz_human_seq, desc = "Downloading Human Genome Sequence")
Download.getFile(url_human_annotation, gz_human_annotation, desc = "Downloading Human Genome Annotation")

def decompress(gzFile, out, **kwargs):
    import gzip, io, functools
    from tqdm.auto import tqdm
    with gzip.open(gzFile, 'rb') as srcFile, open(out, 'wb') as dstFile:
        # We determine the size of the file by seeking from the beginning => end
        size = srcFile.seek(0, io.SEEK_END)
        # Only so we can have a nice progress bar to see what is happening
        srcFile.seek(0)
        buffer_size = 512
        with tqdm(total=size, unit='B', unit_scale=True, **kwargs) as progress:
            for byte in iter(functools.partial(srcFile.read, buffer_size), b''):
                dstFile.write(byte)
                progress.update(len(byte))

decompress(gz_human_seq, fna_human_seq, desc = "Extracting Human Genome Sequence")
decompress(gz_human_annotation, gff3_human_annotation, desc = "Extracting Human Genome Annotation")

print("Creating Human Annotation Database: This takes a while")
annotations = gffutils.create_db(
    force   = True,
    verbose = True,
    data    = gff3_human_annotation,
    dbfn    = db_human_annotation,
    id_spec = "ID",
    merge_strategy = "create_unique",
    keep_order     = True,
    disable_infer_genes = True, disable_infer_transcripts = True
)
print("Creating Human Annotation Database: Done")

import pandas as pd

gffutils.constants.always_return_list = False

# **Rational**
# We want to pick the most representative CDS for each gene
# much like with the largest genome of each genus for prokaryotes
# to avoid genes with many CDSs potentially weighing the results down

# For reference: http://daler.github.io/gffutils/database-schema.html

from tqdm.auto import tqdm

# Get all the CDSs which are protein_coding

CDSs = [
    CDS for CDS in tqdm(
        annotations.features_of_type("CDS"),
        total = annotations.count_features_of_type("CDS")
    ) if CDS.attributes["transcript_type"] == "protein_coding"
]

# Then group the CDSs by gene
pairs = dict.fromkeys(set([ CDS.attributes["gene_id"] for CDS in CDSs ]), [])
for CDS in CDSs:
    pairs[CDS.attributes["gene_id"]].append(CDS)

# Then select the best cds from each gene

def selectCDS(x):
    gene, CDSs = x
    CDSs = sorted(CDSs, key = lambda x: x.end - x.start)
    return (gene, CDSs[0])

final = dict(Parallel.loadParallel(
    selectCDS, pairs.items(),
    len(pairs)
))
