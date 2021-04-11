
import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from IPython import get_ipython
if(get_ipython() != None):
    import os
    os.chdir("../..")

import gffutils
from common import Download

url_human_seq = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"
gz_human_seq  = "data/genomes/build/download/GRCh38_latest_genomic.fna.gz"
fna_human_seq = "data/genomes/build/out/GRCh38_latest_genomic.fna"

url_human_annotation = "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz"
gz_human_annotation  = "data/genomes/build/download/GRCh38_latest_genomic.gff3.gz"
fna_human_annotation = "data/genomes/build/out/GRCh38_latest_genomic.gff3"
db_human_annotation  = "data/genomes/build/out/GRCh38_latest_genomic.gff3.sqlite"

Download.getFile(url_human_seq, gz_human_seq, desc = "Downloading Human Genome Sequence")
Download.getFile(url_human_annotation, gz_human_annotation, desc = "Downloading Human Genome Annotation")

def decompress(gzFile, out, **kwargs):
    import gzip, io, functools
    from tqdm.auto import tqdm
    with gzip.open(gzFile, 'rb') as srcFile, open(out, 'wb') as dstFile:
        # We determine the size of the file by seeking from the beginning => end
        size = srcFile.seek(0, io.SEEK_END)

        # Only so we can have a nice progress bar to see what is happening
        buffer_size = 512
        with tqdm(total=size, unit='B', unit_scale=True, **kwargs) as progress:
            for byte in tqdm(iter(functools.partial(srcFile.read, buffer_size), b'')):
                dstFile.write(byte)
                progress.update(len(byte))

decompress(gz_human_seq, fna_human_seq, desc = "Extracting Human Genome Sequence")
decompress(gz_human_annotation, fna_human_annotation, desc = "Extracting Human Genome Annotation")

print("Creating Human Annotation Database: Building... (no progress bar avaliable)", end = "")
gffutils.create_db(
    force = True,
    data  = 'data/genomes/GRCh38_latest_genomic.gff3',
    dbfn  = 'data/genomes/GRCh38_latest_genomic.gff3.sqlite',
    id_spec="ID",
    merge_strategy="create_unique",
    keep_order=True,
    disable_infer_genes=True, disable_infer_transcripts=True
)
print("\rCreating Human Annotation Database: Done")
