
import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from IPython import get_ipython
if(get_ipython() != None):
    import os
    os.chdir("../..")

from common import Download, Parallel
import gffutils
import pandas as pd
from tqdm.auto import tqdm


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

decompress(gz_human_seq, fna_human_seq, desc = "Decompressing Human Genome Sequence")
decompress(gz_human_annotation, gff3_human_annotation, desc = "Decompressing Human Genome Annotation")

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

print("Creating Human Annotation Database: Annotation Tables Extension")

gffutils.constants.always_return_list = False

# **Rational**
# We want to pick the most representative CDS for each gene
# much like with the largest genome of each genus for prokaryotes
# to avoid genes with many CDSs potentially weighing the results down

# For reference: http://daler.github.io/gffutils/database-schema.html

# To enable queries totally within SQLite (for faster execution) we need
# to extract the avaliable attributes and embedded them within the SQLite database
# This doesn't effect the existing information produced by gffutils, as these
# stored in seperate tables and the original tables are left untouched

def isolateAttribute(feature):
    def extractAttribute(attributes, name):
        return attributes[name] if name in attributes else None
    return {
       "ID":                       feature.id,
       "Parent":                   extractAttribute(feature.attributes, "Parent"),
       "ccdsid":                   extractAttribute(feature.attributes, "ccds"),
       "gene_id":                  extractAttribute(feature.attributes, "gene_id"),
       "gene_type":                extractAttribute(feature.attributes, "gene_type"),
       "gene_name":                extractAttribute(feature.attributes, "gene_name"),
       "level":                    extractAttribute(feature.attributes, "level"),
       "hgnc_id":                  extractAttribute(feature.attributes, "hgnc_id"),
       "havana_gene":              extractAttribute(feature.attributes, "havana_gene"),
       "transcript_id":            extractAttribute(feature.attributes, "transcript_id"),
       "transcript_type":          extractAttribute(feature.attributes, "transcript_type"),
       "transcript_name":          extractAttribute(feature.attributes, "transcript_name"),
       "transcript_support_level": extractAttribute(feature.attributes, "transcript_support_level"),
       "havana_transcript":        extractAttribute(feature.attributes, "havana_transcript"),
       "exon_number":              extractAttribute(feature.attributes, "exon_number"),
       "exon_id":                  extractAttribute(feature.attributes, "exon_id")
    }

def isolateLists(feature, name):
    return [{
        "featureId": feature.id,
        "ordering": i,
        "value": val,
    } for i,val in enumerate(feature.attributes[name])]


with annotations.conn as connection:
    # Create the attribute tables
    def listTable(name): return f"""
        CREATE TABLE IF NOT EXISTS {name} (
            featureId TEXT,
            ordering INT,
            value TEXT,
            PRIMARY KEY (featureId, ordering),
            FOREIGN KEY (featureId) REFERENCES features(id)
        );
        """
    connection.executescript(
        f"""
        CREATE TABLE IF NOT EXISTS attributes (
            id TEXT,
            ccdsid TEXT,
            gene_id TEXT,
            gene_type TEXT,
            gene_name TEXT,
            level TEXT,
            hgnc_id TEXT,
            havana_gene TEXT,
            Parent TEXT,
            transcript_id TEXT,
            transcript_type TEXT,
            transcript_name TEXT,
            transcript_support_level TEXT,
            havana_transcript TEXT,
            exon_number TEXT,
            exon_id TEXT,
            PRIMARY KEY (id),
            FOREIGN KEY (id) REFERENCES features(id)
        );
        {listTable("tags")}
        {listTable("ont")}
        """
    )

    attributesList = []
    tagsList = []
    ontList = []

    # Extract the attributes
    for feature in tqdm(
        annotations.all_features(),
        total = annotations.count_features_of_type()
        desc  = "Extracting attribute information"
    ):
        attributesList.append(isolateAttribute(feature))
        if("tag" in feature.attributes): tagsList.extend(isolateLists(feature, "tag"))
        if("ont" in feature.attributes): ontList.extend(isolateLists(feature, "ont"))

    print("Creating attribute tables")
    # Insert the attributes
    connection.executemany(
        """
        INSERT INTO attributes VALUES (
            :ID,
            :ccdsid,
            :gene_id,
            :gene_type,
            :gene_name,
            :level,
            :hgnc_id,
            :havana_gene,
            :Parent,
            :transcript_id,
            :transcript_type,
            :transcript_name,
            :transcript_support_level,
            :havana_transcript,
            :exon_number,
            :exon_id
        );
        """, attributesList
    )
    connection.executemany(
        """
        INSERT INTO tags VALUES (
            :featureId,
            :ordering,
            :value
        );
        """, tagsList
    )
    connection.executemany(
        """
        INSERT INTO ont VALUES (
            :featureId,
            :ordering,
            :value
        );
        """, ontList
    )

print("Creating Human Annotation Database: Done")

# Brief explaination of this SQL query:
# ------------------------------------
# We want to find the ID of the longest CDS in each gene to store it for
# downstream analysis
#
# The gffutils relations table stores this hierarchical information but
# nothing else. To find out information about what those IDs refer to, we
# find the relevent features and attribute table rows (for both parent and
# child) and use these to filter the results to include children which are
# CDSs with the protein_coding transcript type (thus removing psuedo genes
# and the like). We then need to specify that of these only pick those whose
# parent is a gene (the relations table keeps an entry for each possible
# hierarchical connection)
#
# Finally, once we have our list of protein coding CDS and gene pairings, we
# group the CDSs by gene, using the CDS start and end columns to determine
# the length and selecting the attribute with the longest sequence
#
query = """
    SELECT DISTINCT
        relations.parent,
        relations.child,
        MAX(childFeat.end - childFeat.start) AS length
    FROM
        relations
        INNER JOIN features childFeat      ON childFeat.id    == relations.child
        INNER JOIN features parentFeat     ON parentFeat.id   == relations.parent
        INNER JOIN attributes childAttrib  ON childAttrib.id  == relations.child
        INNER JOIN attributes parentAttrib ON parentAttrib.id == relations.parent
    WHERE
        childFeat.featuretype       == "CDS" AND
        parentFeat.featuretype      == "gene" AND
        childAttrib.transcript_type == "protein_coding"
    GROUP BY relations.parent
"""

# We count the number of results we will get (for tqdm)
# Then we extract these values
count  = annotations.execute(f"SELECT COUNT() FROM ({query})").fetchone()[0]
result = annotations.execute(query)

result = [ dict(res) for res in tqdm(
    result,
    total = count,
    desc  = "Extracting Representative CDS IDs"
) ]
result = pd.DataFrame.from_records(result)
result.columns = ["gene", "cds", "length"]

result.to_csv("data/qc/proteins/human.csv", index=False)
