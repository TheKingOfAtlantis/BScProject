
# Want to check what ID types are avaliable to keep a record of genes
# Plan:
#  - Pull all the ID types possible from each genome
#  - Shrink the list available in each genome down (w/i the seperate processes)
#  - Shrink the list available across all genomes

def pullGeneIds(feature):
    idTypes = []
    qualifiers = feature.qualifiers

    # Assume we have nothing
    # So we have to check each is present (can't take it for granted)

    if("locus_tag" in qualifiers):  idTypes.append("locus_tag")
    if("protein_id" in qualifiers): idTypes.append("protein_id")
    if("db_xref" in qualifiers):
        # With external refs need to seperate ID from ID type
        for id in qualifiers["db_xref"]:
            idTypes.append(id.split(":")[0])

    return set(idTypes)

def shrinkGenomeIds(file):
    # First we store everything we get (even repeats)
    # Then we reduce
    from Bio import SeqIO

    idTypes = []
    record = next(SeqIO.parse(file, "embl"))
    for feature in filter(lambda x: x.type in "CDS", record.features):
        idTypes.append(pullGeneIds(feature))

    return functools.reduce(lambda x,y: x & y, idTypes)

from common import loadGlob
import functools
if __name__ == "__main__":
    idTypes = loadGlob("data/genomes/*/*.embl", shrinkGenomeIds)
    print(idTypes)
    idTypes = functools.reduce(lambda x,y: x & y, idTypes)
    print(idTypes)
