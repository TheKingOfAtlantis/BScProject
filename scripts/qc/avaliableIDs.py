
# Want to check what ID types are avaliable to keep a record of genes
# Plan:
#  - Pull all the ID types possible from each genome
#  - Shrink the list available in each genome down (w/i the seperate processes)
#  - Shrink the list available across all genomes

import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from common import Filesystem

import functools, collections, itertools

def pullGeneIds(feature):
    idTypes = []
    qualifiers = feature.qualifiers

    # Assume we have nothing
    # So we have to check each is present (can't take it for granted)

    if("locus_tag" in qualifiers):  idTypes.append("locus_tag")
    if("protein_id" in qualifiers): idTypes.append("protein_id")
    if("gene" in qualifiers): idTypes.append("gene")
    if("db_xref" in qualifiers):
        # With external refs need to seperate ID from ID type
        for id in qualifiers["db_xref"]:
            idTypes.append("db_xref:" + id.split(":")[0])

    if(idTypes == []):
        idTypes.extend(list(qualifiers.keys()))

    # We can't really use the "note" qualifier as an ID
    if("note" in idTypes):
        idTypes.remove("note")

    return set(idTypes) # Return as set type to ensure we unique set of ID types

def shrinkGenomeIds(file, toFind):
    # First we store everything we get (even repeats)
    # Then we reduce
    from Bio import SeqIO

    record = next(SeqIO.parse(file, "embl"))

    ids = []
    for feature in filter(lambda x: x.type in toFind, record.features):
        res = pullGeneIds(feature)
        if(res == set()): print(record.id, feature)
        else: ids.append(res)
    return ids

if __name__ == "__main__":
    idTypes = Filesystem.loadGlob("data/genomes/*/*.embl", shrinkGenomeIds, toFind=["CDS", "tRNA"], desc = "Checking ID usage")

    # Plan to find minimal set
    # 1) Count frequence of each Id types usage
    # 2) Select highest used ID type
    # 3) Remove all which contain the ID type (these can be addressed with an ID of that type)
    # 4) Repeat on set of unaddressable genes, until all can be addressed


    # Use an OrderedDict (as an OrderedSet doesn't exist)
    # Use the order to keep record of how many elements can be addressed
    # With the first able to address the most - Order will be used to determine
    # in what order to check for each ID

    print("Determining minimum ID set")
    minimal = collections.OrderedDict()
    remains = list(itertools.chain.from_iterable(idTypes))

    while len(remains) != 0:
        count = collections.Counter(itertools.chain.from_iterable(remains))
        minimal[count.most_common(1)[0][0]] = None

        remains = list(filter(lambda x: set(minimal.keys()).isdisjoint(x), remains))

    with open("data/qc/minIds.txt", "w") as file:
        file.write("\n".join(list(minimal.keys())))
