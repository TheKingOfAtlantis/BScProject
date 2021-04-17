# Checks if all the genes which are used in gene-GC.py are well-behaved protein
# coding sequences:
# - Checking sequence length multiple of 3
# - Check no internal stop codons
# - Check ends with a stop codon

import sys, pathlib
sys.path.append(str(pathlib.Path(__file__).parent.parent.absolute()))
from common import Filesystem, getID

from Bio import SeqIO
from Bio.Data import CodonTable

# Utility Function(s)

def getCodonTable(feature):
    return CodonTable.generic_by_id[int(feature.qualifiers["transl_table"][0])]

# The checks to be done

def lengthCheck(feature):
    '''
        Checks whether or not the sequence length of this CDS is a multiple of three
        Which indicates it has N codons
    '''
    return len(feature) % 3 == 0
def startCheck(seq, feature):
    '''
        Checks against the codon table for the CDS that it indeed starts with the one
        of its avaliable start codons
    '''
    codonTable = getCodonTable(feature)
    featureSeq = feature.extract(seq)
    startCodon = featureSeq[:3]
    return startCodon in codonTable.start_codons
def endCheck(seq, feature):
    '''
        Checks against the codon table for the CDS that it indeed ends with the one
        of its avaliable end codons
    '''
    featureSeq = feature.extract(seq)
    endCodon   = featureSeq[-3:]
    codonTable = getCodonTable(feature)
    return endCodon in codonTable.stop_codons
def internalCheck(seq, feature):
    '''
        Checks the internal segment of the sequence to ensure that no stop codons are
        present within the sequence beyond the stop codon used to terminate it
    '''
    featureSeq = feature.extract(seq)
    codonTable = getCodonTable(feature)
    for i in range(3, len(feature), 3):
        if(featureSeq[i-3:i] in codonTable.stop_codons):
            # If we encounter a stop then it must be internal
            return False
    return True

def __performCheck(file):
    ''' Collects the results of all the checks '''
    record = next(SeqIO.parse(file, "embl"))
    return (record.id, pd.DataFrame([{
        "id":       getID(feature),                     # Retrieves a suitable ID for identifying the gene
        "length":   lengthCheck(feature),               # Result of checking length of sequence (If sequence multiple of three == True)
        "start":    startCheck(record.seq, feature),    # Result of checking start is start codon (If start wtih start == True)
        "end":      endCheck(record.seq, feature),      # Result of checking end is stop codon (If ends with stop == True)
        "internal": internalCheck(record.seq, feature), # Result of checking for internal stop codons (If none found == True)
    } for feature in filter(lambda x: x.type in "CDS", record.features)]))

def performCheck(pattern, **kwargs):
    data = pd.concat(dict(Filesystem.loadGlob(pattern, __performCheck, **kwargs)))
    data = data.droplevel(1)\
               .set_index("id", append=True)\
               .rename_axis(index=["genome", "id"])
    return data

if __name__ == "__main__":
    import pandas as pd

    Filesystem.mkdir("data/qc/proteins/")

    result = performCheck("data/genomes/archaea/*", desc = "Checking Archaea CDSs")
    result.to_json("data/qc/proteins/archaea.json", orient = "table")

    result = performCheck("data/genomes/bacteria/*", desc = "Checking Bacteria CDSs")
    result.to_json("data/qc/proteins/bacteria.json", orient = "table")
