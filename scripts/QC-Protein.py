# Checks if all the genes which are used in gene-GC.py are well-behaved protein
# coding sequences:
# - Checking sequence length multiple of 3
# - Check no internal stop codons
# - Check ends with a stop codon

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

def check(file):
    ''' Collects the results of all the checks '''
    record = next(SeqIO.parse(file, "embl"))
    return (record.id, [{
        "position": i, # Position in filtered features list (to find culprits)
        "length":   lengthCheck(feature),               # Result of checking length of sequence (If sequence multiple of three == True)
        "start":    startCheck(record.seq, feature),    # Result of checking start is start codon (If start wtih start == True)
        "end":      endCheck(record.seq, feature),      # Result of checking end is stop codon (If ends with stop == True)
        "internal": internalCheck(record.seq, feature), # Result of checking for internal stop codons (If none found == True)
    } for i, feature in enumerate(filter(lambda x: x.type in "CDS", record.features))])


if __name__ == "__main__":
    from common import loadGlob
    import pandas as pd
    import itertools, pathlib, json

    pathlib.Path("data/qc/proteins/").mkdir(parents=True, exist_ok=True)

    result = dict(loadGlob("data/genomes/archaea/*", check))
    with open("data/qc/proteins/archaea.json", "w") as file: json.dump(result, file)

    result = dict(loadGlob("data/genomes/bacteria/*", check))
    with open("data/qc/proteins/bacteria.json", "w") as file: json.dump(result, file)
