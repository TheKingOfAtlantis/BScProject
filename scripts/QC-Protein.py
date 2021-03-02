# Checks if all the genes which are used in gene-GC.py are well-behaved protein
# coding sequences:
# - Checking sequence length multiple of 3
# - Check no internal stop codons
# - Check ends with a stop codon

from Bio import SeqIO
from Bio.Data import CodonTable


def getCodonTable(feature):
    return CodonTable.generic_by_id[int(feature.qualifiers["transl_table"][0])]

# TODO: Consider moving this into common - Would be useful in many other scripts
def isStop(codon, feature):
    '''
        Given a codon checks if it is a stop codon
    '''
    codonTable = getCodonTable(feature)
    return codon in codonTable.stop_codons

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
    return isStop(endCodon, feature)
def internalCheck(seq, feature):
    '''
        Checks the internal segment of the sequence to ensure that no stop codons are
        present within the sequence beyond the stop codon used to terminate it
    '''
    featureSeq = feature.extract(seq)
    for i in range(3, len(feature), 3):
        if(isStop(featureSeq[i-3:i], feature)):
            # If we encounter a stop then it must be internal
            return False
    return True

def check(file):
    ''' Collects the results of all the checks '''
    record = next(SeqIO.parse(file, "embl"))
    for feature in filter(lambda x: x.type in "CDS", record.features):
        return (record.id, {
            "length":   lengthCheck(feature),               # Result of checking length of sequence (If sequence multiple of three == True)
            "start":    startCheck(record.seq, feature),    # Result of checking start is start codon (If start wtih start == True)
            "end":      endCheck(record.seq, feature),      # Result of checking end is stop codon (If ends with stop == True)
            "internal": internalCheck(record.seq, feature), # Result of checking for internal stop codons (If none found == True)
        })

if __name__ == "__main__":
    from common import loadGlob
    import pandas as pd
    import itertools, pathlib

    pathlib.Path("data/qc/proteins/").mkdir(parents=True, exist_ok=True)

    # Generate results
    # Then check results for any which failed any 1 of the criteria

    result = dict(loadGlob("data/genomes/archaea/*", check))
    result = pd.DataFrame(result).T
    result.to_csv("data/qc/proteins/archaea.csv")

    result = dict(loadGlob("data/genomes/bacteria/*", check))
    result = pd.DataFrame(result).T
    result.to_csv("data/qc/proteins/bacteria.csv")
