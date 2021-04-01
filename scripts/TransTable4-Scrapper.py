
# Pull IDs from
# https://www.ebi.ac.uk/genomes/bacteria.details.txt
# https://www.ebi.ac.uk/genomes/bacteria.txt

# Pull data with
# http://www.ebi.ac.uk/Tools/dbfetch/emblfetch?db=embl&id=$accNo&format=embl&style=raw&Retrieve=Retrieve

###

# Can pul from the new database:
# Assemblies: https://www.ebi.ac.uk/ena/portal/api/search?query=genome_representation=%22full%22&result=assembly&fields=accession,tax_id,base_count&format=tsv

import pandas as pd
import requests, io, os, tqdm

from taxadb.taxid import TaxID

ids_url = "https://www.ebi.ac.uk/ena/portal/api/search?result=assembly&query=genome_representation%3D%22full%22&fields=accession%2Ctax_id%2Cbase_count&format=tsv"
ids = requests.get(ids_url)
ids = pd.read_csv(io.StringIO(ids.content.decode("utf-8")), sep="\t")

taxaDB = TaxID(dbtype='sqlite', dbname='data/taxadb.sqlite')

def filterSuperkingdom(x):
    taxid, name = x
    return __filterSuperkingdom(taxid, name)
def __filterSuperkingdom(taxid, name):
    superkingdom = {
        "bacteria": 2,
        "archaea": 2157,
        "eukaryota": 2759
    }
    if(taxaDB.has_parent(taxid, superkingdom[name])):
        return taxid

def getLineage(taxaid):
    return dict(taxaDB.lineage_id(taxaid, ranks = True))

from multiprocessing import Pool
from common import concat
import itertools
if __name__ == "__main__":
    with Pool(os.cpu_count()) as pool:
        taxaIds = pd.unique(ids["tax_id"])

        bacteria_taxa = list(tqdm.tqdm(pool.imap(
            filterSuperkingdom, zip(
            taxaIds,
            itertools.repeat("bacteria")
        )), total = len(taxaIds)))
        bacteria_taxa = list(filter(lambda x: x != None, bacteria_taxa))

        bacteria_taxa = pd.DataFrame({
            "tax_id": bacteria_taxa
        })
        bacteria_taxa["lineage"] = list(tqdm.tqdm(pool.imap(
            getLineage,
            bacteria_taxa["tax_id"]
        ), total = len(bacteria_taxa)))

    def preprocessor(block):
        return {x[0]:pd.Series(x[1]) for x in block}
    bacteria_taxa = concat(bacteria_taxa.values, preprocessor = preprocessor, axis = 1).T
    bacteria_taxa = bacteria_taxa.rename_axis("tax_id").reset_index()
    bacteria_taxa = bacteria_taxa[[
        "tax_id",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "clade",
        "superkingdom",
        "no rank"
    ]]

    bacteria_ids = ids.merge(
        bacteria_taxa,
        on = "tax_id",
        validate = "m:1" # Ensure its a many-to-one relation
    )

