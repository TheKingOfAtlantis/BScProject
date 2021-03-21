
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

def filterBacteria(taxid):
    if(taxaDB.has_parent(taxid, 2)): # Bacteria => Taxonomy ID: 2
        return taxid

from multiprocessing import Pool
if __name__ == "__main__":
    with Pool(os.cpu_count()) as pool:
        taxaIds = pd.unique(ids["tax_id"])
        bacteria_taxa = list(tqdm.tqdm(pool.imap(
            filterBacteria,
            taxaIds
        ), total = len(taxaIds)))
bacteria_taxa = list(filter(lambda x: x != None, bacteria_taxa))
bacteria_ids = ids[ids["tax_id"].isin(bacteria_taxa)]
