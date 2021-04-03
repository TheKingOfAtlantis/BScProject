
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

# Url to ebi API endpoint to retrieve list of all the assemblies
# Specifically we ask for:
ids_url = (
    'https://www.ebi.ac.uk/ena/portal/api/search?'
    'result=sequence'                          # We want to accession IDs for sequences (not the assemblies)
    '&limit=0'                                 # We want every ID we can get a hold of (default: 100000)
    '&query='
        f'tax_tree({superkingdom[name]}) AND ' # We want to retrieve accession numbers for specific superkingdom(s)
        'mol_type="genomic DNA"'               # Ensure that we are retrieving genomic DNA
    '&fields=accession,tax_id,base_count'      # We want to have: accession id, NCBI taxanomic id and base count
    '&format=tsv'                              # Gives our data back as a tsv
)

# Make a GET request to retrieve the data before converting to a DataFrame
ids = requests.get(ids_url)
ids = pd.read_csv(io.StringIO(ids.content.decode("utf-8")), sep="\t")

# We also want to retrieve the taxaDB we generated
# 1) We use to determine the lineage of each taxa ID
# 2) Identify the superkingdom of each taxa ID
# 3) To filter each genus down to a single represenative genome
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
    # Right now we want a list of taxa which fulfill our previous requirements
    # We can then perform a reverse search on the table of IDs to find relevent accession IDs
    with Pool(os.cpu_count()) as pool:
        # Some taxa ids referred to by multiple accession no.
        # So we filter our list to a unique set of taxa IDs
        # Before filtering that list to just those in the (eu)bacteria superkingdom
        taxaIds = pd.unique(ids["tax_id"])

        bacteria_taxa = list(tqdm.tqdm(pool.imap(
            filterSuperkingdom, zip(
            taxaIds,
            itertools.repeat("bacteria")
        )), total = len(taxaIds)))
        bacteria_taxa = list(filter(lambda x: x != None, bacteria_taxa))

        # Now we need to group by genus and select just one species to use
        # First lets get the lineage of each bacteria taxa ID
        bacteria_taxa = pd.DataFrame({
            "tax_id": bacteria_taxa
        })
        bacteria_taxa["lineage"] = list(tqdm.tqdm(pool.imap(
            getLineage,
            bacteria_taxa.tax_id
        ), total = len(bacteria_taxa)))

    # Combine those results in a dataframe
    def preprocessor(block):
        return {x[0]:pd.Series(x[1]) for x in block}
    bacteria_taxa = concat(bacteria_taxa.values, preprocessor = preprocessor, axis = 1).T
    bacteria_taxa = bacteria_taxa.rename_axis("tax_id").reset_index()

    # Since not all tax IDs have an associated genus tax ID
    # We must either:
    # 1) Find another method to reduce phylogenetic influence from the list
    #    this may include filtering at a higher rank (e.g. family) whereever genus is not present
    # 2) Or remove them as our current method works by filtering for one representative species for a genus
    # For now we will simply remove them potentially adjusting our algorithm later

    # Finally we merge our list of lineages with the original ids list
    # This takes the intersection of both tables so will will only keep rows where we have entries in both
    bacteria_ids = pd.merge(
        ids, bacteria_taxa[bacteria_taxa.genus.isnull() == False],
        on = "tax_id",
        validate = "m:1" # Ensure its a many-to-one relation
    )

    # Select the species in each genus with the largest base count
    bacteria = bacteria_ids.loc[
        bacteria_ids.groupby("genus", sort=False)["base_count"].idxmax()
    ]

    # As a sanity check we'll assert that we want only 1 species per genus
    assert not (bacteria.groupby("genus").count()["species"] != 1).any()

    # As a reference lets keep a list of the accession no. of each bacteria
    bacteria = bacteria[["accession", "tax_id", "base_count", "genus"]]
    bacteria["sci_name"] = bacteria["tax_id"].apply(taxaDB.sci_name)
    bacteria.to_csv("data/genomes/bacteria.csv",index=False)
