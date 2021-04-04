
# Pull IDs from
# https://www.ebi.ac.uk/genomes/bacteria.details.txt
# https://www.ebi.ac.uk/genomes/bacteria.txt

# Pull data with
# http://www.ebi.ac.uk/Tools/dbfetch/emblfetch?db=embl&id=$accNo&format=embl&style=raw&Retrieve=Retrieve

###

# Can pul from the new database:
# Assemblies: https://www.ebi.ac.uk/ena/portal/api/search?query=genome_representation=%22full%22&result=assembly&fields=accession,tax_id,base_count&format=tsv


import requests
import io, os, tqdm, pathlib
import more_itertools, itertools
from multiprocessing import Pool

import pandas as pd

from taxadb.taxid import TaxID

def getDB(): return TaxID(dbtype='sqlite', dbname='data/taxadb.sqlite')

def getEndpoint(name):
    superkingdom = {
        "bacteria": 2,
        "archaea": 2157,
        "eukaryota": 2759
    }

    # Url to EBI ENA API endpoint to retrieve list of all the sequencies
    # Specifically we ask for:
    return (
        'https://www.ebi.ac.uk/ena/portal/api/search?'
        'result=assembly'
        '&query='
            'assembly_level="chromosome" AND '      # Get assemblies for chromosomes only
            'genome_representation="full" AND '     # We want the whole genome not just subdivisions of the genome (e.g. chromosomes)
            f'tax_tree({superkingdom[name]})'       # We want to retrieve accession numbers for specific superkingdom(s)
        '&limit=0'                                  # We want every ID we can get a hold of (default: 100000)
        '&fields=accession,tax_id,base_count'       # We want to have: accession id, NCBI taxanomic id and base count
        '&format=tsv'                               # Gives our data back as a tsv
    )
    # These were useful for working out what filters to apply:
    # - List of the Data Classes: https://ena-docs.readthedocs.io/en/latest/retrieval/general-guide/data-classes.html
    # - List of the Taxonomic Divisions: https://www.ncbi.nlm.nih.gov/genbank/htgs/table1/

# We also want to retrieve the taxaDB we generated
# 1) We use to determine the lineage of each taxa ID
# 2) Identify the superkingdom of each taxa ID
# 3) To filter each genus down to a single represenative genome

def filterSuperkingdom(x):
    taxid, name = x
    return __filterSuperkingdom(taxid, name)
def __filterSuperkingdom(taxid, name):
    taxaDB = getDB()

    superkingdom = {
        "bacteria": 2,
        "archaea": 2157,
        "eukaryota": 2759
    }
    if(taxaDB.has_parent(taxid, superkingdom[name])):
        return taxid

def getLineage(taxid):
    taxaDB = getDB()
    out = taxaDB.lineage_id(taxid, ranks = True)
    if(out is None): return { "None": taxid }
    else: return dict(out)

from common import concat
def lineageConcat(block):
    return {x[0]:pd.Series(x[1], dtype='int') for x in block}

if __name__ == "__main__":
    for superkingdom in ["bacteria"]: #["bacteria", "archaea"]:
        # Make a GET request to retrieve the data before converting to a DataFrame
        ena_result = requests.get(getEndpoint(superkingdom))
        ena_result = pd.read_csv(io.StringIO(ena_result.content.decode("utf-8")), sep="\t")

        # Right now we want a list of taxa which fulfill our previous requirements
        # We can then perform a reverse search on the table of IDs to find relevent accession IDs
        with Pool(os.cpu_count()) as pool:
            # Some taxa ids referred to by multiple accession no.
            # So we filter our list to a unique set of taxa IDs
            # Before filtering that list to just those in the (eu)bacteria superkingdom
        taxaIds = pd.unique(ena_result.tax_id)

            taxa_result = list(tqdm.tqdm(pool.imap(
                filterSuperkingdom, zip(
                taxaIds,
                itertools.repeat("bacteria")
            )), total = len(taxaIds)))
            taxa_result = list(filter(lambda x: x != None, taxa_result))

            # Now we need to group by genus and select just one species to use
            # First lets get the lineage of each bacteria taxa ID
            taxa_result = pd.DataFrame({
                "tax_id": taxa_result
            })
            taxa_result["lineage"] = list(tqdm.tqdm(pool.imap(
                    getLineage,
                taxa_result.tax_id
            ), total = len(taxa_result)))

        # Combine those results in a dataframe
        taxa = concat(
            taxa_result.values,
            lineageConcat, axis = 1
        ).T.rename_axis("tax_id").reset_index()

        # Since not all tax IDs have an associated genus tax ID
        # We must either:
        # 1) Find another method to reduce phylogenetic influence from the list
        #    this may include filtering at a higher rank (e.g. family) whereever genus is not present
        # 2) Or remove them as our current method works by filtering for one representative species for a genus
        # For now we will simply remove them potentially adjusting our algorithm later

        # Finally we merge our list of lineages with the original ids list
        # This takes the intersection of both tables so will will only keep rows where we have entries in both
        ids = pd.merge(
            ena_result, taxa[taxa.genus.isnull() == False],
            on       = "tax_id",
            validate = "m:1" # Ensure its a many-to-one relation
        )

        # Select the species in each genus with the largest base count
        largest = ids.loc[ids.groupby("genus", sort=False)["base_count"].idxmax()]

        # As a sanity check we'll assert that we want only 1 species per genus
        assert not (largest.groupby("genus").count()["species"] != 1).any()

        # As a reference lets keep a list of the accession no. of each bacteria
        taxaDB = getDB()
        out = largest[["accession", "tax_id", "base_count", "genus"]]
        out["sci_name"] = out["tax_id"].apply(taxaDB.sci_name)

        pathlib.Path("data/genomes/").mkdir(parents=True, exist_ok=True)
        out.to_csv(f"data/genomes/{superkingdom}.csv",index=False)
