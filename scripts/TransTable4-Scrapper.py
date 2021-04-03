
# Pull IDs from
# https://www.ebi.ac.uk/genomes/bacteria.details.txt
# https://www.ebi.ac.uk/genomes/bacteria.txt

# Pull data with
# http://www.ebi.ac.uk/Tools/dbfetch/emblfetch?db=embl&id=$accNo&format=embl&style=raw&Retrieve=Retrieve

###

# Can pul from the new database:
# Assemblies: https://www.ebi.ac.uk/ena/portal/api/search?query=genome_representation=%22full%22&result=assembly&fields=accession,tax_id,base_count&format=tsv


import requests
import io, os, tqdm
import more_itertools
from multiprocessing import Pool

import pandas as pd

from taxadb.taxid import TaxID

def getEndpoint(name):
    superkingdom = {
        "bacteria": 2,
        "archaea": 2157,
        "eukaryota": 2759
    }

    # Url to EBI ENA API endpoint to retrieve list of all the sequencies
    return (
        'https://www.ebi.ac.uk/ena/portal/api/search?'
        'result=sequence'                          # We want to accession IDs for sequences (not the assemblies)
    '&limit=0'                                 # We want every ID we can get a hold of (default: 100000)
        '&query='
            f'tax_tree({superkingdom[name]}) AND ' # We want to retrieve accession numbers for specific superkingdom(s)
            'mol_type="genomic DNA"'               # Ensure that we are retrieving genomic DNA
        '&fields=accession,tax_id,base_count'      # We want to have: accession id, NCBI taxanomic id and base count
        '&format=tsv'                              # Gives our data back as a tsv
    )
    # These were useful for working out what filters to apply:
    # - List of the Data Classes: https://ena-docs.readthedocs.io/en/latest/retrieval/general-guide/data-classes.html
    # - List of the Taxonomic Divisions: https://www.ncbi.nlm.nih.gov/genbank/htgs/table1/


import xml.etree.ElementTree as xml
def getLineage(chunk):
    taxaUrl = "https://www.ebi.ac.uk/ena/browser/api/xml/" + ",".join(chunk)
    taxa_result = requests.get(taxaUrl)
    taxa_result.raise_for_status()
    taxa_result = taxa_result.content.decode("utf-8")

    taxa_result = xml.fromstring(taxa_result)

    return pd.DataFrame.from_dict([{
        "tax_id":      taxa.attrib["taxId"],              # We need to this to merge the data properly
        "sci_name":    taxa.attrib["scientificName"],        # Want to keep a record of this
        "trans_table": taxa.attrib["geneticCode"],        # Since its avaliable and something we care about we'll keep this too
        "lineage": pd.DataFrame([                         # Finally we retrieve the lineage for this tax id
            taxon.attrib for taxon in taxa.find("lineage")
        ]).drop(columns=["hidden", "commonName"])         # We don't need these columns
            .rename_axis("order")                         # Just so we have a reference to the rank order
            .reset_index()                                # Reset the ID since we want it as a column
            .dropna()                                     # Some ranks have no name (but we don't need these, so we can lose them)
    } for taxa in taxa_result])

from common import concat
def lineageConcat(block):
    return {x[0]:x[1][["rank","taxId"]].set_index("rank") for x in block}

if __name__ == "__main__":
    for superkingdom in ["bacteria"]: #["bacteria", "archaea"]:
        # Make a GET request to retrieve the data before converting to a DataFrame
        ena_result = requests.get(getEndpoint(superkingdom))
        ena_result = pd.read_csv(io.StringIO(ena_result.content.decode("utf-8")), sep="\t")

        # To minimise the number of requests needed as well as to avoid duplication of requests
        # We reduce the list of taxa ids to be unique
        taxaIds = pd.unique(ena_result.tax_id)

        # Once we have our reduced list we start pooling the ENA database for taxonomic data
        # To avoid 400 errors the requests are made in chunks of 1000 ids
        with Pool(os.cpu_count()) as pool:
            chunk_size  = 1000 # Seems to be the largest no. of ids we can send at once
            taxa_result = pd.concat(
                tqdm.tqdm(pool.imap(
                    getLineage,
                    more_itertools.chunked(map(str, taxaIds), chunk_size)
                ), total = len(taxaIds)//chunk_size + 1),
                ignore_index=True
            ).astype({'tax_id': 'int'})

        taxa = concat(
            taxa_result[["tax_id", "lineage"]].values,
            lineageConcat, axis = 1
        ).T.droplevel(1)\
           .rename_axis("tax_id")\
           .reset_index()

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

        out = pd.merge(
            taxa_result[["tax_id", "sci_name"]],
            largest,
            on = "tax_id"
        )[["accession", "tax_id", "genus", "sci_name", "base_count"]]
        out.to_csv(f"data/genomes/{superkingdom}.csv",index=False)
