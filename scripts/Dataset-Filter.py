
import os, pathlib
import itertools, more_itertools

import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool

from common import loadGlob, concat

def getDB():
    from taxadb.taxid import TaxID
    return TaxID(dbtype='sqlite', dbname='data/taxadb.sqlite')

def getLineage(taxid):
    taxaDB = getDB()
    out = taxaDB.lineage_id(taxid, ranks = True)
    if(out is None): return { "missing": taxid }
    else: return dict(out)

def lineageConcat(block):
    return {x[0]:pd.Series(x[1], dtype='int') for x in block}

def process(superkingdom):
    print(superkingdom)
    ena_result = pd.read_csv(f"data/genomes/{superkingdom}-ids.tsv", delimiter='\t')
    with Pool(os.cpu_count()) as pool:
        # Some taxa ids referred to by multiple accession no.
        # So we filter our list to a unique set of taxa IDs

        # Now we need to group by genus and select just one species to use
        # First lets get the lineage of each bacteria taxa ID
        taxa_result = pd.DataFrame({ "tax_id": pd.unique(ena_result.tax_id) })
        taxa_result["lineage"] = list(tqdm(pool.imap(
            getLineage,
            taxa_result.tax_id
        ), total = len(taxa_result), desc="Retrieving lineages"))

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

    # As a sanity check we'll assert that we want only 1 accession per genus
    assert not (largest.groupby("genus").count()["accession"] != 1).any()

    # As a reference lets keep a list of the accession no. of each bacteria
    taxaDB = getDB()
    out = largest.loc[,["accession", "tax_id", "base_count", "genus"]]
    out["sci_name"] = out["tax_id"].apply(taxaDB.sci_name)

    pathlib.Path("data/genomes/").mkdir(parents=True, exist_ok=True)
    out.to_csv(f"data/genomes/{superkingdom}.csv",index=False)

for domain in ["archaea", "bacteria"]:
    process(domain)