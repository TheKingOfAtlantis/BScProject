import pandas as pd
from EntrezCommon import Entrez

import pathlib

Database = "nucleotide"
entrez_search = Entrez.get(
    Entrez.Action.search,
    db   = Database,
    term = (
        '("Bacteria"[Organism] OR "Archaea"[Organism]) AND '
        'refseq[filter] AND '
        '"complete genome"'
    ),
    useHistory = "y"
)

entrez_summary = Entrez.getAll(
    Entrez.Action.summary,
    output    = "native",
    db        = Database,
    count     = entrez_search["Count"],
    query_key = entrez_search["QueryKey"],
    WebEnv    = entrez_search["WebEnv"],
    version   = "2.0"
)

summary = pd.DataFrame([
    {
        "uid":         item.findtext("Gi"), # Get the list of UIDs used by NCBI
        "accession":   item.findtext("AccessionVersion"),
        "tax_id":      item.findtext("TaxId"),
        "base_count":  item.findtext("Slen"),
        "trans_table": item.findtext("GeneticCode"),
        "genome":      item.findtext("Genome"),
    } for result in entrez_summary for item in result.findall(".//DocumentSummary")
])
pathlib.Path("data/genomes/build/").mkdir(parents=True, exist_ok=True)
summary.to_csv("data/genomes/build/ncbi.csv", index=False)
