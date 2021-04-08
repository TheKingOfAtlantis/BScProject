import pandas as pd
from EntrezCommon import Entrez
from common import Filesystem

import pathlib

Database = "assembly"
entrez_search = Entrez.get(
    Entrez.Action.search,
    db   = Database,
    term = (
        'prokaryota[orgn] AND '
        '"representative genome"[refseq category] AND '
        '((latest[filter] OR "latest refseq"[filter]) AND all[filter] NOT anomalous[filter]) AND'
        '"complete genome"[filter]'
    ),
    useHistory = "y"
)

entrez_summary = Entrez.getAll(
    Entrez.Action.summary,
    output    = "native",
    db        = Database,
    count     = entrez_search["Count"],
    query_key = entrez_search["QueryKey"],
    WebEnv    = entrez_search["WebEnv"]
)

def LengthFromCDATA(cdata):
    import xml.etree.ElementTree as xml
    import re

    return int(re.findall('<Stat category="total_length" sequence_tag="all">(\d+)</Stat>', cdata)[0])

summary = pd.DataFrame([{
    "uid":               item.attrib["uid"],
    "accession":         item.find("Synonym").findtext("Genbank").split(".")[0],
    "accession_genbank": item.find("Synonym").findtext("Genbank"),
    "accession_refseq":  item.find("Synonym").findtext("RefSeq"),
    "status":            item.findtext("AssemblyStatus"),
    "tax_id":            item.findtext("Taxid"),
    "length":            LengthFromCDATA(item.findtext("Meta")),
    "link":              item.findtext("FtpPath_RefSeq")
} for result in entrez_summary for item in result.findall(".//DocumentSummary")])

# Save space by keeping a truncated list
summary.link = summary.link.str.replace("ftp://ftp.ncbi.nlm.nih.gov/genomes", "", regex=False)

Filesystem.mkdir("data/genomes/build/")
summary.to_csv("data/genomes/build/ncbi.csv", index=False)
