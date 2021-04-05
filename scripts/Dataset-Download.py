from enum import Enum
from tqdm.auto import tqdm

import numpy as np

import pathlib

class EntrezXML:
    class Action(Enum):
        search  = "esearch"
        summary = "esummary"
        info    = "einfo"

    def __init__(self, email, api = None):
        self._email = email
        self._api   = api


    lastRequest = 0    # When we last requested data
    delay       = 0.37 # Can send 3 requests per sec w/o API (so delay for little)

    @staticmethod
    def __runAction(action, progress, **params):
        """
            Sends the request to Entrez

        Parameters:
            action -- Action by which to process the request
            params -- Parameters to pass to Entrez
        """

        import requests
        import time
        # Can only make 3 requests per second
        # So if we are making requests too quickly we need to slow down

        if(EntrezXML.lastRequest - time.time() > EntrezXML.delay):
            wait = (EntrezXML.lastRequest + EntrezXML.delay) - time.time()
            time.sleep(wait)
        EntrezXML.lastRequest = time.time()

        base      = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
        actionUrl = base + f"{action.value}.fcgi"

        if(progress):
            response   = requests.get(actionUrl, params=params, stream=True)
            block_size = 1024
            file_size  = int(response.headers.get('Content-Length', 0))

            progress = tqdm(
                desc       = f"Entrez {action.name}",
                total      = file_size if (file_size > 0) else None,
                unit       = 'iB',
                unit_scale = True
            )

            result = b''
            for data in response.iter_content(block_size):
                progress.update(len(data))
                result += data

            return result.decode("utf-8")
        else: return requests.get(actionUrl, params=params).text

    @staticmethod
    def __parseNode(node):
        ''' Used to recursively parse the XML DOM '''
        if(len(list(node)) == 0):
            try:    return {node.tag: int(node.text)}
            except: return {node.tag: node.text}
        else:
            children = [EntrezXML.__parseNode(child) for child in node]
            attrib   = [] if node.attrib == {} else [node.attrib]

            return { node.tag: attrib + children }

    @staticmethod
    def __toXML(response):
        import xml.etree.ElementTree as xml
        return xml.fromstring(response)

    @staticmethod
    def __parse(response):
        ''' Converts the result into a JSON-like python object'''
        response_xml = EntrezXML.__toXML(response)
        as_dict = EntrezXML.__parseNode(response_xml)
        return as_dict

    @staticmethod
    def __process(action, result):
        ''' Performs additional processing on parsed results'''
        import collections
        if(action == EntrezXML.Action.search):
            [(rootName, values)] = result.items()
            return dict(collections.ChainMap(*values))
        else: return result

    def get(self, action, progress = True, output = "python", **params):
        """
        Sends a request to performa specific action using the Entrez Utility APIs

        Paramaters:
            action -- The action to be used for this request
                      The actions include:
                        - Entrez.Action.search:  Send a eSearch GET request
                        - Entrez.Action.summary: Sends a esummary GET request
                        - Entrez.Action.info:    Sends a einfo GET request
            progress -- Whether or not to display the tqdm progress bar
                        while downloading the results
            output -- Determines how the result is outputted:
                        - raw => Raw result received are returned
                        - native => Skips the post-processing of the results into a JSON-like python object
                        - python => Creates a nested JSON-like python object
        """

        # Add additional parameters required by the Entrez API
        if("tool" not in params):
            params["tool"] = "python"
        if("email" not in params and self._email is not None):
            params["email"] = self._email
        if("api_key" not in params and self._api is not None):
            params["api_key"] = self._api

        response = self.__runAction(action, progress, **params)

        if(output == "raw"): return response
        elif(output == "native"): return self.__toXML(response)
        elif(output == "python"):
            xml = self.__parse(response)
            return self.__process(action, xml)

Entrez = EntrezXML("ss2980@bath.ac.uk")

entrez_search = Entrez.get(
    EntrezXML.Action.search,
    db   = "assembly",
    term = (
        '("Bacteria"[Organism] OR "Archaea"[Organism]) AND '
        'latest_refseq[filter] AND '
        '(latest[filter] AND all[filter] NOT anomalous[filter])'
    ),
    useHistory = "y"
)

# Roughly 52.4 MiB download
entrez_summary = Entrez.get(
    EntrezXML.Action.summary,
    output    = "native",
    query_key = entrez_search["QueryKey"],
    WebEnv    = entrez_search["WebEnv"],
)

import pandas as pd

summary = pd.DataFrame(
    [{
        "uid": item.attrib["uid"], # Get the list of UIDs used by NCBI
        "accession_genbank": item.find("Synonym").findtext("Genbank"),
        "accession_refseq":  item.find("Synonym").findtext("RefSeq")
    } for item in entrez_summary.findall(".//DocumentSummary")]
)
pathlib.Path("data/genomes/build/").mkdir(parents=True, exist_ok=True)
summary.to_csv("data/genomes/build/ncbi.csv", index=False)

# We nearly have everything, we just need sequence length
# We can already do this with ENA so will just pull that and merge the results

def getENAEndpoint(name):
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
        f'&query=tax_tree({superkingdom[name]})'    # We want to retrieve accession numbers for specific domain/superkingdom(s)
        '&limit=0'                                  # We want every ID we can get a hold of (default: 100000)
        '&fields=accession,tax_id,base_count'       # We want to have: accession id, NCBI taxanomic id and base count
        '&format=tsv'                               # Gives our data back as a tsv
    )
    # These were useful for working out what filters to apply:
    # - List of the Data Classes: https://ena-docs.readthedocs.io/en/latest/retrieval/general-guide/data-classes.html
    # - List of the Taxonomic Divisions: https://www.ncbi.nlm.nih.gov/genbank/htgs/table1/

def downloadFile(domain):
    import requests, pathlib
    from tqdm.contrib import tenumerate

    url      = getENAEndpoint(domain)
    response = requests.get(url, stream=True)
    # Estimates the number of bar updates
    block_size = 1024
    file_size  = int(response.headers.get('Content-Length', 0))

    progress = tqdm(
            position = 0 if(domain == "archaea") else 1,
            leave = True,
            total = file_size,
            desc = f"Dowloading {domain} IDs",
            unit = 'iB',
            unit_scale=True
    )

    with open(f"data/genomes/build/{domain}-ena.tsv", 'wb') as file:
        for data in response.iter_content(block_size):
            progress.update(len(data))
            file.write(data)

from multiprocessing import Pool, RLock, freeze_support
from tqdm import tqdm
import os
if __name__ == "__main__":
    freeze_support()        # for Windows support
    tqdm.set_lock(RLock())  # for managing output contention
    with Pool(
        processes=2,
        initializer=tqdm.set_lock,
        initargs=(tqdm.get_lock(),)
    ) as pool: pool.map(downloadFile, ["bacteria", "archaea"])

summary["accession"] = summary["accession_genbank"].str.split(".", expand = True)[0]

for domain in ["bacteria", "archaea"]:
    ena_result = pd.read_csv(f"data/genomes/build/{domain}-ena.tsv", delimiter='\t')
    result = pd.merge(
        ena_result, summary,
        on = "accession"
    )
    result.to_csv(f"data/genomes/build/{domain}-result.csv")
