
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
        'result=sequence'
        '&query='
            ''
            'mol_type="genomic DNA" AND '
            f'tax_tree({superkingdom[name]})'       # We want to retrieve accession numbers for specific superkingdom(s)
        '&limit=0'                                  # We want every ID we can get a hold of (default: 100000)
        '&fields=accession,tax_id,base_count'       # We want to have: accession id, NCBI taxanomic id and base count
        '&format=tsv'                               # Gives our data back as a tsv
    )
    # return (
    #     'https://www.ebi.ac.uk/ena/portal/api/search?'
    #     'result=assembly'
    #     '&query='
    #         'assembly_level="chromosome" AND '      # Get assemblies for chromosomes only
    #         'genome_representation="full" AND '     # We want the whole genome not just subdivisions of the genome (e.g. chromosomes)
    #         f'tax_tree({superkingdom[name]})'       # We want to retrieve accession numbers for specific superkingdom(s)
    #     '&limit=0'                                  # We want every ID we can get a hold of (default: 100000)
    #     '&fields=accession,tax_id,base_count'       # We want to have: accession id, NCBI taxanomic id and base count
    #     '&format=tsv'                               # Gives our data back as a tsv
    # )
    # These were useful for working out what filters to apply:
    # - List of the Data Classes: https://ena-docs.readthedocs.io/en/latest/retrieval/general-guide/data-classes.html
    # - List of the Taxonomic Divisions: https://www.ncbi.nlm.nih.gov/genbank/htgs/table1/


def downloadFile(domain):
    import requests
    from tqdm.contrib import tenumerate

    url      = getEndpoint(domain)
    response = requests.get(url, stream=True)
    # Estimates the number of bar updates
    block_size = 1024
    file_size  = int(response.headers.get('Content-Length', 0))

    with open(f"data/genomes/{domain}-ids.tsv", 'wb') as file:
        for i, data in tenumerate(
            response.iter_content(block_size),
            position = 0 if(domain == "archaea") else 1,
            leave = True,
            total = file_size,
            desc = f"Dowloading {domain} IDs",
            unit = 'iB',
            unit_scale=True
        ): file.write(data)

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
