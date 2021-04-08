#!/bin/bash

python scripts/DatasetDownload.py # Fetchs the metadata needed to select the assemblies we need
python scripts/DatasetFilter.py   # Filters the assemblies based on the metadata
python scripts/DatasetFetch.py    # Downloads the assemblies which have been selected

cp -r data/genomes/build/out/* data/genomes/

if [[ $1 =~ [Nn]$ ]]; then
    echo "Build files have been kept"
elif [[ $1 =~ [Yy]$ ]]; then
    REPLY=y
else
    read -p "Do you want to delete build files? " -n 1 -r
    echo
fi
if [[ $REPLY =~ [Yy]$ ]]; then
    echo Deleting build files
    rm -r data/genomes/build/
fi
