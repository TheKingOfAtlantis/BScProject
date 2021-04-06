
python scripts/DatasetDownload.py # Fetchs the metadata needed to select the assemblies we need
python scripts/DatasetFilter.py # Filters the assemblies based on the metadata
python scripts/DatasetFetch.py # Downloads the assemblies which have been selected

cp data/genomes/build/out/* data/genomes/

read -p "Do you want to delete data/genomes/build/? " -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    rm -r data/genomes/build/
fi
