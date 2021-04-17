
mkdir -p data/genomes/build
cd data/genomes/build

taxadb download -t taxa -o download/taxadb/
taxadb create -d taxa -i download/taxadb/ -n taxadb.sqlite
