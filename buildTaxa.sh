cd data
taxadb download -t full -o taxadb
taxadb create -i taxadb --dbname taxadb.sqlite
rm -r taxadb
