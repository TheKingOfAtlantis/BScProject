cd data
taxadb download -t taxa -o taxadb
taxadb create -i taxadb --dbname taxadb.sqlite
rm -r taxadb
