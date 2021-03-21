cd data
taxadb download -t taxa -o taxadb
taxadb create -d taxa -i taxadb -n taxadb.sqlite
rm -r taxadb
