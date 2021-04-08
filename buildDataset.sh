#!/bin/bash

header=$(tput bold; tput smul)
headerend=$(tput sgr0)

echo ${header} Building Prokaryote Dataset${headerend}
echo

python scripts/dataset/FetchMetadata.py    # Fetchs the metadata needed to select the assemblies we need
python scripts/dataset/Filter.py           # Filters the assemblies based on the metadata
python scripts/dataset/DownloadAssembly.py # Downloads the assemblies which have been selected

echo
echo "${header}Building Human Dataset${headerend}"
echo

url_human_seq="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz"
gz_human_seq="data/genomes/build/download/GRCh38_latest_genomic.fna.gz"
fna_human_seq="data/genomes/build/out/GRCh38_latest_genomic.fna"

url_human_annotation="https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.gff.gz"
gz_human_annotation="data/genomes/build/download/GRCh38_latest_genomic.gff3.gz"
fna_human_annotation="data/genomes/build/out/GRCh38_latest_genomic.gff3"

echo "Downloading Human Genome Sequence"
curl -# $url_human_seq -o $gz_human_seq
echo "Downloading Human Genome Annotation"
curl -# $url_human_annotation -o $gz_human_annotation

size=$(gzip -l $gz_human_seq  | grep -Eo '\b[[:digit:]]+\b ' | tail -n1)
gzip -dc $gz_human_seq | tqdm --desc "Extracting Human Genome Sequence" --bytes --total $size > $fna_human_seq

size=$(gzip -l $gz_human_annotation  | grep -Eo '\b[[:digit:]]+\b ' | tail -n1)
gzip -dc $gz_human_annotation | tqdm --desc "Extracting Human Genome Annotation" --bytes --total $size > $fna_human_annotation


echo
echo "${header}Finalising Dataset${headerend}"
echo

echo "Copying files to final location"
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
    echo "Deleting build files"
    rm -r data/genomes/build/
fi
