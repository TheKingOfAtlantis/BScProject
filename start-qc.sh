
echo "Gene vs Genome"
python scripts/QC-GenomeVsGenes.py
echo "Translation Table"
python scripts/QC-TranslationTable.py
echo "Proteins"
python scripts/QC-Protein.py
python scripts/QC-Protein-Isolate.py
echo "IDs"
python scripts/QC-AvaliableIDs.py
