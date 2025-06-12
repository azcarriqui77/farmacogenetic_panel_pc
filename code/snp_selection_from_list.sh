# Script to select the IDs and genotype information from an ID list
#!/bin/bash

# Verify that the correct files are passed as arguments
if [ "$#" -ne 2 ]; then
    echo "Uso: $0 archivo_genotipos.tsv lista_ids.txt"
    exit 1
fi

GENOTIPOS="$1"
LISTA_IDS="$2"

# Use awk to filter the lines which are in the list
awk '
NR==FNR { 
    gsub("\r", "", $1);  # Delete posible \r from Windows
    ids[$1] = 1;         # Save IDs in a dictionary
    next;
}
FNR==1 || ($1 in ids)  # Print header and coinciding lines
' "$LISTA_IDS" "$GENOTIPOS"