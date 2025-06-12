# Script which extract the ID (SNP) and Genotype information from a multisample VCF and
# saves it in a .tsv file
#!/bin/bash

# Verify that there is a real file as argument
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 archivo.vcf"
    exit 1
fi

VCF="$1"

# Extract the header, which contains the column names
header=$(grep -m1 "^#CHROM" "$VCF")

# Entract the ID and each sample's genotype columns
awk -v header="$header" '
BEGIN {
    FS="\t"; OFS="\t";
}
# Print personalized header
/^#CHROM/ {
    split(header, cols, "\t");
    printf "ID\t";
    for (i = 10; i <= NF; i++) {
        printf "%s%s", cols[i], (i < NF ? "\t" : "\n");
    }
    next;
}
# Extract ID and genotypes
!/^#/ {
    printf "%s\t", $3;  # Columna ID
    for (i = 10; i <= NF; i++) {
        split($i, gt, ":");
        printf "%s%s", gt[1], (i < NF ? "\t" : "\n");  # Solo el campo de genotipo (GT)
    }
}
' "$VCF"