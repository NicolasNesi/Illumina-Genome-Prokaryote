#!/usr/bin/bash

# ---------------------------------
# add_quast_reference.sh
# Add a reference genome (FASTA + GFF) to QUASTDB, from the Parsnp genome folder.
#
# Usage: bash add_quast_reference.sh Genus_species
#   e.g. bash add_quast_reference.sh Pseudomonas_koreensis
#        bash add_quast_reference.sh Pseudomonas_moraviensis_subsp__stanleyae
# A full Parsnp file name is also accepted (the __GCF_... suffix is stripped):
#        bash add_quast_reference.sh Pseudomonas_koreensis__GCF_001654515.1_ASM165451v1_genomic.fna
#
# - Reads  PARSNPDB/<Genus>/<Genus_species>__<GCF_accession>_<assembly>_genomic.fna
# - Writes QUASTDB/<Genus_species>_genomic.fna  (renamed, no GCF code)
# - Downloads the matching GFF from NCBI (same assembly, so sequence IDs match the FASTA)
# - Skips species already present in QUASTDB; nothing is written if the GFF download fails
# - Afterwards run rename_quast_with_assembly_level.sh to add the assembly level to the new file names
# Run on a node with internet access (login node).
# ---------------------------------

# ---------------------------------
set -euo pipefail

PARSNPDB="${PARSNPDB:-/dlocal/home/2019013/Databases/ParSNP}"
QUASTDB="${QUASTDB:-/dlocal/home/2019013/Databases/QUAST/Reference_genomes}"

SPECIES="${1:?Usage: $0 Genus_species}"
SPECIES="$(basename "${SPECIES}")"
SPECIES="${SPECIES%.fna}"
SPECIES="${SPECIES%__GCF_*}"
GENUS="${SPECIES%%_*}"

shopt -s nullglob
present=( "${QUASTDB}/${SPECIES}"__*_genomic.fna )      # renamed files (with assembly level)
shopt -u nullglob
if [[ -s "${QUASTDB}/${SPECIES}_genomic.fna" || ${#present[@]} -gt 0 ]]; then
    echo "SKIP: ${SPECIES} already in QUASTDB"
    exit 0
fi
shopt -s nullglob
refs=( "${PARSNPDB}/${GENUS}/${SPECIES}"__*_genomic.fna )
shopt -u nullglob
[[ ${#refs[@]} -eq 1 ]] || { echo "ERROR: expected 1 file for ${SPECIES}, found ${#refs[@]}" >&2; exit 1; }

base="$(basename "${refs[0]}" _genomic.fna)"
asm="${base##*__}"                      # e.g. GCF_001654515.1_ASM165451v1
path="$(sed -E 's/^GCF_([0-9]{3})([0-9]{3})([0-9]{3})\..*/\1\/\2\/\3/' <<< "${asm}")"
url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${path}/${asm}/${asm}_genomic.gff.gz"

tmp="$(mktemp)"
trap 'rm -f "${tmp}"' EXIT
if ! wget -q -O "${tmp}" "${url}"; then
    echo "ERROR: GFF download failed for ${SPECIES} (${url})" >&2
    exit 1
fi

gunzip -c "${tmp}" > "${QUASTDB}/${SPECIES}_genomic.gff"
cp "${refs[0]}" "${QUASTDB}/${SPECIES}_genomic.fna"

echo "OK: ${QUASTDB}/${SPECIES}_genomic.{fna,gff}"
# ---------------------------------