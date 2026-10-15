#!/usr/bin/bash

# ---------------------------------
# rename_quast_with_assembly_level.sh
# Add the NCBI assembly level to the reference file names in QUASTDB:
#   Genus_species_genomic.fna/.gff  ->  Genus_species__AssemblyLevel_genomic.fna/.gff
#
# For each GFF:
#   1. read the accession (GCF_...) from the header line "#!genome-build-accession NCBI_Assembly:GCF_..."
#   2. look it up in assembly_summary.txt (column 1 = accession, column 12 = assembly_level)
#   3. rename the FASTA and the GFF with the level (spaces removed, e.g. "Complete Genome" -> CompleteGenome)
#
# Usage: bash rename_quast_with_assembly_level.sh [--dry-run] [assembly_summary.txt]
#   If assembly_summary.txt is missing it is downloaded (RefSeq bacteria): run on a node with internet access.
# Files already carrying a level are skipped, so the script can be re-run safely.
# A summary of what was done is appended to quast_reference_levels.tsv (current directory).
# ---------------------------------

# ---------------------------------
set -euo pipefail

QUASTDB="${QUASTDB:-/dlocal/home/2019013/Databases/QUAST/Reference_genomes}"

DRY=0
if [[ "${1:-}" == "--dry-run" ]]; then DRY=1; shift; fi
SUMMARY="${1:-assembly_summary_bacteria.txt}"

if [[ ! -s "${SUMMARY}" ]]; then
    echo "Downloading assembly_summary.txt (RefSeq bacteria)..."
    wget -q -O "${SUMMARY}" "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
fi

OUT="quast_reference_levels.tsv"
[[ -s "${OUT}" ]] || printf "species\taccession\tassembly_level\n" > "${OUT}"

shopt -s nullglob
for gff in "${QUASTDB}"/*_genomic.gff; do
    base="$(basename "${gff}" _genomic.gff)"          # Genus_species

    # Already renamed?
    if [[ "${base}" =~ __(CompleteGenome|Chromosome|Scaffold|Contig)$ ]]; then
        continue
    fi

    fna="${QUASTDB}/${base}_genomic.fna"
    [[ -s "${fna}" ]] || { echo "WARN ${base}: FASTA missing, skipped" >&2; continue; }

    acc="$(head -n 20 "${gff}" | grep -m1 '^#!genome-build-accession' | sed -E 's/.*:[[:space:]]*//' | tr -d '[:space:]' || true)"
    [[ -n "${acc}" ]] || { echo "WARN ${base}: no #!genome-build-accession in GFF header, skipped" >&2; continue; }

    level="$(awk -F'\t' -v a="${acc}" '$1==a {print $12; exit}' "${SUMMARY}" | tr -d ' ')"
    [[ -n "${level}" ]] || { echo "WARN ${base}: ${acc} not found in ${SUMMARY}, skipped" >&2; continue; }

    new="${base}__${level}"
    echo "${base}: ${acc} -> ${level}"
    printf "%s\t%s\t%s\n" "${base}" "${acc}" "${level}" >> "${OUT}"

    (( DRY )) && continue
    mv -n "${fna}" "${QUASTDB}/${new}_genomic.fna"
    mv -n "${gff}" "${QUASTDB}/${new}_genomic.gff"
done

echo "Done. Summary in ${OUT}"
# ---------------------------------