#!/usr/bin/bash

# ---------------------------------
# build_reference_db.sh
# Build the reference genome databases for one genus in a single run:
#   1. download the RefSeq reference genomes of the genus (NCBI datasets) into PARSNPDB/<Genus>/
#        file names: <Genus_species>__<GCF_accession>_<assembly>_genomic.fna
#   2. for the species you choose, add FASTA + GFF to QUASTDB   (add_quast_reference.sh)
#   3. add the NCBI assembly level to the QUASTDB file names     (rename_quast_with_assembly_level.sh)
#        file names: <Genus_species>__<AssemblyLevel>_genomic.fna/.gff
#
# Usage:
#   bash build_reference_db.sh -g GENUS [-s sp1,sp2,...] [-l species_list.txt] [-A] [-a assembly_summary.txt]
#     -g  Genus (e.g. Pseudomonas)                                            [mandatory]
#     -s  Comma-separated Genus_species names to add to QUASTDB
#     -l  File with one Genus_species per line (lines starting with # ignored)
#     -A  Add ALL species of the genus (every genome in PARSNPDB/<Genus>/) to QUASTDB
#     -a  Local assembly_summary.txt (downloaded if missing)
#     -h  Show this help
#   -s, -l and -A can be combined. Without any of them only step 1 is run.
#
# Requirements: ncbi-datasets-cli (datasets, dataformat) ncbi-datasets_env, unzip, wget; internet access (login node).
# add_quast_reference.sh and rename_quast_with_assembly_level.sh must be in the same folder as this script.
# Safe to re-run: existing files are never overwritten.
# Example:
# bash build_reference_db.sh -g Pseudomonas                       # step 1 only (download the genomes)
# bash build_reference_db.sh -g Pseudomonas -s Pseudomonas_koreensis,Pseudomonas_putida_nbrc_14164
# bash build_reference_db.sh -g Pseudomonas -l especes_quast.txt  # one Genus_species per line
# bash build_reference_db.sh -g Citrobacter -A                    # all species of the genus into QUASTDB
# ---------------------------------

# ---------------------------------
set -euo pipefail

PARSNPDB="${PARSNPDB:-/dlocal/home/2019013/Databases/ParSNP}"
QUASTDB="${QUASTDB:-/dlocal/home/2019013/Databases/QUAST/Reference_genomes}"
export PARSNPDB QUASTDB
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Print the header block (between the first two "# ----" lines)
usage() { awk '/^# -+$/ {n++; if (n==2) exit; next} n==1 {sub(/^# ?/, ""); print}' "${BASH_SOURCE[0]}" >&2; exit "${1:-1}"; }

GENUS=""; SPLIST=""; SPFILE=""; ASUM=""; ALL=0
while getopts "g:s:l:a:Ah" opt; do
    case "${opt}" in
        g) GENUS="${OPTARG}" ;;
        s) SPLIST="${OPTARG}" ;;
        l) SPFILE="${OPTARG}" ;;
        a) ASUM="${OPTARG}" ;;
        A) ALL=1 ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done
[[ -n "${GENUS}" ]] || { echo "ERROR: -g GENUS is mandatory" >&2; usage 1; }

# Folder names are case-sensitive: always "Genus"
GENUS="${GENUS,,}"; GENUS="${GENUS^}"
# ---------------------------------

# ---------------------------------
# 1. Download genomes into PARSNPDB/<Genus>/
# ---------------------------------
OUTDIR="${PARSNPDB}/${GENUS}"
mkdir -p "${OUTDIR}"
WORK="$(mktemp -d -p "${PWD}" refdb_tmp.XXXXXX)"
trap 'rm -rf "${WORK}"' EXIT

echo "[1/3] Downloading RefSeq reference genomes for ${GENUS}..."
datasets download genome taxon "${GENUS}" \
    --reference --assembly-source RefSeq \
    --include genome --filename "${WORK}/ref.zip"
unzip -oq "${WORK}/ref.zip" -d "${WORK}/ref"
D="${WORK}/ref/ncbi_dataset/data"

# accession -> organism name
dataformat tsv genome --inputfile "${D}/assembly_data_report.jsonl" \
    --fields accession,organism-name --elide-header > "${WORK}/acc_names.tsv"

new=0
shopt -s nullglob
while IFS=$'\t' read -r acc name; do
    sp="$(echo "${name}" | tr '[:upper:]' '[:lower:]' | tr -c 'a-z0-9\n' '_' | sed 's/_*$//')"
    sp="${sp^}"                                   # Genus_species with a capital G
    for f in "${D}/${acc}"/*_genomic.fna; do
        dest="${OUTDIR}/${sp}__$(basename "${f}")"
        [[ -e "${dest}" ]] && continue
        cp "${f}" "${dest}"
        new=$((new + 1))
    done
done < "${WORK}/acc_names.tsv"
shopt -u nullglob

echo "    ${new} new genome(s) copied; ${OUTDIR} now holds $(ls "${OUTDIR}"/*.fna 2>/dev/null | wc -l) genome(s)"
# ---------------------------------

# ---------------------------------
ls "${OUTDIR}"/*fna | xargs -n 1 basename > "${GENUS}_Quast.txt"
# ---------------------------------

# ---------------------------------
# 2. Species to add to QUASTDB
# ---------------------------------
SPECIES=()
if [[ -n "${SPLIST}" ]]; then
    IFS=',' read -r -a SPECIES <<< "${SPLIST}"
fi
if [[ -n "${SPFILE}" ]]; then
    while IFS= read -r line; do
        line="${line//[[:space:]]/}"
        if [[ -n "${line}" && "${line}" != \#* ]]; then
            SPECIES+=("${line}")
        fi
    done < "${SPFILE}"
fi
if (( ALL )); then
    # Genus_species__GCF_...  ->  Genus_species
    shopt -s nullglob
    for f in "${OUTDIR}"/*__GCF_*_genomic.fna; do
        b="$(basename "${f}")"
        SPECIES+=("${b%%__GCF_*}")
    done
    shopt -u nullglob
fi

if [[ ${#SPECIES[@]} -eq 0 ]]; then
    echo "No species given (-s / -l / -A): QUASTDB not modified."
    exit 0
fi

# Remove duplicates
mapfile -t SPECIES < <(printf '%s\n' "${SPECIES[@]}" | sort -u)

echo "[2/3] Adding ${#SPECIES[@]} species to ${QUASTDB}..."
failed=()
for sp in "${SPECIES[@]}"; do
    bash "${HERE}/add_quast_reference.sh" "${sp}" || failed+=("${sp}")
done
# ---------------------------------

# ---------------------------------
# 3. Assembly level in the QUASTDB file names
# ---------------------------------
echo "[3/3] Adding the assembly level to the QUASTDB file names..."
bash "${HERE}/rename_quast_with_assembly_level.sh" ${ASUM:+"${ASUM}"}

if [[ ${#failed[@]} -gt 0 ]]; then
    echo "WARNING: ${#failed[@]} species could not be added: ${failed[*]}" >&2
    exit 1
fi
echo "Done."
# ---------------------------------
