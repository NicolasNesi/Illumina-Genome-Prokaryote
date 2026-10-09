#!/usr/bin/bash
#SBATCH --job-name=P01_DeNovoAssembly_Annotation              # Job name
#SBATCH --output=01_DeNovoAssembly_Annotation.%a.%A.out    # File to which stdout will be written
#SBATCH --error=P01_DeNovoAssembly_Annotation.%a.%A.err    # File to which stderr will be written
#SBATCH --partition=tcourt                                     # Partition
#SBATCH --cpus-per-task=48                        # Number of cores/cpus, 1 core has 2 threads (96 cores = 192 threads)
#SBATCH --time=01-00:00                                   # Runtime in DD-HH:MM
#SBATCH --mem=100G                                    # Memory for all cores in Gbytes (or --mem 752000 in Mbytes)
#SBATCH --mail-type=ALL                               # BEGIN,END,ALL
#SBATCH --mail-user=nicolas.nesi@unicaen.fr     # Email address

# ---------------------------------
# Version v2.1
# ---------------------------------
# author: Benedicte Langlois, Nicolas Nesi
# University Caen Normandy, INSERM UMR 1311 DYNAMICURE
# Date: 06/10/2026
# ---------------------------------

# ---------------------------------
# Submission:
# sbatch --array=1-x P01_DeNovoAssembly_Annotation_Template.sh 
# ---------------------------------

# ---------------------------------
# Variables and Paths
# ---------------------------------
# safety for failing steps
set -euo pipefail

# conda switch environments
conda_switch() {
    set +u
    conda deactivate 2>/dev/null || true
    conda activate "$1"
    set -u
}

# Positional parameters
GROUP="BacterioCaen"
PROJECT="FlorieGenomes"
RUN=""
QSCORE="20"
GENUS=""

usage() {
   local exit_code="${1:-1}"
   {
    echo "Usage: sbatch --array=1-x P01_DeNovoAssembly_Annotation_Template.sh -r RUN -g GENUS [-q QSCORE]"
	echo "   -r, --run    Mandatory. Name of run directory"
	echo "   -g, --genus  Mandatory. Name of the genus for Bakta annotation (e.g. Pseudomonas)"
	echo "   -q, --qscore Optional. Qscore threshold for reads quality filter [default 20]"
	echo "   -h, --help   Show this help"
    } >&2
	exit "$exit_code"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
    -r|--run) RUN="$2"; shift 2 ;;
	-q|--qscore) QSCORE="$2"; shift 2 ;;
	-g|--genus) GENUS="$2"; shift 2 ;;
	-h|--help)  usage 0 ;;
    *) 
         echo "Invalid option: $1" >&2
         usage 1
         ;;
    esac
done

if [[ -z "$RUN" || -z "$GENUS" ]]; then
    echo "Erreur: -r and -g are mandatory" >&2
    usage 1
fi

# Paths
FOLDER="/dlocal/home/2019013/Data/$GROUP/$PROJECT"
PATH1="${FOLDER}/Output_Illumina/${RUN}"
SCRIPT="/home/2019013/PARTAGE/Various_Scripts_for_various_purposes"

# Databases
KRAKENDB="/dlocal/home/2019013/Databases/Kraken2_PlusPFP-2026"
BAKTADB="/dlocal/home/2019013/Databases/Bakta/db"
QUASTDB="/dlocal/home/2019013/Databases/QUAST/Reference_genomes"
BUSCODB="/dlocal/home/2019013/Databases/BUSCO/lineages"
PARSNPDB"/dlocal/home/2019013/Databases/ParSNP"

# Array job
SAMPLE="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" "$FOLDER/Scripts/List_samples_$RUN.txt" | awk '{print $1}')"
REFGENOME="$(sed -n -e "${SLURM_ARRAY_TASK_ID}p" "$FOLDER/Scripts/List_samples_$RUN.txt" | awk '{print $2}')"

# Variables
RUNDATE="$(date '+%Y-%m-%d')"

# Threads
export THREADS="${SLURM_CPUS_PER_TASK}"
echo "number of threads used ${THREADS}"
# ---------------------------------

# ---------------------------------
# Copy input data to LOCAL_WORK_DIR
# ---------------------------------
cp "${PATH1}"/Fastq/"${SAMPLE}"_R*.fastq.gz "${LOCAL_WORK_DIR}"
cd "${LOCAL_WORK_DIR}"/
gunzip *.fastq.gz

# Generate setting file of the analysis:
printf "Analysis run the: $RUNDATE \n
Working directory: $PWD \n
Paths variables are: group $GROUP project $PROJECT and analyses of the run $RUN \n
Filtering using fastq with a minimum qscore of $QSCORE \n
De novo assembly using Unicycler \n
genomes completeness and contamination with CheckM and BUSCO \n
Annotation using Bakta for genus $GENUS \n"
# ---------------------------------

# ---------------------------------
# Environments
# ---------------------------------
module load py_env/miniconda3/25.7.0
set +u
eval "$(conda shell.bash hook)"
conda activate orthoproka_env
set -u
# Include:
# bowtie2=2.5.4
# checkm-genome=1.2.3
# fastp=0.23.4
# unicycler=0.5.1
# SPAdes=4.0.0
# samtools=1.20
# ---------------------------------

# ---------------------------------
# Short read QC
# ---------------------------------
echo "QC and filtering of short reads for sample: $SAMPLE"

fastp \
--in1 "${SAMPLE}_R1.fastq" \
--in2 "${SAMPLE}_R2.fastq" \
--detect_adapter_for_pe \
--length_required 50 \
--qualified_quality_phred "${QSCORE}" \
--low_complexity_filter \
--thread 16 \
--json "${SAMPLE}"_ReportFastp.json \
--html "${SAMPLE}"_ReportFastp.html \
--report_title "$SAMPLE fastq report" \
--out1 "${SAMPLE}_CleanedReads_R1.fastq" \
--out2 "${SAMPLE}_CleanedReads_R2.fastq" \
--unpaired1 "${SAMPLE}_CleanedReads_unpaired.fastq" \
--unpaired2 "${SAMPLE}_CleanedReads_unpaired.fastq" 2>&1 | tee -a "${SAMPLE}_fastp_run${RUNDATE}.log"
# ---------------------------------

# ---------------------------------
# De novo assembly from Illumina reads
# ---------------------------------
echo "De novo assembly with Unicycler for sample $SAMPLE"

unicycler \
--short1 "${SAMPLE}_CleanedReads_R1.fastq" \
--short2 "${SAMPLE}_CleanedReads_R2.fastq" \
--unpaired "${SAMPLE}_CleanedReads_unpaired.fastq" \
--out "${SAMPLE}"_unicycler_assembly/ \
--min_fasta_length 500 \
--depth_filter 0.25 \
--keep 2 \
--mode normal \
--threads "${THREADS}" 2>&1 | tee -a "${SAMPLE}_Unicycler_run${RUNDATE}.log"

mv "${SAMPLE}"_unicycler_assembly/assembly.fasta "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta"
mv "${SAMPLE}"_unicycler_assembly/assembly.gfa "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.gfa"
mv "${SAMPLE}"_unicycler_assembly/unicycler.log "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_unicycler.log"
# ---------------------------------

# ---------------------------------
# Quality average coverage
# ---------------------------------
echo "Mapping back reads against assembly for sample $SAMPLE"

bowtie2-build "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta" ./"${SAMPLE}"_assembly

bowtie2 \
-1 "${SAMPLE}_CleanedReads_R1.fastq" \
-2 "${SAMPLE}_CleanedReads_R2.fastq" \
-U "${SAMPLE}_CleanedReads_unpaired.fastq" \
-x ./"${SAMPLE}"_assembly \
--very-sensitive \
--dovetail \
--time \
--threads "${THREADS}" \
-S "${SAMPLE}_assembly.sam" 2>&1 | tee -a "${SAMPLE}_Bowtie2Assembly_run${RUNDATE}.log"

echo "Generate coverage statistics"

samtools view --threads 96 -S --bam "${SAMPLE}_assembly.sam" -o "${SAMPLE}_assembly.bam"
samtools sort --threads 96 "${SAMPLE}_assembly.bam" -o "${SAMPLE}_assembly.sorted.bam"
samtools index "${SAMPLE}_assembly.sorted.bam"
samtools depth "${SAMPLE}_assembly.sorted.bam" | \
awk '{sum+=$3; count++} END {if (count>0) print sum/count; else print "No data"}' > "${SAMPLE}"_MeanDepthCoverage.txt
# ---------------------------------

# ---------------------------------
# Taxonomic Profile
# ---------------------------------
echo "Starting taxonomic profile for sample $SAMPLE"

conda_switch MTX_env
# include:
# Kraken2=2.1.3
# KrakenTools=1.2
# Krona=2.8.1

kraken2 \
--db "${KRAKENDB}" \
--memory-mapping \
--threads "${THREADS}" \
--minimum-hit-groups 3 \
--output "${SAMPLE}_Kraken2.tsv" \
--report "${SAMPLE}.k2report" \
"${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta" 2>&1 | tee -a "${SAMPLE}_Kraken2.log"

# Kraken2 Krona
kreport2krona.py \
-r "${SAMPLE}.k2report" \
-o "${SAMPLE}_Kraken2_Krona.txt"

ktImportText "${SAMPLE}_Kraken2_Krona.txt" \
-o "${SAMPLE}_Kraken2_Krona.html"
# ---------------------------------

# ---------------------------------
# Plasmid contigs detection et suppression
# ---------------------------------
echo "Plasmids contigs detection for sample $SAMPLE"

conda_switch plasme_env
# include:
# PLASMe=1.1
# seqkit=2.9.0

sed -i.bak "s/>/>contig_/" "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta"

python /soft/2019013/Logiciels/PLASMe/PLASMe.py \
"${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta" \
"${SAMPLE}"_predicted_plasmids.fasta \
--mode balance \
--thread "${THREADS}" \
--database /soft/2019013/Logiciels/PLASMe/DB/ \
--temp "${SAMPLE}"_temp

awk -F '\t' '{print $1}' "${SAMPLE}"_predicted_plasmids.fasta_report.csv | sed '1d' > "${SAMPLE}"_list_contig_plasme.txt

seqkit grep -c -v -f "${SAMPLE}"_list_contig_plasme.txt "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta" > "${SAMPLE}"_unicycler_assembly/"${SAMPLE}"_assembly_chrom.fasta
# ---------------------------------

# ---------------------------------
# Genome annotation
# ---------------------------------
conda_switch bakta_env
# Include:
# bakta=1.12.1

echo "Annotation of chromosome for sample $SAMPLE"

bakta \
--db "${BAKTADB}" \
"${SAMPLE}_unicycler_assembly/${SAMPLE}_assembly_chrom.fasta" \
--prefix "${SAMPLE}"_chrom \
--genus "${GENUS}" \
--output "${SAMPLE}"_chrom_bakta/ 2>&1 | tee -a "${SAMPLE}_chrom_bakta_run${RUNDATE}.log"

python "${SCRIPT}"/convert_gbff_to_gbk.py -input "${SAMPLE}"_chrom_bakta/"${SAMPLE}"_chrom.gbff -output "${SAMPLE}"_chrom_bakta/"${SAMPLE}"_chrom.gbk

echo "Add pseudogene tag in sequences header"

awk -F'\t|;' '$0 ~ /pseudo=True/ {for (i=1; i<=NF; i++) if ($i ~ /^ID=/) print substr($i, 4)}' "${SAMPLE}"_chrom_bakta/"${SAMPLE}"_chrom.gff3 > "${SAMPLE}"_chrom_bakta/"${SAMPLE}"_pseudogene_ids.txt

if [[ ! -s "${SAMPLE}_chrom_bakta/${SAMPLE}_pseudogene_ids.txt" ]]; then
    echo "No chromosomic pseudogenes found by Bakta"
else
  awk 'BEGIN {while ((getline < "'"${SAMPLE}_chrom_bakta/${SAMPLE}_pseudogene_ids.txt"'") > 0) ids[$1] = 1; close("'"${SAMPLE}_chrom_bakta/${SAMPLE}_pseudogene_ids.txt"'")}{if ($0 ~ /^>/) {match($0, /^>([^ ]+)/, arr);if (arr[1] in ids) print $0 " [pseudogene]"; else print $0;} else { print $0;}}' "${SAMPLE}"_chrom_bakta/"${SAMPLE}"_chrom.ffn > "${SAMPLE}"_chrom_bakta/"${SAMPLE}"_chrom_edited.ffn
fi

if [[ -s "${SAMPLE}"_predicted_plasmids.fasta  ]]; then
   echo "Annotation of plasmids for sample $SAMPLE"
   
   bakta \
   --db "${BAKTADB}" \
   "${SAMPLE}"_predicted_plasmids.fasta \
   --prefix "${SAMPLE}"_plasmid \
   --genus "${GENUS}" \
   --output "${SAMPLE}"_plasmid_bakta/ 2>&1 | tee -a "${SAMPLE}_plasmid_bakta_run${RUNDATE}.log"
   
   python "${SCRIPT}"/convert_gbff_to_gbk.py -input "${SAMPLE}"_plasmid_bakta/"${SAMPLE}"_plasmid.gbff -output "${SAMPLE}"_plasmid_bakta/"${SAMPLE}"_plasmid.gbk
   
   echo "Add pseudogene tag in sequences header"
   
   awk -F'\t|;' '$0 ~ /pseudo=True/ {for (i=1; i<=NF; i++) if ($i ~ /^ID=/) print substr($i, 4)}' "${SAMPLE}"_plasmid_bakta/"${SAMPLE}"_plasmid.gff3 > "${SAMPLE}"_plasmid_bakta/"${SAMPLE}"_pseudogene_ids.txt
   
   if [[ ! -s "${SAMPLE}_chrom_bakta/${SAMPLE}_pseudogene_ids.txt" ]]; then
       echo "No plasmid pseuodgenes found by Bakta"
   else
     awk 'BEGIN {while ((getline < "'"${SAMPLE}_plasmid_bakta/${SAMPLE}_pseudogene_ids.txt"'") > 0) ids[$1] = 1; close("'"${SAMPLE}_plasmid_bakta/${SAMPLE}_pseudogene_ids.txt"'")}{if ($0 ~ /^>/) {match($0, /^>([^ ]+)/, arr);if (arr[1] in ids) print $0 " [pseudogene]"; else print $0;} else { print $0;}}' "${SAMPLE}"_plasmid_bakta/"${SAMPLE}"_plasmid.ffn > "${SAMPLE}"_plasmid_bakta/"${SAMPLE}"_plasmid_edited.ffn
   fi
else
  echo "No plasmids found"
fi
# ---------------------------------

# ---------------------------------
# Genome completeness contamination
# ---------------------------------
echo "Genome completeness and contamination for sample $SAMPLE"

conda_switch orthoproka_env

checkm lineage_wf \
--extension fasta \
"${SAMPLE}"_unicycler_assembly/ \
--tab_table \
--threads "${THREADS}" \
"${SAMPLE}"_checkm/

conda_switch BUSCO_env
# Include:
# busco=5.8.3

busco \
--in "${SAMPLE}_unicycler_assembly/${SAMPLE}_assembly_chrom.fasta" \
--mode genome \
--lineage_dataset "${BUSCODB}/enterobacterales_odb12/" \
--out "${SAMPLE}"_BUSCO/ \
--cpu "${THREADS}"

python3 /soft/2019013/Conda_env/BUSCO_env/bin/generate_plot.py --working_directory "${SAMPLE}"_BUSCO/
mv "${SAMPLE}"_BUSCO/busco_figure.png "${SAMPLE}"_BUSCO/"${SAMPLE}"_busco_figure.png

conda_switch denovo_env
# Include:
# quast=5.3.0

echo "Running QUAST for sample $SAMPLE"

mkdir -p "QUAST_${SAMPLE}"

quast \
"${SAMPLE}_unicycler_assembly/${SAMPLE}_assembly_chrom.fasta" \
--output-dir "QUAST_${SAMPLE}" \
--pe1 "${SAMPLE}_CleanedReads_R1.fastq" \
--pe2 "${SAMPLE}_CleanedReads_R1.fastq" \
-r "${QUASTDB}/${REFGENOME}_genomic.fna" \
--features gene:"${QUASTDB}/${REFERENCE}_genomic.gff" \
--gene-finding \
--plots-format svg \
--threads "${THREADS}" 2>&1 | tee -a "${FOLDER}/Log/${RUN}_run${RUNDATE}/${SAMPLE}_QUAST_run${RUNDATE}.log"
# ---------------------------------

# ---------------------------------
# Phylogeny Parsnp, fastANI and IQtree
# ---------------------------------
echo "phylogeny using Parsnp, fastANI and IQtree for sample $SAMPLE"

conda_switch parsnp_env
# Include:
# fastANI=1.34
# parsnp=2.1.6

mkdir -p genomes_parsnp

cp "${SAMPLE}_unicycler_assembly/${SAMPLE}_assembly_chrom.fasta" genomes_parsnp/
cp "${PARSNPDB}"/Pseudomonas/*.fna genomes_parsnp/

parsnp \
--reference "${PARSNPDB}/Pseudomonas/Pseudomonas_koreensis__GCF_001654515.1_ASM165451v1_genomic.fna" \
--sequences genomes_parsnp/ \
--threads "${THREADS}" \
--output-dir parsnp_out

harvesttools -x parsnp_out/parsnp.xmfa -M parsnp_out/core.aln.fa
harvesttools -x parsnp_out/parsnp.xmfa -S parsnp_out/core.snps.fa

fastANI \
--query "${SAMPLE}_unicycler_assembly/${SAMPLE}_assembly_chrom.fasta" \
--refList <(ls "${PARSNPDB}"/Pseudomonas/*.fna) \
--output "${SAMPLE}_fastANI.tsv" \
--threads "${THREADS}"

best=$(sort -k3,3nr ${SAMPLE}_fastANI.tsv | head -n1)
echo "Closest genome to $SAMPLE based on ANI: $(echo "$best" | cut -f2) (ANI = $(echo "$best" | cut -f3)%)"

conda_switch phylo_env
# Include:
# iqtree=2.4.0

iqtree2 \
-s parsnp_out/core.aln.fa \
-m MFP \
-B 1000 \
-alrt 1000 \
-T AUTO \
--prefix parsnp_out/iqtree_core

iqtree2 \
-s parsnp_out/core.snps.fa \
-m MFP+ASC \
-B 1000 \
-alrt 1000 \
-T AUTO \
--prefix parsnp_out/iqtree_snps
# ---------------------------------

# ---------------------------------
# Antibiotics and Secondary Metabolite
# ---------------------------------
echo "antiSMASH for sample $SAMPLE"

conda_switch antismash_env
# Include:
# antismash=8.0.4

antismash \
"${SAMPLE}"_chrom_bakta/"${SAMPLE}_chrom.gbk" \
--databases /dlocal/home/2019013/Databases/antismash \
--taxon bacteria \
--output-dir "antismash_${SAMPLE}" \
--output-basename "${SAMPLE}" \
--fullhmmer \
--cc-mibig \
--asf \
--rre \
--tfbs \
--cb-knownclusters \
--cb-subclusters \
--genefinding-tool none \
--html-title "${SAMPLE}" \
--html-start-compact \
--cpus "${THREADS}"
# ---------------------------------

# ---------------------------------
# Output directories and Move data
# ---------------------------------
if [[ ! -d "${FOLDER}/CleanedReads_Illumina" ]]; then
   mkdir -p "${FOLDER}/CleanedReads_Illumina"
fi

mv "${SAMPLE}"_CleanedReads_*.fastq "${FOLDER}/CleanedReads_Illumina/"

conda_switch NN_env
pigz --best "${FOLDER}"/CleanedReads_Illumina/"${SAMPLE}"_CleanedReads_*.fastq

mv "${SAMPLE}"_ReportFastp.html "${FOLDER}/CleanedReads_Illumina/"

if [[ ! -d "${FOLDER}/Kraken2" ]]; then
   mkdir -p "${FOLDER}/Kraken2"
fi

mv "${SAMPLE}_Kraken2.tsv" "${FOLDER}/Kraken2"
mv "${SAMPLE}.k2report" "${FOLDER}/Kraken2"
mv "${SAMPLE}_Kraken2_Krona.html" "${FOLDER}/Kraken2"
mv "${SAMPLE}_Kraken2.log" "${FOLDER}/Kraken2"

if [[ ! -d "${FOLDER}/Assembly_Unicycler" ]]; then
   mkdir -p "${FOLDER}/Assembly_Unicycler"
fi

mv "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta" "${FOLDER}/Assembly_Unicycler"
mv "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.gfa" "${FOLDER}/Assembly_Unicycler"
mv "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_unicycler.log" "${FOLDER}/Assembly_Unicycler"
mv "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly_chrom.fasta" "${FOLDER}/Assembly_Unicycler"
mv "${SAMPLE}"_unicycler_assembly/"${SAMPLE}_assembly.fasta.bak" "${FOLDER}/Assembly_Unicycler"

if [[ ! -d "${FOLDER}/PLASMe" ]]; then
   mkdir -p "${FOLDER}/PLASMe"
fi

mv "${SAMPLE}_predicted_plasmids.fasta" "${FOLDER}/PLASMe/"
mv "${SAMPLE}_predicted_plasmids.fasta_report.csv" "${FOLDER}/PLASMe/"

if [[ ! -d "${FOLDER}/CheckM" ]]; then
   mkdir -p "${FOLDER}/CheckM"
fi

mv "${SAMPLE}"_checkm/storage/bin_stats_ext.tsv "${FOLDER}"/CheckM/"${SAMPLE}"_qa_statistics.txt
mv "${SAMPLE}"_checkm/storage/bin_stats.analyze.tsv "${FOLDER}"/CheckM/"${SAMPLE}"_general_statistics.txt

if [[ ! -d "${FOLDER}/BUSCO" ]]; then
   mkdir -p "${FOLDER}/BUSCO"
fi

mv "${SAMPLE}"_BUSCO/short_summary*.txt "${FOLDER}"/BUSCO/"${SAMPLE}"_short_summary.txt
mv "${SAMPLE}"_BUSCO/"${SAMPLE}"_busco_figure.png "${FOLDER}"/BUSCO/"${SAMPLE}"_busco_figure.png

if [[ ! -d "${FOLDER}/QUAST" ]]; then
   mkdir -p "${FOLDER}/QUAST"
fi

cp -R "QUAST_${SAMPLE}" "${FOLDER}/QUAST"

if [[ ! -d "${FOLDER}/Bakta" ]]; then
   mkdir -p "${FOLDER}/Bakta"
fi

if [[ -d "${SAMPLE}"_chrom_bakta ]]; then
   mv "${SAMPLE}"_chrom_bakta "${FOLDER}/Bakta"
fi

if [[ -d "${SAMPLE}"_plasmid_bakta ]]; then
    mv "${SAMPLE}"_plasmid_bakta "${FOLDER}/Bakta"
fi

if [[ ! -d "${FOLDER}/Bowtie2" ]]; then
   mkdir -p "${FOLDER}/Bowtie2"
fi

mv "${SAMPLE}"_assembly.sorted.ba* "${FOLDER}/Bowtie2"
mv "${SAMPLE}_MeanDepthCoverage.txt" "${FOLDER}/Bowtie2"
mv "${SAMPLE}_Bowtie2Assembly_run${RUNDATE}.log" "${FOLDER}/Bowtie2"

if [[ ! -d "${FOLDER}/Antismash" ]]; then
   mkdir -p "${FOLDER}/Antismash"
fi

cp -R "antismash_${SAMPLE}" "${FOLDER}/Antismash/"

if [[ ! -d "${FOLDER}/Parsnp" ]]; then
   mkdir -p "${FOLDER}/Parsnp"
fi

mv parsnp_out/parsnp.xmfa "${FOLDER}/Parsnp/${SAMPLE}_parsnp.xmfa"
mv parsnp_out/parsnp.ggr "${FOLDER}/Parsnp/${SAMPLE}_parsnp.ggr"
mv parsnp_out/core.aln.fa "${FOLDER}/Parsnp/${SAMPLE}_core.aln.fa"
mv parsnp_out/core.snps.fa "${FOLDER}/Parsnp/${SAMPLE}_core.snps.fa"
mv parsnp_out/parsnp.tree "${FOLDER}/Parsnp/${SAMPLE}_parsnp.tree"

if [[ ! -d "${FOLDER}/Iqtree" ]]; then
   mkdir -p "${FOLDER}/Iqtree"
fi

mv parsnp_out/iqtree-core.treefile "${FOLDER}/Iqtree/${SAMPLE}_iqtree-core.treefile"
mv parsnp_out/iqtree-core.iqtree "${FOLDER}/Iqtree/${SAMPLE}_iqtree-core.iqtree"
mv parsnp_out/iqtree-snps.treefile "${FOLDER}/Iqtree/${SAMPLE}_iqtree-snps.treefile"
mv parsnp_out/iqtree-snps.iqtree "${FOLDER}/Iqtree/${SAMPLE}_iqtree-snps.iqtree"
# ---------------------------------

# ---------------------------------
# Graph and Statistic assembly
# ---------------------------------
echo "Calculate statistics and make an image of Unicycler assembly for sample $SAMPLE"

conda_switch bandage_env
# include:
# bandage=0.9.0

cd "${FOLDER}/Assembly_Unicycler"

Bandage info "${SAMPLE}_assembly.gfa" > "${SAMPLE}_Bandage_statistics.txt"
grep -E 'Node count:|Edge count:|Smallest edge overlap \(bp\):|Largest edge overlap \(bp\):|Total length \(bp\):|Total length no overlaps \(bp\):|Dead ends:|Percentage dead ends:|Connected components:|Largest component \(bp\):|Total length orphaned nodes \(bp\):|N50 \(bp\):|Shortest node \(bp\):|Lower quartile node \(bp\):|Median node \(bp\):|Upper quartile node \(bp\):|Longest node \(bp\):|Median depth:|Estimated sequence length \(bp\):' "${SAMPLE}_Bandage_statistics.txt" | awk '{print $NF}' > "${SAMPLE}_Bandage_statistics.temp"

Bandage image "${SAMPLE}_assembly.gfa" "${SAMPLE}_assembly.svg"
# ---------------------------------

# ---------------------------------
# Create statistics report QC genome
# ---------------------------------
echo "Generate a QC report of assembly for sample $SAMPLE"

REPORT="${FOLDER}/QCStats/${SAMPLE}_QC_Assembly_Report.tsv"
CHECKM="${FOLDER}/CheckM/${SAMPLE}_qa_statistics.txt"
BANDAGE="${FOLDER}/Assembly_Unicycler/${SAMPLE}_Bandage_statistics.txt"
DEPTH="${FOLDER}/Bowtie2/${SAMPLE}_MeanDepthCoverage.txt"
DEC_COMMA=1   # 1 = virgule décimale (Excel en français), 0 = point

mkdir -p "${FOLDER}/QCStats"

TMP="$(mktemp)"
trap 'rm -f "$TMP"' EXIT

for f in "$CHECKM" "$BANDAGE" "$DEPTH"; do
  [ -f "$f" ] || { echo "ERROR Missing files $f" >&2; exit 1; }
done

awk -F'\t' '$1 ~ /_chrom$/' "$CHECKM" \
  | grep -oE "'(# contigs|Completeness|Contamination|N50 \(contigs\)|Genome size)': [0-9.]+" \
  | sed -E "s/'//g; s/: /\t/;
            s/^# contigs/chrom_n_contigs/;
            s/^Completeness/chrom_completeness/;
            s/^Contamination/chrom_contamination/;
            s/^N50 \(contigs\)/chrom_N50/;
            s/^Genome size/chrom_genome_size/" > "$TMP" || true

[ -s "$TMP" ] || { echo "ERROR: no line _chrom in CheckM for ${SAMPLE}" >&2; exit 1; }

grep -E "^(Dead ends|Percentage dead ends):" "$BANDAGE" \
  | sed -E 's/%//; s/:[[:space:]]+/\t/;
            s/^Dead ends/assembly_dead_ends/;
            s/^Percentage dead ends/assembly_dead_ends_pct/' >> "$TMP" || true

printf "assembly_mean_depth\t%s\n" "$(tr -d '[:space:]' < "$DEPTH")" >> "$TMP"

[ "$(wc -l < "$TMP")" -eq 8 ] || { echo "ERROR: Incomplete repport for ${SAMPLE} ($(wc -l < "$TMP")/8 lines)" >&2; exit 1; }

LC_ALL=C awk -F'\t' -v comma="$DEC_COMMA" 'BEGIN{OFS="\t"}
  { if ($2 ~ /\./) { $2 = sprintf("%.2f", $2); if (comma) sub(/\./, ",", $2) } print }' \
  "$TMP" > "$REPORT"
  
# Klebsiella genome OK: - contigs <300
#                       - dead ends <200
#                       - N50 >50,000
#                       - average coverage >30x
#                       - genome size 4,969,898 to 6,132,846 bp
#                       - completeness >95%
#                       - contamination <1.5%
# ---------------------------------

# ---------------------------------
# Cleaning Working Directory
# ---------------------------------
ls -R > "${FOLDER}/Scripts/List_files_remaining_in_LocalWorkingDir_for_${SAMPLE}_run${RUNDATE}.txt"
rm -R "${LOCAL_WORK_DIR}"/*
# ---------------------------------
