# P01 – Bacterial De Novo Assembly & Annotation (Illumina, SLURM)

A SLURM array-job script that takes paired-end Illumina reads from a bacterial isolate and runs a complete pipeline: read QC, de novo assembly, reference-based consensus, taxonomic profiling, plasmid/chromosome separation, annotation, quality assessment, phylogeny and secondary-metabolite prediction. It processes **one sample per array task**.

- **Version:** 2.3
- **Authors:** Benedicte Langlois, Nicolas Nesi
- **Affiliation:** University of Caen Normandy, INSERM UMR 1311 DYNAMICURE
- **Target environment:** CRIANN HPC cluster (SLURM + conda/mamba)

---

## Pipeline overview

| # | Step | Tool(s) | Main output |
|---|------|---------|-------------|
| 1 | Short-read QC and filtering | `fastp` | Cleaned R1/R2/unpaired FASTQ, HTML/JSON report |
| 2 | De novo assembly | `Unicycler` (normal mode) | `<SAMPLE>_assembly.fasta` / `.gfa` |
| 3 | Mapping back to assembly, mean depth | `bowtie2`, `samtools` | Sorted BAM, mean depth |
| 4 | Reference-based consensus (min. depth 10x) | `bowtie2`, `bcftools`, `bedtools` | `_consensus.fa`, `_calls.vcf.gz`, `_low_cov.bed` |
| 5 | Taxonomic profiling of the assembly | `Kraken2`, `KrakenTools`, `Krona` | Kraken2 report, Krona HTML |
| 6 | Plasmid contig detection and removal | `PLASMe`, `seqkit` | `_predicted_plasmids.fasta`, `_assembly_chrom.fasta` |
| 7 | Annotation (chromosome and plasmids separately) | `Bakta` | GFF3/GBFF/GBK/FFN (pseudogenes tagged in FFN headers) |
| 8 | Completeness / contamination (BUSCO lineage chosen from `--genus`) | `CheckM`, `BUSCO` | CheckM tables, BUSCO summary and figure |
| 9 | Assembly metrics vs. reference | `QUAST` | QUAST report |
| 10 | Phylogeny and closest genome | `Parsnp`, `fastANI`, `IQ-TREE 2` | Core alignment/SNP trees, best ANI hit |
| 11 | Secondary metabolites | `antiSMASH` | antiSMASH HTML report |
| 12 | Assembly graph stats and image | `Bandage` | Bandage statistics, SVG graph |
| 13 | Consolidated QC report | bash/awk | `<SAMPLE>_QC_Assembly_Report.tsv` |

Intermediate files are written to `$LOCAL_WORK_DIR` (node-local scratch), results are moved to the project folder, and the scratch directory is cleaned at the end.

---

## Requirements

### Conda environments

Each step activates its own environment (switched with the `conda_switch` helper):

| Environment | Contents |
|-------------|----------|
| `orthoproka_env` | bcftools 1.21, bedtools 2.31.1, bowtie2 2.5.4, checkm-genome 1.2.3, fastp 0.23.4, unicycler 0.5.1, SPAdes 4.0.0, samtools 1.20 |
| `MTX_env` | Kraken2 2.1.3, KrakenTools 1.2, Krona 2.8.1 |
| `plasme_env` | PLASMe 1.1, seqkit 2.9.0 |
| `bakta_env` | bakta 1.12.2 |
| `BUSCO_env` | busco 5.8.3 |
| `denovo_env` | quast 5.3.0 |
| `parsnp_env` | fastANI 1.34, parsnp 2.1.6 (+ harvesttools) |
| `phylo_env` | iqtree 2.4.0 |
| `antismash_env` | antismash 8.0.4 |
| `bandage_env` | bandage 0.9.0 |
| `NN_env` | pigz |

### Databases

Set in the *Variables and Paths* section of the script:

| Variable | Content |
|----------|---------|
| `KRAKENDB` | Kraken2 PlusPFP database |
| `BAKTADB` | Bakta database |
| `QUASTDB` | Reference genomes `<Genus_species>__<AssemblyLevel>_genomic.fna` and `.gff` (QUAST and Bowtie2 consensus) |
| `BUSCODB` | BUSCO lineage datasets |
| `PARSNPDB` | Genomes for Parsnp / fastANI, one folder per genus: `<GENUS>/<Genus_species>__<GCF_accession>_<assembly>_genomic.fna` |
| `ANTISMASHDB` | antiSMASH databases |

PLASMe is installed outside conda (`/soft/.../PLASMe`) with its own database.

### Helper script

`convert_gbff_to_gbk.py` (in `$SCRIPT`) converts Bakta GBFF output to GBK for antiSMASH.

---

## Input

### Directory layout

```
/dlocal/home/<account>/Data/<GROUP>/<PROJECT>/
├── Output_Illumina/<RUN>/Fastq/<SAMPLE>_R1.fastq.gz
│                              /<SAMPLE>_R2.fastq.gz
├── Scripts/List_samples_<RUN>.txt
└── Log/<RUN>_run<DATE>/
```

### Sample list

`Scripts/List_samples_<RUN>.txt` is a tab-separated file, **one line per sample** (the line number is the SLURM array index):

```
<SAMPLE>    <REFGENOME>
```

- `SAMPLE`: sample name (prefix of the FASTQ files)
- `REFGENOME`: reference species as `Genus_species` (e.g. `Pseudomonas_koreensis`). Subspecies/pathovars use the full prefix (e.g. `Pseudomonas_moraviensis_subsp__stanleyae`). It must exist in both `QUASTDB` and `PARSNPDB/<GENUS>/` (see *Reference genomes*)

## Reference genomes

`REFGENOME` (column 2 of the sample list, `Genus_species`) is used in two databases with two naming conventions:

| Database | Used by | File naming |
|----------|---------|-------------|
| `QUASTDB` | Bowtie2 consensus, bcftools, QUAST | `<Genus_species>__<AssemblyLevel>_genomic.fna` and `.gff`, with `AssemblyLevel` = `CompleteGenome`, `Chromosome`, `Scaffold` or `Contig` |
| `PARSNPDB/<GENUS>/` | Parsnp reference, Parsnp tree, fastANI | `<Genus_species>__<GCF_accession>_<assembly>_genomic.fna` |

The Parsnp reference is found with the pattern `<Genus_species>__*_genomic.fna`; the script stops if it matches zero or several files. The double underscore means `Pseudomonas_syringae` does not match pathovar entries such as `Pseudomonas_syringae_pv__tomato_...`.

To add a species to `QUASTDB`, use the helper script (run on a node with internet access). It copies the FASTA from `PARSNPDB` under the QUAST name and downloads the GFF of the same assembly from NCBI, so sequence IDs match:

```bash
bash add_quast_reference.sh Pseudomonas_koreensis
```

The NCBI assembly level is then added to the file names with `rename_quast_with_assembly_level.sh`. For each GFF it reads the accession from the `#!genome-build-accession` header line, looks it up in NCBI's `assembly_summary.txt` (column 12, `assembly_level`) and renames the FASTA and GFF to `<Genus_species>__<AssemblyLevel>_genomic.fna/.gff`. Use `--dry-run` to preview. The script can be re-run safely and writes a summary to `quast_reference_levels.tsv`.

```bash
bash rename_quast_with_assembly_level.sh --dry-run   # preview
bash rename_quast_with_assembly_level.sh             # rename
```

The pipeline finds the reference in `QUASTDB` with the pattern `<Genus_species>__*_genomic.fna` and logs its assembly level. Prefer `CompleteGenome` references; `Scaffold` and `Contig` references are fragmented and give less reliable QUAST metrics and consensus.

---

## Usage

```bash
sbatch --array=1-<N> P01_DeNovoAssembly_Annotation_Template.sh -r <RUN> -g <GENUS> [-q <QSCORE>]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-r`, `--run` | **Mandatory.** Run directory name | – |
| `-g`, `--genus` | **Mandatory.** Genus used by Bakta (e.g. `Pseudomonas`) | – |
| `-q`, `--qscore` | Phred threshold for fastp read filtering | `20` |
| `-h`, `--help` | Show help | – |

Example for a 12-sample run:

```bash
sbatch --array=1-12 P01_DeNovoAssembly_Annotation_Template.sh -r Run_2026_10 -g Pseudomonas
```

### Hard-coded settings to adapt

| Setting | Value | Where |
|---------|-------|-------|
| `GROUP`, `PROJECT` | `BacterioCaen`, `FlorieGenomes` | Variables section |
| `MINDEPTH` | `10` | Minimum depth for consensus |
| SLURM resources | 48 CPUs, 100 GB, 24 h, partition `tcourt` | `#SBATCH` header |
| BUSCO lineage | Chosen from `--genus` in a `case` block: *Pseudomonas* → `pseudomonadales_odb12`; *Klebsiella, Citrobacter, Enterobacter, Escherichia, Cronobacter, Salmonella* → `enterobacterales_odb12`; *Enterococcus* → `lactobacillales_odb12`; otherwise `bacteria_odb12`. The script stops if the lineage folder is missing in `BUSCODB` | BUSCO step |
| Parsnp / fastANI genome set | `${PARSNPDB}/<GENUS>/*.fna`. `--genus` must match the folder name exactly, including case | Phylogeny step |

---

## Output

Results are collected in the project folder:

| Directory | Content |
|-----------|---------|
| `CleanedReads_Illumina/` | Filtered reads (pigz-compressed) and fastp HTML report |
| `Assembly_Unicycler/` | Assembly (full and chromosome-only FASTA), GFA, Unicycler log, Bandage stats and SVG |
| `Kraken2/` | Kraken2 output, report, Krona HTML |
| `PLASMe/` | Predicted plasmids and PLASMe report |
| `Bakta/` | Chromosome and plasmid annotations |
| `CheckM/`, `BUSCO/`, `QUAST/` | Quality assessment |
| `Bowtie2/` | Mapping files, mean depth, consensus, VCF, low-coverage BED |
| `Parsnp/`, `Iqtree/` | Core alignment, SNP alignment and phylogenetic trees |
| `Antismash/` | antiSMASH results |
| `QCStats/` | `<SAMPLE>_QC_Assembly_Report.tsv` |

### QC report

`<SAMPLE>_QC_Assembly_Report.tsv` summarises 8 metrics per sample: number of contigs, completeness, contamination, N50, genome size (from CheckM, chromosome bin), dead ends and % dead ends (Bandage), and mean depth. Decimal separator is a comma by default for French Excel (`DEC_COMMA=1`; set to `0` for a dot).

Example acceptance thresholds (defined for *Klebsiella*, adapt per species):

| Metric | Threshold |
|--------|-----------|
| Contigs | < 300 |
| Dead ends | < 200 |
| N50 | > 50,000 bp |
| Mean coverage | > 30x |
| Genome size | 4,969,898 – 6,132,846 bp |
| Completeness | > 95 % |
| Contamination | < 1.5 % |

---

## Notes

- The script uses `set -euo pipefail`: any failing step stops the task for that sample.
- `LOCAL_WORK_DIR` must be provided by the cluster environment.
- Contigs predicted as plasmids by PLASMe are removed from the assembly before chromosome annotation, CheckM bin selection, BUSCO and QUAST.
- If PLASMe predicts no plasmid, all contigs are treated as chromosome (`_assembly_chrom.fasta` is a copy of the assembly) and no plasmid annotation is run.
- Pseudogenes flagged by Bakta are tagged `[pseudogene]` in the FFN headers.

## Citation

If you use this pipeline, please cite the underlying tools (fastp, Unicycler, Bowtie2, SAMtools/BCFtools, Kraken2, PLASMe, Bakta, CheckM, BUSCO, QUAST, Parsnp, fastANI, IQ-TREE 2, antiSMASH, Bandage).

## Contact

Nicolas Nesi, INSERM UMR 1311 DYNAMICURE, University of Caen Normandy
