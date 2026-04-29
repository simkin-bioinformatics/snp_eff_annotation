# snp_eff_annotation

A Snakemake workflow for annotating variant calls with snpEff and extracting per-sample amino acid allele depth tables for targeted mutations. Designed for *Plasmodium falciparum* drug resistance surveillance but works with any organism given the appropriate genome files.

## Overview

The workflow:
1. Builds a custom snpEff database from a genome FASTA and GFF annotation
2. Annotates a merged VCF with snpEff
3. Parses the annotated VCF into four per-sample CSV tables, one each for coverage, reference allele depth, alternate allele depth, and other allele depth

Targeted mutations (defined in `config/targets.tsv`) are always included in the output tables even if they are not called in the VCF. The example targets cover common *P. falciparum* drug resistance loci: *crt*, *dhfr-ts*, *dhps*, *k13*, and *mdr1*.

## Requirements

[pixi](https://pixi.sh) — all other dependencies (Snakemake, snpEff, rsync) are managed by pixi automatically.

## Setup

1. Clone the repository
2. Edit `config/config.yaml` to point to your genome FASTA, GFF, and VCF files
3. Edit `config/targets.tsv` to define your mutations of interest (or leave the provided *P. falciparum* targets as-is)

## Running

```bash
pixi run snakemake
```

This installs all environments on first run and then executes the workflow with 8 cores.

## Configuration

`config/config.yaml`:

| Field | Description |
|---|---|
| `genome_fasta` | Path to the reference genome FASTA |
| `genome_gff` | Path to the genome annotation GFF |
| `merged_vcf` | Path to the genotyped cohort VCF (`.vcf` or `.vcf.gz`) |
| `genome_database_name` | Name for the snpEff database to be created |
| `database_description` | Human-readable description for the database |
| `targets_tsv` | Path to the targets file |
| `output_directory` | Directory where all results are written |

## Targets TSV

`config/targets.tsv` defines mutations of particular interest. Required columns:

| Column | Description |
|---|---|
| `CHROM` | Chromosome name matching the VCF |
| `POS` | 1-based position |
| `REF` / `ALT` | Reference and alternate alleles |
| `gene_name` | Short gene name (e.g. `k13`) |
| `aa_change_position` | Amino acid position label (e.g. `C580`) |
| `aminoacid_change` | Full AA change in 3-letter code (e.g. `Cys580Tyr`) |
| `gene_id` | Gene ID matching the GFF (e.g. `PF3D7_1343700`) |
| `mutation_name` | Display name for the mutation |

## Outputs

All outputs are written to `output_directory` (default: `results/`).

```
results/
├── config/               # Copy of config files used for this run
├── snpEff_processing/    # snpEff database and config
├── snpEff_output/        # Annotated VCF and snpEff summary HTML
└── AA_tables/
    ├── coverage_AA_table.csv    # Total read depth per sample per mutation
    ├── reference_AA_table.csv   # Reference allele depth
    ├── alternate_AA_table.csv   # Alternate allele depth
    └── other_AA_table.csv       # Depth from all other alleles at the site
```

Each CSV has mutations as columns and samples as rows, with 6 header rows (Gene ID, Gene, Mutation Name, ExonicFunc, AA Change, Targeted). Targeted mutations are flagged with `Yes` in the Targeted row; all other annotated variants are `No`.
