# Helixer vs Expert Annotation Comparison Pipeline

**Evaluating Helixer Deep-Learning Gene Predictions Against Reference Annotations**
*Internship — GBOT, Université Paris-Saclay*

-------
This project implements a comparison pipeline between gene annotations predicted by **Helixer** (a deep-learning gene prediction tool) and expert reference annotations (TAIR12, ARAPORT11, B73, and other plant genomes), built on top of the **GBOT** framework.

The pipeline queries a FLAGdb instance for expert and Helixer-predicted features (CDS, mRNA, transposable elements, PFAM domains), classifies them into **Added / Missed / Common** categories per species, and produces statistics, BED/GFF3 files, protein-level functional annotation (via BLASTp against UniProt/SwissProt), UTR difference reports, transposable element overlap analysis, promoter motif (PLM) extraction, and BUSCO completeness scores — along with the R-generated figures used in the final report.

-----------------------------------------------------
## Repository structure

```
Projet/
├── Assemblies/               # Reference genome assemblies + expert TE annotations (TES/)
├── Rapport/                  # Report generation
│   ├── CodesFig/             # R scripts used to produce the figures
│   ├── Figures/              # Figures included in the report
│   ├── commands.txt          # Commands used to run each species analysis
│   └── resume.txt            # Aggregated per-species summary
├── Species/                  # Per-species output directories (created by the script)
│   └── <species>/            # seqs, figs, Stats, Bed, TES, UTRAnalysis, Analysis, PLM, BlastAnalysisADDED
├── Uniprot/                  # Local SwissProt BLAST database (downloaded on --makeDB)
├── Assemblies.zip            # Archived copy of the Assemblies/ folder
└── find_differences.py       # Main analysis script (entry point)
```

-------------------------------------------------------------------------------------

## How it works

For a given species, `find_differences.py`:

1. Pulls expert and `helixer_*` features from FLAGdb (CDS, mRNA, TEs, PFAM).
2. Overlaps them to classify each feature as **Added** (Helixer-only), **Missed** (expert-only), or **Common**.
3. Runs whichever analyses are requested (see flags below), writing BED/TSV outputs, FASTA files, and calling the R scripts in `Rapport/CodesFig/` to produce figures.
4. Appends a one-line summary per species to `Rapport/resume.txt`.

## How to use

1. **Install the required dependencies**

Python packages (FLAGdb server modules, plus the usual bioinformatics stack), and the following external tools on `$PATH`:
- `blastp` / `blastn` / `makeblastdb` (BLAST+)
- `bedtools`, `samtools`
- `R` / `Rscript` (for the figures in `Rapport/CodesFig/`)
- `busco`
- A Chrome/Chromedriver setup for Selenium (only needed for `--PLMView`)

2. **Run an analysis for a species**

```bash
python find_differences.py --dir Species/TAIR12 --gender <schema> --speId 34 --doAll True
```

Or run only specific analyses instead of `--doAll`:

```bash
python find_differences.py --dir Species/B73 --gender <schema> --speId 23 \
    --CatAnalysis True --TEAnalysis True --BUSCO True
```

### Command-line options

| Flag | Description |
|---|---|
| `--dir` | Destination output directory (required) |
| `--gender` | FLAGdb schema to query (required) |
| `--speId` | Source species ID in FLAGdb (required) |
| `--doAll` | Run the full analysis (all flags below) |
| `--UTRAnalysis` | Compare 5'/3' UTR regions between Helixer and expert mRNAs |
| `--RNAAnalysis` | Compare mRNA length and exon counts, Helixer vs reference |
| `--CatAnalysis` | Added/Missed/Common functional analysis (BLASTp vs UniProt) |
| `--TEAnalysis` | Analyze transposable elements predicted by Helixer |
| `--PLMView` | Submit Helixer promoter sequences to PLMView for motif detection |
| `--Email` | Email address required by PLMView |
| `--makeDB` | (Re)download and build the local UniProt/SwissProt BLAST database |
| `--proteinCompare` | Compare Helixer vs reference proteins |
| `--BUSCO` | Run BUSCO on the Added/Missed/Common protein sets |
| `--generate` | Generate the per-category protein `.fasta` files |

### Reference species IDs

| ID | Species |
|---|---|
| 34 | *Arabidopsis thaliana* 
| 20 | *Arabidopsis thaliana* |
| 23 | *Zea mays* | 
| 30 | *Arabidopsis lyrata* | 
| 22 | *Brachypodium distachyon* | 
| 26 | *Malus domestica* | 
| 28 | *Rosa chinensis* | 
| 24 | *Solanum lycopersicum* | 
