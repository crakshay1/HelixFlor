# Helixflor
![Python](https://img.shields.io/badge/-Python-3776AB?logo=python&logoColor=white&style=flat-square) 
![R](https://img.shields.io/badge/-R-276DC3?logo=r&logoColor=white&style=flat-square) 
![SQL](https://img.shields.io/badge/-SQL-003B57?logo=sqlite&logoColor=white&style=flat-square) 

## 🔍 Context

Although high-throughput sequencing provides access to complete genomes, the structural annotation of genes in these genomes still remains a key step, especially in plants that have complex genomes (polyploidy, numerous transposable elements). The recent application of deep learning in annotation tools will surely make it possible to go faster in proposing annotations on both the structural and functional sides.

## Objective of the Internship

The GBOT database contains 6 plant genomes, all of which contain the official annotation of these genomes plus for some the annotation generated after the use of Helixer (Stiehler et al. 2021), an annotation tool that combines deep neural networks and HMM-type models to predict gene models from the genomic sequence alone. The internship consists in :

1. Applying and understanding Helixer on genomes not yet annotated

2. Redoing a global comparison on each genome

3. Targeting new genes defined by Helixer and highlight their characteristics on the structural side (gene size, number of exons) and functional side (generated protein and functional annotation)

4. Targeting genes corresponding to known genes but without 5’UTR, and take stock of the properties of these genes, check if TATA-box near the new annotated UTR



⚠️ Special thanks to **Franck SAMSON** 

---

## 🔍 Features

* **Interactive GUI**: Easy-to-use selection for species identifiers and input files via `easygui`.
* **Transposable Element (TE) Analysis**:
    * Local alignment against `TrepDB` using local BLASTn.
    * Local alignment against `URGIDB` using local BLASTn.
    * Automated web-submission to **Censor** (Giri Institute) for soft-masking and TE detection.
* **Gene & Protein Prediction**:
    * **Augustus**: Local prediction executed on both raw and masked sequences.
    * **FGENESH**: Automated web-based prediction (optimized for plant genomes) executed on both raw and masked sequences.
* **Automated Validation**:
    * **BLASTn**: Validation of mRNA against the TSA database.
    * **BLASTp**: Protein validation against NR and SwissProt databases.
    * **BLASTx**: Final check of the full sequence against SwissProt to catch missed coding regions.
* **Prediction Visualization**:  .gff3 files are generated for Artemis visualization.
* **Anonymity & Access**: Integrates with the **Tor** network to manage web requests.

---

## 📖 Prerequisites

### System Requirements
* **Linux/Unix** environment (recommended), if you don't have Linux, WSL can be used.
* **Mozilla Firefox** and **Geckodriver** (for Selenium automation).
* **Tor Service** (installed and configurable via `sudo service tor`).

### External Tools
The following must be installed and available in your system `$PATH`:
* **EMBOSS** (makeblastdb, blastn, blastp, blastx).
* **Augustus**.
* **Artemis**.

### Python Dependencies
```bash
pip install -r requirements.txt
```

---

## 📁 Setup & Directory Structure

To function correctly, the script requires a specific environment setup and will generate organized results as follows:

### 1. Pre-run Requirements
Every folder in this git has to be downloaded, except for Results.

### 2. Output Organization
The script automatically creates a `Results/` directory. Within that folder, a sub-directory is created for each run using the format `[Species]_[Filename]`.

### 3. Folder Contents
Inside each result folder, you will find:
* **`ARTEMIS/`**: Contains the images obtained via Artemis, and the .gff3 files of the predictions.
* **`TE/`**: Contains the local BLAST outputs against TrepDB and URGIDB + the `result.pdf` masking report from Censor + Dotplots and .fasta files of the TEs found by Censor.
* **`pred_by_augustus/`**: Includes the full Augustus text output, as well as sub-folders for predicted proteins and mRNA in `.fasta` format.
* **`pred_by_fgenesh/`**: (If selected) Contains the parsed PDF results and corresponding `.fasta` predictions.
* **`blastp_results/`**: Organized by database (`nr/` and `swiss/`), containing the alignment reports for predicted proteins.
* **`blastn_results/`**: Contains the validation reports for predicted mRNA against the TSA database.
* **`blastx_result/`**: Contains a final validation check of the original sequence against protein databases.
