# FASTA Sequence Identity Analysis

A Python 3 command-line tool for computing pairwise k-mer-based sequence identity from FASTA files. Given a set of nucleotide sequences, the script extracts overlapping k-mers from each sequence and computes a symmetric identity matrix using multiset intersection.

This script was developed for the analysis of partial coding sequences of the NaV1.4 alpha subunit (*SCN4A* gene) from four *Phyllobates* poison-dart frog species, retrieved from NCBI GenBank.

---

## Installation

### 1. Clone or download the repository

```bash
git clone https://github.com/TahmidA139/fasta_analysis.git
cd fasta_analysis
```

### 2. Create and activate the Conda environment

```bash
conda env create -f environment.yml
conda activate fasta_analysis
```

This installs Python 3.11 and NumPy from the `conda-forge` channel.

---

## Usage

```
python analyze_fasta.py --input <fasta_file> [--kmer <k>] [--output <output_file>]
```

### Arguments

| Argument         | Short | Required | Default | Description                                      |
|------------------|-------|----------|---------|--------------------------------------------------|
| `--input`        | `-i`  | Yes      | —       | Path to the input FASTA file                     |
| `--kmer`         | `-k`  | No       | `1`     | K-mer size for comparison (integer >= 1)         |
| `--output`       | `-o`  | No       | stdout  | Path to output file; omit to print to stdout     |

### Example: k = 1

```bash
python analyze_fasta.py --input sequences.fasta --kmer 1
```

```
ID          MK278857.1  MK278840.1  MK278831.1  MK278830.1
MK278857.1  1.0000      0.9982      1.0000      1.0000
MK278840.1  0.9982      1.0000      0.9963      0.9963
MK278831.1  1.0000      0.9963      1.0000      1.0000
MK278830.1  1.0000      0.9963      1.0000      1.0000
```

### Example: k = 2

```bash
python analyze_fasta.py --input sequences.fasta --kmer 2
```

```
ID          MK278857.1  MK278840.1  MK278831.1  MK278830.1
MK278857.1  1.0000      0.9945      1.0000      1.0000
MK278840.1  0.9945      1.0000      0.9927      0.9927
MK278831.1  1.0000      0.9927      1.0000      1.0000
MK278830.1  1.0000      0.9927      1.0000      1.0000
```

### Example: k = 3

```bash
python analyze_fasta.py --input sequences.fasta --kmer 3
```

```
ID          MK278857.1  MK278840.1  MK278831.1  MK278830.1
MK278857.1  1.0000      0.9872      1.0000      1.0000
MK278840.1  0.9872      1.0000      0.9853      0.9853
MK278831.1  1.0000      0.9853      1.0000      1.0000
MK278830.1  1.0000      0.9853      1.0000      1.0000
```

### Save output to a file

```bash
python analyze_fasta.py --input sequences.fasta --kmer 3 --output results.tsv
```

---

## Output Format

The output is a **tab-separated table** printed to stdout (or written to a file if `--output` is provided):

- The first row is a header containing the string `ID` followed by all sequence identifiers.
- Each subsequent row begins with a sequence identifier followed by its identity values against all sequences (including itself).
- Values are formatted to **four decimal places**.
- The diagonal is always `1.0000` (a sequence is identical to itself).
- The matrix is **symmetric**: identity(i, j) == identity(j, i).

### Identity calculation

For a pair of sequences *A* and *B*, identity is computed as:

```
identity(A, B) = |kmers(A) ∩ kmers(B)|  /  min(|kmers(A)|, |kmers(B)|)
```

where `∩` denotes **multiset intersection** (repeated k-mers count up to the number of times they appear in the shorter sequence), and `|kmers(X)|` is the total k-mer count for sequence *X* (= sequence length − k + 1).

---

## Project Structure

```
fasta_analysis/
├── analyze_fasta.py      # Main Python script
├── sequences.fasta       # Input FASTA file (4 Phyllobates SCN4A sequences)
├── README.md             # This file
├── LICENSE               # MIT License
└── environment.yml       # Conda environment specification
```

---

## References

- Pearson, W. R., & Lipman, D. J. (1988). Improved tools for biological sequence comparison. *Proceedings of the National Academy of Sciences*, 85(8), 2444–2448. https://doi.org/10.1073/pnas.85.8.2444
- Compeau, P. E. C., Pevzner, P. A., & Tesler, G. (2011). How to apply de Bruijn graphs to genome assembly. *Nature Biotechnology*, 29(11), 987–991. https://doi.org/10.1038/nbt.2023
- Tarvin, R. D., et al. (2017). Interacting amino acid replacements allow poison frogs to evolve epibatidine resistance. *Science*, 357(6357), 1261–1266. https://doi.org/10.1126/science.aan5061
- NCBI GenBank accessions: MK278857.1, MK278840.1, MK278831.1, MK278830.1
