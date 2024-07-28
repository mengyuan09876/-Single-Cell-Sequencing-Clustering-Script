# Single-Cell Sequencing Clustering Script

## Overview

This script is designed to process single-cell sequencing data by reading multiple FASTA files, combining them, aligning the sequences, calculating a similarity matrix, performing clustering, and saving the results. The script uses Clustal Omega for sequence alignment and scikit-learn for clustering.

## Requirements

- Python 3.x
- Biopython
- Clustal Omega
- scikit-learn
- Matplotlib
- NumPy

## Installation

1. **Install Biopython**:
   ```bash
   pip install biopython
   ```

2. **Install scikit-learn**:
   ```bash
   pip install scikit-learn
   ```

3. **Install Matplotlib**:
   ```bash
   pip install matplotlib
   ```

4. **Install Clustal Omega**:
   - On Ubuntu:
     ```bash
     sudo apt-get install clustalo
     ```
   - On macOS:
     ```bash
     brew install clustal-omega
     ```
   - Or download from the [Clustal Omega website](http://www.clustal.org/omega/).

## Usage

1. **Run the script with command-line arguments**:
   ```bash
   python cluster_sequences_with_saving.py --fasta_dir /path/to/fasta/files --results_dir /path/to/results --threads 40
   ```

   - `--fasta_dir`: Directory containing the FASTA files.
   - `--results_dir`: Directory to save the results.
   - `--threads`: Number of threads to use for alignment (default: 40).

## Script Details

### Step-by-Step Process

1. **Read FASTA Files**:
   - The script reads all FASTA files from the specified directory and combines the sequences into a single list.

2. **Combine Sequences**:
   - All sequences are written to a combined FASTA file for alignment.

3. **Align Sequences**:
   - Clustal Omega is used to align the combined sequences using the specified number of threads. The aligned sequences are saved to a file.

4. **Calculate Similarity Matrix**:
   - A similarity matrix is generated based on the alignment scores of the sequences.

5. **Perform Clustering**:
   - Hierarchical clustering is performed using the similarity matrix. Cluster labels are assigned to each sequence.

6. **Save Results**:
   - The combined FASTA file, aligned sequences file, similarity matrix, and cluster labels are saved to the results directory.
   - A dendrogram plot is generated and saved, visualizing the hierarchical clustering of the sequences.

### Output Files

- `combined_sequences.fasta`: Combined FASTA file of all sequences.
- `aligned.fasta`: Aligned sequences file.
- `similarity_matrix.npy`: Similarity matrix file (NumPy array).
- `cluster_labels.pkl`: Cluster labels file (Pickle format).
- `dendrogram.png`: Dendrogram plot image (PNG format).

### Example Output

After running the script, the results directory will contain the following files:

```
single_cell_results/
├── combined_sequences.fasta
├── aligned.fasta
├── similarity_matrix.npy
├── cluster_labels.pkl
└── dendrogram.png
```

These files include all necessary data to understand the clustering of your single-cell sequencing results.

## License

This script is provided "as is", without warranty of any kind. Use it at your own risk.
