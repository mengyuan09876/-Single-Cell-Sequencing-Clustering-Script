import argparse
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import Align
from Bio.Align import substitution_matrices
import numpy as np
import os
import subprocess
import pickle
from sklearn.cluster import AgglomerativeClustering
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

def read_fasta_files(fasta_dir):
    sequences = []
    for filename in os.listdir(fasta_dir):
        if filename.endswith('.fasta'):
            filepath = os.path.join(fasta_dir, filename)
            for record in SeqIO.parse(filepath, 'fasta'):
                sequences.append(record)
    print(f"Total sequences read: {len(sequences)}")
    return sequences

def write_combined_fasta(sequences, output_file):
    SeqIO.write(sequences, output_file, 'fasta')
    print(f"Combined FASTA file created: {output_file}")

def align_sequences(input_file, output_file, threads):
    clustalomega_cline = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True, threads=threads)
    subprocess.call(str(clustalomega_cline), shell=True)
    print("Sequences aligned using Clustal Omega.")

def calculate_similarity_matrix(aligned_sequences):
    matrix = substitution_matrices.load("BLOSUM62")
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = matrix
    num_sequences = len(aligned_sequences)
    similarity_matrix = np.zeros((num_sequences, num_sequences))
    for i in range(num_sequences):
        for j in range(i, num_sequences):
            score = aligner.score(aligned_sequences[i].seq, aligned_sequences[j].seq)
            similarity_matrix[i, j] = score
            similarity_matrix[j, i] = score
    print("Similarity matrix calculated.")
    return similarity_matrix

def perform_clustering(similarity_matrix):
    clustering = AgglomerativeClustering(affinity='precomputed', linkage='average', distance_threshold=0, n_clusters=None)
    clustering.fit(similarity_matrix)
    return clustering.labels_

def save_similarity_matrix(similarity_matrix, output_file):
    np.save(output_file, similarity_matrix)
    print(f"Similarity matrix saved: {output_file}")

def save_cluster_labels(labels, output_file):
    with open(output_file, 'wb') as f:
        pickle.dump(labels, f)
    print(f"Cluster labels saved: {output_file}")

def plot_dendrogram(similarity_matrix, labels, output_file):
    linked = linkage(similarity_matrix, 'single')
    plt.figure(figsize=(10, 7))
    dendrogram(linked, labels=labels)
    plt.title('Dendrogram of Clustered Sequences')
    plt.xlabel('Sequences')
    plt.ylabel('Distance')
    plt.savefig(output_file, dpi=200)
    plt.show()
    print(f"Dendrogram saved: {output_file}")

def main(fasta_dir, results_dir, threads):
    os.makedirs(results_dir, exist_ok=True)

    # Step 1: Read FASTA files
    sequences = read_fasta_files(fasta_dir)

    # Step 2: Combine sequences into one FASTA file
    combined_fasta = os.path.join(results_dir, 'combined_sequences.fasta')
    write_combined_fasta(sequences, combined_fasta)

    # Step 3: Align sequences
    aligned_fasta = os.path.join(results_dir, 'aligned.fasta')
    align_sequences(combined_fasta, aligned_fasta, threads)

    # Step 4: Calculate similarity matrix
    aligned_sequences = list(SeqIO.parse(aligned_fasta, 'fasta'))
    similarity_matrix = calculate_similarity_matrix(aligned_sequences)

    # Step 5: Save similarity matrix
    similarity_matrix_file = os.path.join(results_dir, 'similarity_matrix.npy')
    save_similarity_matrix(similarity_matrix, similarity_matrix_file)

    # Step 6: Perform clustering
    labels = perform_clustering(similarity_matrix)

    # Step 7: Save cluster labels
    cluster_labels_file = os.path.join(results_dir, 'cluster_labels.pkl')
    save_cluster_labels(labels, cluster_labels_file)

    # Step 8: Plot and save dendrogram
    dendrogram_file = os.path.join(results_dir, 'dendrogram.png')
    plot_dendrogram(similarity_matrix, [seq.id for seq in aligned_sequences], dendrogram_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Cluster single-cell sequencing data based on similarity.')
    parser.add_argument('--fasta_dir', type=str, required=True, help='Directory containing FASTA files')
    parser.add_argument('--results_dir', type=str, required=True, help='Directory to save results')
    parser.add_argument('--threads', type=int, default=40, help='Number of threads to use for alignment (default: 40)')
    args = parser.parse_args()

    main(args.fasta_dir, args.results_dir, args.threads)
