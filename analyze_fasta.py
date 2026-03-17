#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
analyze_fasta.py - Pairwise k-mer-based sequence identity analysis for FASTA files.
"""

import re
import sys
import argparse
import numpy as np


class FastaAnalyzer:
    """Encapsulates FASTA parsing and pairwise k-mer sequence identity analysis."""

    def __init__(self, input_path: str, kmer_size: int, output_path: str | None) -> None:
        """
        Initialize the FastaAnalyzer and validate input arguments.

        Parameters
        ----------
        input_path : str
            Path to the input FASTA file.
        kmer_size : int
            Length of k-mers used for comparison (>= 1).
        output_path : str or None
            Path for output file; None means stdout.

        Raises
        ------
        SystemExit
            If any validation check fails.
        """
        import os
        if not os.path.isfile(input_path):
            sys.stderr.write(f"ERROR: Input file not found or not readable: {input_path}\n")
            sys.exit(1)
        if kmer_size < 1:
            sys.stderr.write(f"ERROR: k-mer size must be a positive integer >= 1, got {kmer_size}\n")
            sys.exit(1)
        self.input_path = input_path
        self.kmer_size = kmer_size
        self.output_path = output_path
        self.sequences: dict[str, str] = {}

    def parse_fasta(self) -> dict[str, str]:
        """
        Read and parse a FASTA file using the re module.

        Extracts sequence identifiers and concatenates multi-line sequences.
        Validates that sequence bodies contain only A, C, G, T, N (case-insensitive).
        Stores results in self.sequences.

        Returns
        -------
        dict[str, str]
            Mapping of sequence identifier -> nucleotide sequence (uppercase).

        Raises
        ------
        SystemExit
            If fewer than two sequences are found, or any sequence is shorter than k.
        """
        header_re = re.compile(r'^>(\S+)(.*)$')
        valid_re = re.compile(r'^[ACGTNacgtn]+$')

        sequences: dict[str, str] = {}
        current_id: str | None = None
        current_seq_parts: list[str] = []

        def _finalize(seq_id: str, parts: list[str]) -> None:
            seq = ''.join(parts).upper()
            if not valid_re.match(seq):
                invalid = re.sub(r'[ACGTNacgtn]', '', seq)
                sys.stderr.write(
                    f"WARNING: Unexpected characters in sequence '{seq_id}': "
                    f"{set(invalid)}\n"
                )
            sequences[seq_id] = seq

        with open(self.input_path, 'r', encoding='utf-8') as fh:
            for line in fh:
                line = line.rstrip('\n').strip()
                if not line:
                    continue
                m = header_re.match(line)
                if m:
                    if current_id is not None:
                        _finalize(current_id, current_seq_parts)
                    current_id = m.group(1)
                    current_seq_parts = []
                elif current_id is not None:
                    current_seq_parts.append(line)

        if current_id is not None:
            _finalize(current_id, current_seq_parts)

        if len(sequences) < 2:
            sys.stderr.write(
                "ERROR: FASTA file must contain at least two sequences for pairwise comparison.\n"
            )
            sys.exit(1)

        for seq_id, seq in sequences.items():
            if len(seq) < self.kmer_size:
                sys.stderr.write(
                    f"ERROR: Sequence '{seq_id}' (length {len(seq)}) is shorter than "
                    f"k-mer size {self.kmer_size}.\n"
                )
                sys.exit(1)

        self.sequences = sequences
        return sequences

    def extract_kmers(self, sequence: str) -> np.ndarray:
        """
        Extract all overlapping k-mers of length self.kmer_size from a sequence.

        Uses a sliding window of step 1, yielding L - k + 1 k-mers for a
        sequence of length L.

        Parameters
        ----------
        sequence : str
            The nucleotide sequence string (uppercase).

        Returns
        -------
        np.ndarray
            1-D array of k-mer strings.
        """
        k = self.kmer_size
        return np.array([sequence[i:i + k] for i in range(len(sequence) - k + 1)])

    def _multiset_intersection_count(self, kmers_a: np.ndarray, kmers_b: np.ndarray) -> int:
        """
        Compute the multiset intersection count between two k-mer arrays.

        For each unique k-mer, the contribution is min(count_in_a, count_in_b).

        Parameters
        ----------
        kmers_a : np.ndarray
            K-mers from the first sequence.
        kmers_b : np.ndarray
            K-mers from the second sequence.

        Returns
        -------
        int
            Total number of shared k-mers (multiset intersection size).
        """
        unique_a, counts_a = np.unique(kmers_a, return_counts=True)
        unique_b, counts_b = np.unique(kmers_b, return_counts=True)

        # Find common k-mers
        common_mask_a = np.isin(unique_a, unique_b)
        common_mask_b = np.isin(unique_b, unique_a)

        shared_a = counts_a[common_mask_a]
        shared_b = counts_b[common_mask_b]

        # Sort to align (np.isin preserves order; sort by the k-mer strings)
        order_a = np.argsort(unique_a[common_mask_a])
        order_b = np.argsort(unique_b[common_mask_b])

        return int(np.sum(np.minimum(shared_a[order_a], shared_b[order_b])))

    def compute_identity_matrix(self) -> np.ndarray:
        """
        Compute the pairwise sequence identity matrix using k-mer overlap.

        Identity for a pair (i, j) is defined as the multiset intersection
        count of their k-mers divided by the number of k-mers in the shorter
        sequence.

        Returns
        -------
        np.ndarray
            Symmetric (n x n) float64 matrix of identity values; diagonal is 1.0.
        """
        ids = list(self.sequences.keys())
        n = len(ids)
        matrix = np.zeros((n, n), dtype=np.float64)

        all_kmers = [self.extract_kmers(self.sequences[seq_id]) for seq_id in ids]

        for i in range(n):
            matrix[i, i] = 1.0
            for j in range(i + 1, n):
                shared = self._multiset_intersection_count(all_kmers[i], all_kmers[j])
                shorter = min(len(all_kmers[i]), len(all_kmers[j]))
                identity = shared / shorter if shorter > 0 else 0.0
                matrix[i, j] = identity
                matrix[j, i] = identity

        return matrix

    def write_output(self, matrix: np.ndarray) -> None:
        """
        Format and write the identity matrix as a tab-separated table.

        Sequence identifiers are used as both row and column headers.
        Values are formatted to four decimal places.
        Uses sys.stdout.write() for stdout output.

        Parameters
        ----------
        matrix : np.ndarray
            The (n x n) pairwise identity matrix.

        Returns
        -------
        None
        """
        ids = list(self.sequences.keys())
        header = 'ID\t' + '\t'.join(ids) + '\n'
        rows = []
        for i, seq_id in enumerate(ids):
            values = '\t'.join(f'{matrix[i, j]:.4f}' for j in range(len(ids)))
            rows.append(f'{seq_id}\t{values}\n')

        output = header + ''.join(rows)

        if self.output_path:
            with open(self.output_path, 'w', encoding='utf-8') as fh:
                fh.write(output)
        else:
            sys.stdout.write(output)


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed argument object with attributes: input, kmer, output.
    """
    parser = argparse.ArgumentParser(
        description='Pairwise k-mer-based sequence identity analysis for FASTA files.'
    )
    parser.add_argument(
        '--input', '-i',
        type=str,
        required=True,
        help='Path to the input FASTA file.'
    )
    parser.add_argument(
        '--kmer', '-k',
        type=int,
        default=1,
        help='K-mer size for comparison (default: 1).'
    )
    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Path to output file (default: stdout).'
    )
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    analyzer = FastaAnalyzer(
        input_path=args.input,
        kmer_size=args.kmer,
        output_path=args.output
    )
    analyzer.parse_fasta()
    matrix = analyzer.compute_identity_matrix()
    analyzer.write_output(matrix)
