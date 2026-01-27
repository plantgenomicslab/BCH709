#!/usr/bin/env python3
"""
Trypsin digestion analysis of protein sequences.
Trypsin cleaves after K (lysine) and R (arginine), except when followed by P (proline).
"""

import re
from pathlib import Path


def read_fasta(filepath):
    """Read sequences from a FASTA file."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)

        # Add the last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def trypsin_digest(sequence):
    """
    Digest a protein sequence with trypsin.
    Trypsin cleaves after K and R, unless followed by P.

    Returns a list of tuples: (fragment_sequence, cleavage_site_position)
    """
    fragments = []
    current_fragment = []

    for i, aa in enumerate(sequence):
        current_fragment.append(aa)

        # Check if this is a cleavage site
        # Cleave after K or R, but not if followed by P
        if aa in ['K', 'R']:
            # Check if this is not the last residue and next is not P
            if i < len(sequence) - 1:
                if sequence[i + 1] != 'P':
                    # This is a cleavage site
                    fragment_seq = ''.join(current_fragment)
                    fragments.append((fragment_seq, i + 1))  # Position is 1-indexed
                    current_fragment = []
            else:
                # Last residue, still cleave
                fragment_seq = ''.join(current_fragment)
                fragments.append((fragment_seq, i + 1))
                current_fragment = []

    # Add any remaining sequence as the last fragment
    if current_fragment:
        fragment_seq = ''.join(current_fragment)
        fragments.append((fragment_seq, len(sequence)))

    return fragments


def main():
    # Read the FASTA file
    fasta_file = Path(__file__).parent / 'seqs.fa'
    sequences = read_fasta(fasta_file)

    # Output file
    output_file = Path(__file__).parent / 'trypsin_fragments.txt'

    with open(output_file, 'w') as out:
        # Write header
        out.write('Sequence_Name\tFragment_Sequence\tFragment_Length\tCleavage_Site_Position\n')

        # Process each sequence
        for seq_name, seq in sequences.items():
            fragments = trypsin_digest(seq)

            for fragment_seq, cleavage_pos in fragments:
                out.write(f'{seq_name}\t{fragment_seq}\t{len(fragment_seq)}\t{cleavage_pos}\n')

    print(f'Trypsin digestion complete!')
    print(f'Results written to: {output_file}')
    print(f'\nSummary:')
    print(f'Total sequences processed: {len(sequences)}')

    total_fragments = 0
    for seq_name, seq in sequences.items():
        fragments = trypsin_digest(seq)
        print(f'\n{seq_name}')
        print(f'  Original length: {len(seq)} aa')
        print(f'  Number of fragments: {len(fragments)}')
        total_fragments += len(fragments)

    print(f'\nTotal fragments generated: {total_fragments}')


if __name__ == '__main__':
    main()
