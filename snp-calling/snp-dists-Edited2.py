#!/usr/bin/env python

import sys
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(description='Calculate pairwise SNP distances from a FASTA alignment')
parser.add_argument('infile', help='input FASTA file', type=str)
parser.add_argument('--outfile_simple', help='output TSV file for simple distances (default: stdout)', type=str, default=None)
parser.add_argument('--outfile_evolutionary', help='output TSV file for evolutionary distances', type=str, default=None)
args = parser.parse_args()

# Usage examples:
# Output both matrices to separate files
# python snp_distances.py alignment.fasta --outfile_simple simple_dist.tsv --outfile_evolutionary evo_dist.tsv

# Output simple to stdout, evolutionary to file
# python snp_distances.py alignment.fasta --outfile_evolutionary evo_dist.tsv

# Output just simple to stdout (default)
# python snp_distances.py alignment.fasta


# IUPAC ambiguity codes
IUPAC_CODES = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'},
    '-': set()
}

def read_fasta(filename):
    """Read FASTA file and return dict of sequences"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

# def count_evolutionary_distance(base1, base2):
#     """
#     Calculate minimum evolutionary distance between two bases/genotypes.
    
#     Returns the minimum number of mutations needed to go from base1 to base2.
#     Handles diploid and triploid genotypes.
    
#     Examples:
#     - A vs A = 0 (identical)
#     - A vs G = 1 (one substitution)
#     - A vs R (A/G) = 1 (A/A -> A/G, one mutation)
#     - R (A/G) vs Y (C/T) = 2 (A/G -> C/T, need to change both alleles)
#     - A vs V (A/C/G) = 2 (A/A/A -> A/C/G, two mutations)
#     - V (A/C/G) vs B (C/G/T) = 1 (A/C/G -> C/G/T, change A->T and keep C,G)
#     """
#     bases1 = IUPAC_CODES.get(base1, set())
#     bases2 = IUPAC_CODES.get(base2, set())
    
#     # If either is empty (gap) or N, return None (skip)
#     if not bases1 or not bases2 or base1 == 'N' or base2 == 'N':
#         return None
    
#     # If identical sets, distance = 0
#     if bases1 == bases2:
#         return 0
    
#     # Ploidy of each genotype
#     ploidy1 = len(bases1) # Number of distinct alleles in genotype 1
#     ploidy2 = len(bases2) # Number of distinct alleles in genotype 2
    
#     # Calculate overlap and unique bases
#     overlap = bases1 & bases2
#     unique_to_1 = bases1 - overlap
#     unique_to_2 = bases2 - overlap
    
#     # Minimum mutations needed
#     if overlap:
#         # Share some alleles - only need to change the non-overlapping ones
#         # Example: R (A/G) vs W (A/T) → overlap={A}, need G->T = 1 mutation
#         # Example: V (A/C/G) vs B (C/G/T) → overlap={C,G}, need A->T = 1 mutation
#         # The minimum is the size of the symmetric difference divided by 2, 
#         # but we need to account for ploidy changes
        
#         # If moving from lower to higher ploidy, need to add alleles
#         # If moving from higher to lower, need to remove alleles
#         # Simplification: count non-shared alleles that need changing
        
#         if ploidy1 == ploidy2:
#             # Same ploidy: need to change non-overlapping alleles
#             return len(unique_to_1)  # Same as len(unique_to_2)
#         else:
#             # Different ploidy: more complex
#             # Need max of how many to add/remove
#             return max(len(unique_to_1), len(unique_to_2))
    
#     else:
#         # No overlap - complete replacement needed
#         # Distance = max ploidy (worst case: change all alleles)
#         # Example: R (A/G) vs Y (C/T) → 2 mutations minimum
#         # Example: A (A/A/A) vs V (C/G/T) → 3 mutations
#         return max(ploidy1, ploidy2)

def count_evolutionary_distance(base1, base2):
    """
    Calculate minimum evolutionary distance between two bases/genotypes.
    Returns the minimum number of mutations needed to go from base1 to base2.
    """
    bases1 = IUPAC_CODES.get(base1, set())
    bases2 = IUPAC_CODES.get(base2, set())
    
    # If either is empty (gap) or N, return None (skip)
    if not bases1 or not bases2 or base1 == 'N' or base2 == 'N':
        return None
    
    # If identical sets, distance = 0
    if bases1 == bases2:
        return 0
    
    # Calculate overlap and unique bases
    overlap = bases1 & bases2
    unique_to_1 = bases1 - overlap
    unique_to_2 = bases2 - overlap
    
    if overlap:
        # Share some alleles - need to change the non-overlapping ones
        # Max handles both symmetric (R vs W) and asymmetric (R vs G) cases
        return max(len(unique_to_1), len(unique_to_2))
    else:
        # No overlap - complete replacement needed
        return max(len(bases1), len(bases2))

# def count_snp_differences(seq1, seq2, method='simple'):
#     """
#     Count SNP differences between two sequences.
    
#     method='simple': Any difference counts as 1
#     method='evolutionary': Minimum mutations needed
#     """
#     if len(seq1) != len(seq2):
#         sys.stderr.write("Warning: sequences have different lengths\n")
#         return None
    
#     differences = 0
#     compared_sites = 0
    
#     for i in range(len(seq1)):
#         base1 = seq1[i].upper()
#         base2 = seq2[i].upper()
        
#         # Skip if either is N or gap
#         if base1 in ['N', '-'] or base2 in ['N', '-']:
#             continue
        
#         compared_sites += 1
        
#         if method == 'simple':
#             # Count as difference if bases are not identical
#             if base1 != base2:
#                 differences += 1
#         elif method == 'evolutionary':
#             # Count minimum evolutionary distance
#             dist = count_evolutionary_distance(base1, base2)
#             if dist is not None:
#                 differences += dist
    
#     return differences
def count_snp_differences(seq1, seq2, method='simple'):
    """Count SNP differences between two sequences."""
    if len(seq1) != len(seq2):
        sys.stderr.write("Warning: sequences have different lengths\n")
        return None
    
    differences = 0
    compared_sites = 0
    skipped_simple = 0
    skipped_evo = 0
    
    for i in range(len(seq1)):
        base1 = seq1[i].upper()
        base2 = seq2[i].upper()
        
        # Skip if either is N or gap
        if base1 in ['N', '-'] or base2 in ['N', '-']:
            if method == 'simple':
                skipped_simple += 1
            continue
        
        compared_sites += 1
        
        if method == 'simple':
            if base1 != base2:
                differences += 1
        elif method == 'evolutionary':
            dist = count_evolutionary_distance(base1, base2)
            if dist is None:
                skipped_evo += 1
                compared_sites -= 1  # We thought we'd compare but didn't
                sys.stderr.write(f"Position {i}: {base1} vs {base2} returned None in evolutionary!\n")
            else:
                differences += dist
    
    #if method == 'simple':
    #    sys.stderr.write(f"Simple: compared {compared_sites}, skipped {skipped_simple}\n")
    #else:
    #    sys.stderr.write(f"Evolutionary: compared {compared_sites}, skipped {skipped_evo} extra positions\n")
    
    return differences


def calculate_distance_matrix(sequences, method='simple'):
    """Calculate pairwise SNP distances for all sequences"""
    isolate_names = list(sequences.keys())
    n = len(isolate_names)
    
    # Create distance matrix
    distances = defaultdict(dict)
    
    for i in range(n):
        for j in range(n):
            name1 = isolate_names[i]
            name2 = isolate_names[j]
            
            if i == j:
                distances[name1][name2] = 0
            elif name2 in distances[name1]:
                # Already calculated
                continue
            else:
                diff = count_snp_differences(sequences[name1], sequences[name2], method=method)
                distances[name1][name2] = diff
                distances[name2][name1] = diff
    
    return isolate_names, distances

def write_distance_matrix(isolate_names, distances, outfile=None, label=''):
    """Write distance matrix in TSV format"""
    output = sys.stdout if outfile is None else open(outfile, 'w')
    
    # Header row
    output.write('\t'.join([''] + isolate_names) + '\n')
    
    # Data rows
    for name1 in isolate_names:
        row = [name1]
        for name2 in isolate_names:
            row.append(str(distances[name1][name2]))
        output.write('\t'.join(row) + '\n')
    
    if outfile is not None:
        output.close()
        sys.stderr.write(f"{label} distance matrix written to {outfile}\n")
    else:
        sys.stderr.write(f"{label} distance matrix written to stdout\n")

# Main execution
sys.stderr.write(f"Reading sequences from {args.infile}\n")
sequences = read_fasta(args.infile)
sys.stderr.write(f"Found {len(sequences)} sequences\n")

# Verify all sequences have same length
seq_lengths = set(len(seq) for seq in sequences.values())
if len(seq_lengths) > 1:
    sys.stderr.write(f"ERROR: Sequences have different lengths: {seq_lengths}\n")
    sys.exit(1)

# Calculate simple distances
sys.stderr.write("Calculating simple pairwise distances...\n")
isolate_names, simple_distances = calculate_distance_matrix(sequences, method='simple')

# Calculate evolutionary distances
sys.stderr.write("Calculating evolutionary pairwise distances...\n")
_, evolutionary_distances = calculate_distance_matrix(sequences, method='evolutionary')

# Write outputs
if args.outfile_simple:
    write_distance_matrix(isolate_names, simple_distances, args.outfile_simple, label='Simple')
elif not args.outfile_evolutionary:
    # If no specific output, write simple to stdout
    write_distance_matrix(isolate_names, simple_distances, None, label='Simple')

if args.outfile_evolutionary:
    write_distance_matrix(isolate_names, evolutionary_distances, args.outfile_evolutionary, label='Evolutionary')

sys.stderr.write("Done!\n")

