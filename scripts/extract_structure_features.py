#!/usr/bin/env python3
"""
Extract RNA secondary structure features around informative read pairs.

Uses ViennaRNA command-line tools (RNAfold) instead of Python bindings.
"""

import pysam
import argparse
import sys
from collections import defaultdict
import pandas as pd
import numpy as np
import subprocess
import tempfile
import os

def check_rnafold():
    """Check if RNAfold is available."""
    try:
        result = subprocess.run(['RNAfold', '--version'], 
                              capture_output=True, text=True)
        return True
    except FileNotFoundError:
        return False

def extract_sequence(fasta_file, chrom, start, end, strand='+'):
    """Extract genomic sequence from fasta file."""
    with pysam.FastaFile(fasta_file) as fasta:
        seq = fasta.fetch(chrom, start, end)
        
        # Convert to RNA and handle strand
        seq = seq.upper().replace('T', 'U')
        
        if strand == '-':
            # Reverse complement
            complement = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
            seq = ''.join(complement.get(base, 'N') for base in reversed(seq))
    
    return seq

def calculate_mfe_structure_cmdline(sequence):
    """
    Calculate minimum free energy structure using RNAfold command line.
    
    Returns:
        dict with MFE, structure, ensemble_energy
    """
    # Create temporary file for sequence
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f">seq\n{sequence}\n")
        temp_file = f.name
    
    try:
        # Run RNAfold with partition function
        result = subprocess.run(
            ['RNAfold', '--noPS', '-p'],
            stdin=open(temp_file),
            capture_output=True,
            text=True
        )
        
        # Parse output
        lines = result.stdout.strip().split('\n')
        
        # Line 2: structure and MFE
        # Format: "(((...)))  (-12.34)"
        structure_line = lines[1]
        structure = structure_line.split()[0]
        mfe = float(structure_line.split('(')[1].rstrip(')'))
        
        # Line 3: ensemble energy (from partition function)
        # Format: "(((...))) [-12.45]"
        if len(lines) > 2:
            ensemble_line = lines[2]
            ensemble_energy = float(ensemble_line.split('[')[1].rstrip(']'))
        else:
            ensemble_energy = None
        
        return {
            'mfe': mfe,
            'structure': structure,
            'ensemble_energy': ensemble_energy,
            'sequence_length': len(sequence)
        }
    
    except Exception as e:
        print(f"Warning: RNAfold failed: {e}")
        return {
            'mfe': None,
            'structure': None,
            'ensemble_energy': None,
            'sequence_length': len(sequence)
        }
    finally:
        os.unlink(temp_file)

def calculate_accessibility(sequence, structure):
    """
    Calculate unpaired probabilities (accessibility) from structure.
    Simply count dots (unpaired) vs parentheses (paired).
    
    For more detailed accessibility, would need RNAplfold command line tool.
    """
    unpaired_count = structure.count('.')
    total_bases = len(structure)
    
    if total_bases == 0:
        return 0.0
    
    return unpaired_count / total_bases

def calculate_local_structure_windows(sequence, window_sizes=[50, 100, 200]):
    """
    Calculate MFE in sliding windows of different sizes.
    This captures local structural propensity.
    """
    results = {}
    
    for window_size in window_sizes:
        if len(sequence) < window_size:
            results[f'mfe_window_{window_size}_mean'] = None
            results[f'mfe_window_{window_size}_min'] = None
            results[f'mfe_window_{window_size}_std'] = None
            continue
        
        mfe_values = []
        for i in range(0, len(sequence) - window_size + 1, 10):  # Step by 10 to speed up
            window_seq = sequence[i:i+window_size]
            window_result = calculate_mfe_structure_cmdline(window_seq)
            if window_result['mfe'] is not None:
                mfe_values.append(window_result['mfe'])
        
        if mfe_values:
            # Summary statistics for this window size
            results[f'mfe_window_{window_size}_mean'] = np.mean(mfe_values)
            results[f'mfe_window_{window_size}_min'] = np.min(mfe_values)
            results[f'mfe_window_{window_size}_std'] = np.std(mfe_values)
        else:
            results[f'mfe_window_{window_size}_mean'] = None
            results[f'mfe_window_{window_size}_min'] = None
            results[f'mfe_window_{window_size}_std'] = None
    
    return results

def calculate_bond_energies(sequence, structure):
    """
    Calculate energy contribution of individual base pairs.
    
    Note: For detailed bond energies with command-line tools,
    you would need to parse RNAeval output. For now, we just
    count base pairs and estimate.
    """
    # Get base pairs from structure
    pairs = []
    stack = []
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs.append((j, i))
    
    # Simple statistics
    return {
        'num_base_pairs': len(pairs),
        'pairing_density': len(pairs) / len(sequence) if len(sequence) > 0 else 0
    }

def analyze_read_pair_structure(read1, read2, fasta_file, 
                               flank_size=100, 
                               window_sizes=[50, 100, 200]):
    """
    Analyze RNA structure around an informative read pair.
    
    Args:
        read1, read2: pysam AlignedSegment objects
        fasta_file: path to genome fasta
        flank_size: bp to extract on each side of read
        window_sizes: sizes for sliding window analysis
    
    Returns:
        dict of structure features
    """
    features = {
        'read_name': read1.query_name,
        'chrom': read1.reference_name,
    }
    
    # Extract sequence around each read
    for read_idx, read in enumerate([read1, read2], 1):
        start = max(0, read.reference_start - flank_size)
        end = read.reference_end + flank_size
        
        # Get genomic sequence
        try:
            sequence = extract_sequence(fasta_file, read.reference_name, 
                                       start, end, strand='+')
        except Exception as e:
            print(f"Warning: Could not extract sequence for {read.query_name}: {e}")
            return None
        
        # Calculate MFE for full window
        mfe_results = calculate_mfe_structure_cmdline(sequence)
        features[f'read{read_idx}_mfe'] = mfe_results['mfe']
        features[f'read{read_idx}_ensemble_energy'] = mfe_results['ensemble_energy']
        features[f'read{read_idx}_structure'] = mfe_results['structure']
        
        # Calculate accessibility (simplified version)
        if mfe_results['structure']:
            accessibility = calculate_accessibility(sequence, mfe_results['structure'])
            features[f'read{read_idx}_accessibility'] = accessibility
        else:
            features[f'read{read_idx}_accessibility'] = None
        
        # Calculate local structure in windows
        window_results = calculate_local_structure_windows(sequence, window_sizes)
        for key, value in window_results.items():
            features[f'read{read_idx}_{key}'] = value
        
        # Store sequence info
        features[f'read{read_idx}_seq_length'] = len(sequence)
        features[f'read{read_idx}_gc_content'] = (sequence.count('G') + sequence.count('C')) / len(sequence)
    
    return features

def main():
    parser = argparse.ArgumentParser(
        description='Extract RNA secondary structure features from informative read pairs'
    )
    parser.add_argument('--bam', required=True,
                       help='BAM file with informative read pairs')
    parser.add_argument('--fasta', required=True,
                       help='Genome FASTA file (indexed)')
    parser.add_argument('--output', required=True,
                       help='Output TSV file with structure features')
    parser.add_argument('--flank-size', type=int, default=100,
                       help='Flanking sequence size (bp) around each read (default: 100)')
    parser.add_argument('--window-sizes', type=int, nargs='+', 
                       default=[50, 100, 200],
                       help='Window sizes for local structure analysis')
    parser.add_argument('--max-reads', type=int, default=None,
                       help='Maximum number of read pairs to process (for testing)')
    
    args = parser.parse_args()
    
    # Check if RNAfold is available
    if not check_rnafold():
        print("ERROR: RNAfold not found in PATH")
        print("Make sure ViennaRNA is installed: conda install -c bioconda viennarna")
        sys.exit(1)
    
    print("Using RNAfold command-line tool")
    print(f"Analyzing RNA structure for read pairs in {args.bam}")
    print(f"Flanking size: {args.flank_size} bp")
    print(f"Window sizes: {args.window_sizes}")
    
    # Load BAM file
    bam = pysam.AlignmentFile(args.bam, 'rb')
    
    # Group reads by name to get pairs
    reads_by_name = defaultdict(list)
    total_reads = 0
    
    print("Loading read pairs...")
    for read in bam:
        if not read.is_unmapped:
            reads_by_name[read.query_name].append(read)
            total_reads += 1
            
            if args.max_reads and total_reads >= args.max_reads * 2:
                break
    
    print(f"Found {len(reads_by_name)} unique read pairs")
    
    # Analyze each pair
    all_features = []
    analyzed = 0
    failed = 0
    
    for query_name, reads in reads_by_name.items():
        if len(reads) != 2:
            continue
        
        analyzed += 1
        if analyzed % 1000 == 0:
            print(f"  Processed {analyzed} read pairs...")
        
        read1, read2 = reads
        features = analyze_read_pair_structure(
            read1, read2, args.fasta, 
            args.flank_size, args.window_sizes
        )
        
        if features:
            all_features.append(features)
        else:
            failed += 1
        
        if args.max_reads and analyzed >= args.max_reads:
            break
    
    print(f"\nAnalyzed {analyzed} read pairs")
    print(f"Successfully processed: {len(all_features)}")
    print(f"Failed: {failed}")
    
    # Convert to DataFrame and save
    df = pd.DataFrame(all_features)
    df.to_csv(args.output, sep='\t', index=False)
    
    print(f"\nStructure features saved to {args.output}")
    print(f"Columns: {list(df.columns)}")

if __name__ == "__main__":
    main()