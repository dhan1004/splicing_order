#!/usr/bin/env python3
"""
Extract RNA secondary structure features for individual introns in pairs.
Folds each intron separately to assess splice site accessibility.
"""

import pysam
import argparse
import sys
import pandas as pd
import subprocess
import tempfile
import os

def extract_sequence(fasta_file, chrom, start, end):
    """Extract genomic sequence using pysam."""
    try:
        with pysam.FastaFile(fasta_file) as fasta:
            sequence = fasta.fetch(chrom, start, end)
            sequence = sequence.upper().replace('T', 'U')
            return sequence
    except Exception as e:
        print(f"Warning: Could not extract sequence {chrom}:{start}-{end}: {e}", 
              file=sys.stderr)
        return None

def fold_rna_with_rnafold(sequence):
    """Fold RNA sequence using RNAfold."""
    if not sequence or len(sequence) < 10:
        return None
    
    try:
        # Write directly to stdin
        input_str = f">seq\n{sequence}\n"
        
        result = subprocess.run(
            ['RNAfold', '--noPS', '-p'],
            input=input_str,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode != 0:
            return None
        
        lines = result.stdout.strip().split('\n')
        
        # RNAfold output with -p flag:
        # Line 0: >header
        # Line 1: sequence
        # Line 2: MFE structure
        # Line 3: centroid structure (optional)
        # Line 4+: ensemble info
        
        if len(lines) < 3:
            return None
        
        # MFE structure is on line 2
        structure_line = lines[2].strip()
        
        # Check if this looks like a structure
        if '.' not in structure_line or '(' not in structure_line:
            return None
        
        parts = structure_line.split()
        
        if len(parts) < 1:
            return None
        
        structure = parts[0]
        
        # Extract MFE (in parentheses)
        mfe = None
        if len(parts) >= 2:
            for part in parts[1:]:
                if '(' in part and ')' in part:
                    try:
                        mfe_str = part.strip('()')
                        mfe = float(mfe_str)
                        break
                    except ValueError:
                        continue
        
        # Extract ensemble energy from later lines
        ensemble_energy = None
        for line in lines[3:]:
            if '[' in line and ']' in line:
                try:
                    start_idx = line.rfind('[')
                    end_idx = line.rfind(']')
                    if start_idx != -1 and end_idx != -1:
                        ensemble_str = line[start_idx+1:end_idx].strip()
                        ensemble_energy = float(ensemble_str)
                        break
                except (ValueError, IndexError):
                    continue
        
        if len(structure) != len(sequence):
            return None
        
        return {
            'structure': structure,
            'mfe': mfe,
            'ensemble_energy': ensemble_energy,
            'length': len(structure)
        }
    
    except subprocess.TimeoutExpired:
        return None
    except Exception as e:
        return None

def calculate_accessibility(structure, position, window=10):
    """
    Calculate accessibility around a position.
    Returns fraction of unpaired bases in window.
    """
    if position < 0 or position >= len(structure):
        return None
    
    start = max(0, position - window)
    end = min(len(structure), position + window + 1)
    
    window_structure = structure[start:end]
    unpaired = window_structure.count('.')
    total = len(window_structure)
    
    return unpaired / total if total > 0 else 0.0

def analyze_intron_structure(sequence, intron_length, padding=0):
    """
    Fold an intron sequence and calculate splice site accessibility.
    
    Args:
        sequence: RNA sequence including intron and padding
        intron_length: Length of the actual intron (excluding padding)
        padding: Amount of padding on each side
    
    Returns:
        dict with structure features
    """
    # Fold the sequence
    fold_result = fold_rna_with_rnafold(sequence)

    print("fold result: ", fold_result)
    
    if not fold_result:
        return None
    
    structure = fold_result['structure']
    
    # Calculate splice site positions in the structure
    # 5' splice site is at position 'padding'
    # 3' splice site is at position 'padding + intron_length - 1'
    ss5_pos = padding
    ss3_pos = padding + intron_length - 1
    
    # Calculate accessibility at splice sites
    ss5_acc_10nt = calculate_accessibility(structure, ss5_pos, window=10)
    ss5_acc_25nt = calculate_accessibility(structure, ss5_pos, window=25)
    ss5_acc_50nt = calculate_accessibility(structure, ss5_pos, window=50)
    
    ss3_acc_10nt = calculate_accessibility(structure, ss3_pos, window=10)
    ss3_acc_25nt = calculate_accessibility(structure, ss3_pos, window=25)
    ss3_acc_50nt = calculate_accessibility(structure, ss3_pos, window=50)
    
    # Check if splice sites are directly paired
    ss5_paired = None
    ss3_paired = None
    
    if 0 <= ss5_pos < len(structure):
        ss5_paired = structure[ss5_pos] in '(){}[]<>'
    
    if 0 <= ss3_pos < len(structure):
        ss3_paired = structure[ss3_pos] in '(){}[]<>'
    
    # Calculate overall structural properties
    total_bases = len(structure)
    unpaired_bases = structure.count('.')
    paired_bases = total_bases - unpaired_bases
    
    return {
        'mfe': fold_result['mfe'],
        'ensemble_energy': fold_result['ensemble_energy'],
        'structure_length': len(structure),
        'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence),
        'pairing_fraction': paired_bases / total_bases if total_bases > 0 else 0,
        'ss5_accessibility_10nt': ss5_acc_10nt,
        'ss5_accessibility_25nt': ss5_acc_25nt,
        'ss5_accessibility_50nt': ss5_acc_50nt,
        'ss5_paired': ss5_paired,
        'ss3_accessibility_10nt': ss3_acc_10nt,
        'ss3_accessibility_25nt': ss3_acc_25nt,
        'ss3_accessibility_50nt': ss3_acc_50nt,
        'ss3_paired': ss3_paired
    }

def process_intron_pair(row, fasta_file, padding=50):
    """
    Process one intron pair: fold each intron separately.
    
    Args:
        row: DataFrame row with intron pair information
        fasta_file: Path to reference genome FASTA
        padding: Amount of flanking sequence to include (bp)
    
    Returns:
        dict with features for both introns
    """
    chrom = row['chr']
    gene_id = row['gene_id']
    sample_id = row.get('sample_id', 'unknown')
    
    intron1_start = int(row['intron1_start'])
    intron1_end = int(row['intron1_end'])
    intron1_length = intron1_end - intron1_start
    
    intron2_start = int(row['intron2_start'])
    intron2_end = int(row['intron2_end'])
    intron2_length = intron2_end - intron2_start
    
    # Extract sequences with padding
    intron1_seq = extract_sequence(
        fasta_file, chrom, 
        max(0, intron1_start - padding), 
        intron1_end + padding
    )
    
    intron2_seq = extract_sequence(
        fasta_file, chrom,
        max(0, intron2_start - padding),
        intron2_end + padding
    )
    
    if not intron1_seq or not intron2_seq:
        return None
    
    # Fold each intron
    intron1_structure = analyze_intron_structure(intron1_seq, intron1_length, padding)
    print("intron1 strutcure ", intron1_structure)
    intron2_structure = analyze_intron_structure(intron2_seq, intron2_length, padding)
    
    if not intron1_structure or not intron2_structure:
        return None
    
    # Combine results
    result = {
        'sample_id': sample_id,
        'chr': chrom,
        'gene_id': gene_id,
        'intron1_start': intron1_start,
        'intron1_end': intron1_end,
        'intron1_length': intron1_length,
        'intron2_start': intron2_start,
        'intron2_end': intron2_end,
        'intron2_length': intron2_length,
        'upstream_count': int(row['upstream']),
        'downstream_count': int(row['downstream']),
        'total_reads': int(row['total']),
        'fraction_downstream': float(row['fraction_downstream'])
    }
    
    # Add intron1 features with prefix
    for key, value in intron1_structure.items():
        result[f'intron1_{key}'] = value
    
    # Add intron2 features with prefix
    for key, value in intron2_structure.items():
        result[f'intron2_{key}'] = value
    
    # Calculate comparative features
    if intron1_structure['mfe'] is not None and intron2_structure['mfe'] is not None:
        result['mfe_difference'] = intron1_structure['mfe'] - intron2_structure['mfe']
        result['mfe_ratio'] = intron1_structure['mfe'] / intron2_structure['mfe'] if intron2_structure['mfe'] != 0 else None
    else:
        result['mfe_difference'] = None
        result['mfe_ratio'] = None
    
    # Compare accessibility
    if intron1_structure['ss5_accessibility_25nt'] is not None and intron2_structure['ss5_accessibility_25nt'] is not None:
        result['ss5_accessibility_difference'] = (intron1_structure['ss5_accessibility_25nt'] - 
                                                   intron2_structure['ss5_accessibility_25nt'])
    else:
        result['ss5_accessibility_difference'] = None
    
    if intron1_structure['ss3_accessibility_25nt'] is not None and intron2_structure['ss3_accessibility_25nt'] is not None:
        result['ss3_accessibility_difference'] = (intron1_structure['ss3_accessibility_25nt'] - 
                                                   intron2_structure['ss3_accessibility_25nt'])
    else:
        result['ss3_accessibility_difference'] = None
    
    # Which intron is more structured? (FIX: convert bool to int properly)
    if intron1_structure['pairing_fraction'] is not None and intron2_structure['pairing_fraction'] is not None:
        result['intron1_more_structured'] = int(intron1_structure['pairing_fraction'] > intron2_structure['pairing_fraction'])
    else:
        result['intron1_more_structured'] = None
    
    # Which intron has more accessible splice sites? (FIX: same issue)
    intron1_avg_acc = None
    intron2_avg_acc = None
    
    if intron1_structure['ss5_accessibility_25nt'] is not None and intron1_structure['ss3_accessibility_25nt'] is not None:
        intron1_avg_acc = (intron1_structure['ss5_accessibility_25nt'] + intron1_structure['ss3_accessibility_25nt']) / 2
    
    if intron2_structure['ss5_accessibility_25nt'] is not None and intron2_structure['ss3_accessibility_25nt'] is not None:
        intron2_avg_acc = (intron2_structure['ss5_accessibility_25nt'] + intron2_structure['ss3_accessibility_25nt']) / 2
    
    if intron1_avg_acc is not None and intron2_avg_acc is not None:
        result['intron1_more_accessible'] = int(intron1_avg_acc > intron2_avg_acc)
    else:
        result['intron1_more_accessible'] = None

    print(result)
    
    return result

def main():
    parser = argparse.ArgumentParser(
        description='Extract RNA structure features for individual introns in pairs'
    )
    parser.add_argument('--intron-pairs', required=True,
                       help='TSV file with intron pairs (e.g., filtered pairs)')
    parser.add_argument('--fasta', required=True,
                       help='Reference genome FASTA (indexed)')
    parser.add_argument('--output', required=True,
                       help='Output TSV file')
    parser.add_argument('--padding', type=int, default=50,
                       help='Flanking sequence around each intron (bp, default: 50)')
    parser.add_argument('--max-pairs', type=int, default=None,
                       help='Maximum pairs to process (for testing)')
    
    args = parser.parse_args()
    
    # Check RNAfold
    try:
        subprocess.run(['RNAfold', '--version'], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("ERROR: RNAfold not found", file=sys.stderr)
        sys.exit(1)
    
    # Check FASTA index
    if not os.path.exists(f"{args.fasta}.fai"):
        print("ERROR: FASTA must be indexed", file=sys.stderr)
        sys.exit(1)
    
    print(f"Loading intron pairs from {args.intron_pairs}", file=sys.stderr)
    df = pd.read_csv(args.intron_pairs, sep='\t')
    
    if args.max_pairs:
        df = df.head(args.max_pairs)
    
    print(f"Processing {len(df)} intron pairs", file=sys.stderr)
    print(f"Padding: {args.padding} bp on each side", file=sys.stderr)
    
    all_results = []
    processed = 0
    failed = 0
    
    for idx, row in df.iterrows():
        result = process_intron_pair(row, args.fasta, args.padding)

        # print("Result: ", result)
        
        if result:
            all_results.append(result)
        else:
            failed += 1
        
        processed += 1
        if processed % 100 == 0:
            print(f"  Processed {processed}/{len(df)} pairs, "
                  f"success: {len(all_results)}, failed: {failed}", file=sys.stderr)
    
    print(f"\nProcessing complete:", file=sys.stderr)
    print(f"  Total pairs: {processed}", file=sys.stderr)
    print(f"  Successfully folded: {len(all_results)}", file=sys.stderr)
    print(f"  Failed: {failed}", file=sys.stderr)
    
    if all_results:
        result_df = pd.DataFrame(all_results)
        result_df.to_csv(args.output, sep='\t', index=False)
        
        print(f"\nResults saved to {args.output}", file=sys.stderr)
        print(f"Columns: {len(result_df.columns)}", file=sys.stderr)
        print(f"\nStructure features summary:", file=sys.stderr)
        print(f"  Mean intron1 MFE: {result_df['intron1_mfe'].mean():.2f}", file=sys.stderr)
        print(f"  Mean intron2 MFE: {result_df['intron2_mfe'].mean():.2f}", file=sys.stderr)
        print(f"  Mean MFE difference: {result_df['mfe_difference'].mean():.2f}", file=sys.stderr)
    else:
        print("ERROR: No results generated!", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()