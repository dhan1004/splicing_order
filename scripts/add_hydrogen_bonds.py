#!/usr/bin/env python3
"""
Add hydrogen bond counts to structure features files.
Extracts sequences from genome, folds with RNAfold, and counts H-bonds.
"""

import pandas as pd
import argparse
import sys
from pathlib import Path
import subprocess
import pysam

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
        
        if len(lines) < 3:
            return None
        
        structure_line = lines[2].strip()
        
        if '.' not in structure_line or '(' not in structure_line:
            return None
        
        parts = structure_line.split()
        
        if len(parts) < 1:
            return None
        
        structure = parts[0]
        
        # Extract MFE
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
        
        # Extract ensemble energy
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

def count_hydrogen_bonds(sequence, structure):
    """
    Count hydrogen bonds in an RNA structure.
    
    Watson-Crick pairs:
    - G-C: 3 hydrogen bonds
    - A-U: 2 hydrogen bonds
    
    Wobble pairs:
    - G-U: 2 hydrogen bonds (weaker)
    """
    # Parse structure to find base pairs
    stack = []
    pairs = {}
    
    for i, char in enumerate(structure):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                j = stack.pop()
                pairs[j] = i
                pairs[i] = j
    
    # Count hydrogen bonds
    total_hbonds = 0
    gc_pairs = 0
    au_pairs = 0
    gu_pairs = 0
    other_pairs = 0
    
    for i, j in pairs.items():
        if i < j:  # Count each pair once
            base1 = sequence[i]
            base2 = sequence[j]
            
            # Determine hydrogen bond count
            pair = ''.join(sorted([base1, base2]))
            
            if pair == 'CG':
                hbonds = 3
                gc_pairs += 1
            elif pair == 'AU':
                hbonds = 2
                au_pairs += 1
            elif pair == 'GU':
                hbonds = 2  # Wobble pair
                gu_pairs += 1
            else:
                hbonds = 0  # Non-canonical
                other_pairs += 1
            
            total_hbonds += hbonds
    
    return {
        'total_hbonds': total_hbonds,
        'gc_pairs': gc_pairs,
        'au_pairs': au_pairs,
        'gu_pairs': gu_pairs,
        'other_pairs': other_pairs,
        'total_pairs': len(pairs) // 2
    }

def process_intron_with_hbonds(chrom, intron_start, intron_end, fasta_file, padding=50):
    """
    Extract intron sequence, fold it, and count hydrogen bonds.
    
    Args:
        chrom: Chromosome
        intron_start: Intron start position
        intron_end: Intron end position
        fasta_file: Reference genome FASTA
        padding: Flanking sequence (bp)
    
    Returns:
        dict with sequence, structure, and H-bond counts
    """
    # Extract sequence with padding
    sequence = extract_sequence(
        fasta_file, chrom,
        max(0, intron_start - padding),
        intron_end + padding
    )
    
    if not sequence:
        return None
    
    # Fold the sequence
    fold_result = fold_rna_with_rnafold(sequence)
    
    if not fold_result:
        return None
    
    structure = fold_result['structure']
    
    # Count hydrogen bonds
    hbond_stats = count_hydrogen_bonds(sequence, structure)
    
    # Combine results
    return {
        'sequence': sequence,
        'structure': structure,
        'mfe': fold_result['mfe'],
        'ensemble_energy': fold_result['ensemble_energy'],
        **hbond_stats
    }

def process_structure_file(input_file, output_file, fasta_file, padding=50):
    """
    Process a structure features file and add hydrogen bond counts.
    
    Args:
        input_file: Path to existing structure_features.tsv
        output_file: Path to output file with H-bonds
        fasta_file: Reference genome FASTA
        padding: Flanking sequence around introns (bp)
    """
    print(f"Reading {input_file}...", file=sys.stderr)
    
    try:
        df = pd.read_csv(input_file, sep='\t')
    except Exception as e:
        print(f"Error reading {input_file}: {e}", file=sys.stderr)
        return False
    
    print(f"Loaded {len(df)} intron pairs", file=sys.stderr)
    
    # Initialize new columns
    for suffix in ['1', '2']:
        df[f'intron{suffix}_sequence'] = ''
        df[f'intron{suffix}_structure'] = ''
        df[f'intron{suffix}_total_hbonds'] = 0
        df[f'intron{suffix}_gc_pairs'] = 0
        df[f'intron{suffix}_au_pairs'] = 0
        df[f'intron{suffix}_gu_pairs'] = 0
        df[f'intron{suffix}_other_pairs'] = 0
        df[f'intron{suffix}_total_pairs'] = 0
    
    # Process each row
    for idx, row in df.iterrows():
        if idx % 100 == 0:
            print(f"Processing row {idx}/{len(df)}...", file=sys.stderr)
        
        chrom = row['chr']
        
        # Process intron 1
        intron1_result = process_intron_with_hbonds(
            chrom, 
            int(row['intron1_start']), 
            int(row['intron1_end']),
            fasta_file,
            padding
        )
        
        if intron1_result:
            df.at[idx, 'intron1_sequence'] = intron1_result['sequence']
            df.at[idx, 'intron1_structure'] = intron1_result['structure']
            df.at[idx, 'intron1_total_hbonds'] = intron1_result['total_hbonds']
            df.at[idx, 'intron1_gc_pairs'] = intron1_result['gc_pairs']
            df.at[idx, 'intron1_au_pairs'] = intron1_result['au_pairs']
            df.at[idx, 'intron1_gu_pairs'] = intron1_result['gu_pairs']
            df.at[idx, 'intron1_other_pairs'] = intron1_result['other_pairs']
            df.at[idx, 'intron1_total_pairs'] = intron1_result['total_pairs']
        
        # Process intron 2
        intron2_result = process_intron_with_hbonds(
            chrom,
            int(row['intron2_start']),
            int(row['intron2_end']),
            fasta_file,
            padding
        )
        
        if intron2_result:
            df.at[idx, 'intron2_sequence'] = intron2_result['sequence']
            df.at[idx, 'intron2_structure'] = intron2_result['structure']
            df.at[idx, 'intron2_total_hbonds'] = intron2_result['total_hbonds']
            df.at[idx, 'intron2_gc_pairs'] = intron2_result['gc_pairs']
            df.at[idx, 'intron2_au_pairs'] = intron2_result['au_pairs']
            df.at[idx, 'intron2_gu_pairs'] = intron2_result['gu_pairs']
            df.at[idx, 'intron2_other_pairs'] = intron2_result['other_pairs']
            df.at[idx, 'intron2_total_pairs'] = intron2_result['total_pairs']
    
    # Calculate derived metrics
    df['hbond_ratio'] = df['intron1_total_hbonds'] / (df['intron2_total_hbonds'] + 1)
    df['hbond_difference'] = df['intron1_total_hbonds'] - df['intron2_total_hbonds']
    df['gc_pair_fraction_diff'] = (df['intron1_gc_pairs'] / (df['intron1_total_pairs'] + 1)) - \
                                   (df['intron2_gc_pairs'] / (df['intron2_total_pairs'] + 1))
    
    # Save output
    print(f"Writing to {output_file}...", file=sys.stderr)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Done! Processed {len(df)} intron pairs", file=sys.stderr)
    
    # Print summary statistics
    print("\n=== Hydrogen Bond Statistics ===", file=sys.stderr)
    for intron_num in [1, 2]:
        mean_hbonds = df[f'intron{intron_num}_total_hbonds'].mean()
        mean_gc = df[f'intron{intron_num}_gc_pairs'].mean()
        mean_au = df[f'intron{intron_num}_au_pairs'].mean()
        mean_gu = df[f'intron{intron_num}_gu_pairs'].mean()
        print(f"Intron {intron_num}:", file=sys.stderr)
        print(f"  Mean H-bonds: {mean_hbonds:.1f}", file=sys.stderr)
        print(f"  Mean GC pairs: {mean_gc:.1f}", file=sys.stderr)
        print(f"  Mean AU pairs: {mean_au:.1f}", file=sys.stderr)
        print(f"  Mean GU pairs: {mean_gu:.1f}", file=sys.stderr)
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description='Add hydrogen bond counts to structure features files'
    )
    parser.add_argument('--input', '-i', required=True,
                        help='Input structure_features.tsv file')
    parser.add_argument('--output', '-o', required=True,
                        help='Output file with H-bonds and sequences added')
    parser.add_argument('--fasta', '-f', required=True,
                        help='Reference genome FASTA file (indexed)')
    parser.add_argument('--padding', type=int, default=50,
                        help='Flanking sequence around introns (bp, default: 50)')
    
    args = parser.parse_args()
    
    # Check input exists
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Check FASTA exists and is indexed
    if not Path(args.fasta).exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)
    
    if not Path(f"{args.fasta}.fai").exists():
        print(f"Error: FASTA file must be indexed. Run: samtools faidx {args.fasta}",
              file=sys.stderr)
        sys.exit(1)
    
    # Check RNAfold
    try:
        subprocess.run(['RNAfold', '--version'], capture_output=True, check=True)
    except Exception:
        print("Error: RNAfold not found. Make sure ViennaRNA is installed.",
              file=sys.stderr)
        sys.exit(1)
    
    # Process the file
    success = process_structure_file(args.input, args.output, args.fasta, args.padding)
    
    if not success:
        sys.exit(1)

if __name__ == '__main__':
    main()