#!/usr/bin/env python3
"""
Filter splicing order results to find intron pairs where both introns are small enough to fold.
Creates a comprehensive list across all samples.
"""

import pandas as pd
import argparse
import sys
from pathlib import Path
import glob

def filter_intron_pairs(df, max_intron_size=1000, max_region_size=5000):
    """
    Filter intron pairs based on size criteria.
    
    Args:
        df: DataFrame with splicing order data
        max_intron_size: Maximum size for each individual intron
        max_region_size: Maximum size for the region spanning both introns
    
    Returns:
        Filtered DataFrame
    """
    # Calculate intron sizes
    df['intron1_length'] = df['intron1_end'] - df['intron1_start']
    df['intron2_length'] = df['intron2_end'] - df['intron2_start']
    
    # Calculate region size (from start of first intron to end of second)
    df['region_start'] = df[['intron1_start', 'intron2_start']].min(axis=1)
    df['region_end'] = df[['intron1_end', 'intron2_end']].max(axis=1)
    df['region_length'] = df['region_end'] - df['region_start']
    
    # Filter based on criteria
    filtered = df[
        (df['intron1_length'] <= max_intron_size) &
        (df['intron2_length'] <= max_intron_size) &
        (df['region_length'] <= max_region_size)
    ].copy()
    
    return filtered

def process_sample_file(file_path, max_intron_size, max_region_size, min_reads):
    """Process a single splicing order file."""
    try:
        df = pd.read_csv(file_path, sep='\t')
        
        # Add sample ID from filename
        sample_id = Path(file_path).stem.replace('_pairwise_splicing_order', '')
        df['sample_id'] = sample_id
        
        # Filter by minimum read count
        df = df[df['total'] >= min_reads]
        
        # Filter by size
        filtered = filter_intron_pairs(df, max_intron_size, max_region_size)
        
        return filtered, len(df), len(filtered)
    
    except Exception as e:
        print(f"Error processing {file_path}: {e}", file=sys.stderr)
        return None, 0, 0

def main():
    parser = argparse.ArgumentParser(
        description='Filter splicing order results for foldable intron pairs'
    )
    parser.add_argument('--results-dir', required=True,
                       help='Directory containing sample result folders')
    parser.add_argument('--output', required=True,
                       help='Output TSV file with filtered intron pairs')
    parser.add_argument('--max-intron-size', type=int, default=1000,
                       help='Maximum size for each intron (bp, default: 1000)')
    parser.add_argument('--max-region-size', type=int, default=5000,
                       help='Maximum size for region spanning both introns (bp, default: 5000)')
    parser.add_argument('--min-reads', type=int, default=5,
                       help='Minimum number of supporting reads (default: 5)')
    parser.add_argument('--pattern', default='*_pairwise_splicing_order.tsv',
                       help='File pattern to match (default: *_pairwise_splicing_order.tsv)')
    
    args = parser.parse_args()
    
    # Find all splicing order files
    search_pattern = f"{args.results_dir}/**/{args.pattern}"
    files = glob.glob(search_pattern, recursive=True)
    
    print(f"Found {len(files)} splicing order files", file=sys.stderr)
    
    if len(files) == 0:
        print(f"ERROR: No files found matching pattern: {search_pattern}", file=sys.stderr)
        sys.exit(1)
    
    # Process all files
    all_filtered = []
    total_pairs = 0
    total_filtered = 0
    samples_processed = 0
    samples_with_pairs = 0
    
    for file_path in sorted(files):
        sample_name = Path(file_path).parent.name
        print(f"Processing {sample_name}...", file=sys.stderr)
        
        filtered, n_total, n_filtered = process_sample_file(
            file_path, args.max_intron_size, args.max_region_size, args.min_reads
        )
        
        if filtered is not None:
            samples_processed += 1
            total_pairs += n_total
            total_filtered += n_filtered
            
            if n_filtered > 0:
                all_filtered.append(filtered)
                samples_with_pairs += 1
                print(f"  {n_filtered}/{n_total} pairs passed filter", file=sys.stderr)
            else:
                print(f"  0/{n_total} pairs passed filter", file=sys.stderr)
    
    # Combine all filtered results
    if all_filtered:
        combined_df = pd.concat(all_filtered, ignore_index=True)
        
        # Sort by chromosome and position
        combined_df = combined_df.sort_values(['chr', 'intron1_start'])
        
        # Save to file
        combined_df.to_csv(args.output, sep='\t', index=False)
        
        print(f"\n{'='*80}", file=sys.stderr)
        print(f"SUMMARY", file=sys.stderr)
        print(f"{'='*80}", file=sys.stderr)
        print(f"Samples processed: {samples_processed}", file=sys.stderr)
        print(f"Samples with filtered pairs: {samples_with_pairs}", file=sys.stderr)
        print(f"Total intron pairs: {total_pairs}", file=sys.stderr)
        print(f"Filtered pairs: {total_filtered} ({100*total_filtered/total_pairs:.1f}%)", file=sys.stderr)
        print(f"\nOutput saved to: {args.output}", file=sys.stderr)
        print(f"Unique genes: {combined_df['gene_id'].nunique()}", file=sys.stderr)
        print(f"Unique intron pairs: {len(combined_df)}", file=sys.stderr)
        
        # Show distribution by chromosome
        print(f"\nPairs by chromosome:", file=sys.stderr)
        print(combined_df['chr'].value_counts().head(10).to_string(), file=sys.stderr)
        
        # Show size distribution
        print(f"\nIntron size distribution:", file=sys.stderr)
        print(f"  Mean intron1 length: {combined_df['intron1_length'].mean():.0f} bp", file=sys.stderr)
        print(f"  Mean intron2 length: {combined_df['intron2_length'].mean():.0f} bp", file=sys.stderr)
        print(f"  Mean region length: {combined_df['region_length'].mean():.0f} bp", file=sys.stderr)
        
    else:
        print(f"\nERROR: No intron pairs passed filters!", file=sys.stderr)
        print(f"Try relaxing the filters:", file=sys.stderr)
        print(f"  --max-intron-size {args.max_intron_size * 2}", file=sys.stderr)
        print(f"  --max-region-size {args.max_region_size * 2}", file=sys.stderr)
        print(f"  --min-reads {max(1, args.min_reads - 2)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()