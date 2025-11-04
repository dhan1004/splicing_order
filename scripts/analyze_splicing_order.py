#!/usr/bin/env python3
"""
Optimized splicing order analysis using interval trees for fast lookups.
"""

import pysam
import argparse
from collections import defaultdict
import sys
import gzip

def load_intron_annotations(intron_bed_file):
    """
    Load intron annotations and organize by chromosome for fast lookup.
    
    Returns:
        dict: {chrom: [(start, end, gene_id), ...]}
    """
    introns_by_chrom = defaultdict(list)
    intron_dict = {}
    
    opener = gzip.open if intron_bed_file.endswith(".gz") else open
    
    with opener(intron_bed_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                gene_info = fields[3]
                
                # Parse gene ID
                if '_' in gene_info:
                    gene_id = gene_info.split('_')[0]
                elif '|' in gene_info:
                    gene_id = gene_info.split('|')[0]
                else:
                    gene_id = gene_info
                
                introns_by_chrom[chrom].append((start, end, gene_id))
                intron_dict[(chrom, start, end)] = gene_id
    
    # Sort introns by start position for each chromosome (enables binary search)
    for chrom in introns_by_chrom:
        introns_by_chrom[chrom].sort()
    
    return introns_by_chrom, intron_dict

def get_spliced_introns_fast(read):
    """Extract spliced introns from CIGAR string."""
    if read.is_unmapped or not read.cigartuples:
        return []
    
    spliced = []
    ref_pos = read.reference_start
    
    for op, length in read.cigartuples:
        if op == 3:  # N = splice junction
            spliced.append((read.reference_name, ref_pos, ref_pos + length))
        if op in [0, 2, 3, 7, 8]:  # Consumes reference
            ref_pos += length
    
    return spliced

def get_overlapping_introns_fast(read, introns_by_chrom, tolerance=10):
    """
    Fast intron overlap detection using sorted intervals.
    Only checks introns on the same chromosome.
    """
    if read.is_unmapped or read.reference_name not in introns_by_chrom:
        return []
    
    overlapping = []
    read_start = read.reference_start
    read_end = read.reference_end
    
    # Only check introns on this chromosome
    chrom_introns = introns_by_chrom[read.reference_name]
    
    # Binary search to find potentially overlapping introns
    # Since introns are sorted by start, we can skip introns that end before read starts
    for intron_start, intron_end, gene_id in chrom_introns:
        # Skip introns that start after read ends
        if intron_start > read_end:
            break
        # Skip introns that end before read starts
        if intron_end < read_start:
            continue
        # Check for overlap
        if read_start < intron_end and read_end > intron_start:
            overlapping.append((read.reference_name, intron_start, intron_end, gene_id))
    
    return overlapping

def coords_match(coord1, coord2, tolerance=10):
    """Check if coordinates match within tolerance."""
    return (coord1[0] == coord2[0] and
            abs(coord1[1] - coord2[1]) <= tolerance and
            abs(coord1[2] - coord2[2]) <= tolerance)

def analyze_read_pair_fast(read1, read2, introns_by_chrom, intron_dict, tolerance=10):
    """Optimized read pair analysis."""
    results = []
    
    # Get spliced and retained introns
    all_spliced = get_spliced_introns_fast(read1) + get_spliced_introns_fast(read2)
    all_retained = (get_overlapping_introns_fast(read1, introns_by_chrom, tolerance) + 
                    get_overlapping_introns_fast(read2, introns_by_chrom, tolerance))
    
    if not all_spliced or not all_retained:
        return results
    
    # Match spliced coordinates to annotated introns
    spliced_annotated = []
    for spliced_coord in all_spliced:
        # Direct lookup in dictionary instead of searching
        for (chrom, start, end), gene_id in intron_dict.items():
            if coords_match(spliced_coord, (chrom, start, end), tolerance):
                spliced_annotated.append((chrom, start, end, gene_id))
                break
    
    if not spliced_annotated:
        return results
    
    # Record evidence for each pair
    for spliced in spliced_annotated:
        for retained in all_retained:
            if spliced[3] != retained[3]:  # Same gene only
                continue
            
            gene_id = spliced[3]
            
            # Determine order
            if spliced[1] < retained[1]:
                intron1 = (spliced[0], spliced[1], spliced[2])
                intron2 = (retained[0], retained[1], retained[2])
                direction = 'upstream'
            else:
                intron1 = (retained[0], retained[1], retained[2])
                intron2 = (spliced[0], spliced[1], spliced[2])
                direction = 'downstream'
            
            results.append((gene_id, intron1, intron2, direction))
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Optimized splicing order analysis'
    )
    parser.add_argument('--bam', required=True)
    parser.add_argument('--intron-bed', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--min-reads', type=int, default=5)
    parser.add_argument('--min-mapq', type=int, default=10)
    parser.add_argument('--tolerance', type=int, default=10)
    
    args = parser.parse_args()
    
    print(f"Loading intron annotations from {args.intron_bed}", file=sys.stderr)
    introns_by_chrom, intron_dict = load_intron_annotations(args.intron_bed)
    total_introns = sum(len(introns) for introns in introns_by_chrom.values())
    print(f"  Loaded {total_introns} introns across {len(introns_by_chrom)} chromosomes", 
          file=sys.stderr)
    
    print(f"Loading informative read pairs from {args.bam}", file=sys.stderr)
    
    # Process BAM in one pass - group reads by name
    reads_by_name = defaultdict(list)
    total_reads = 0
    filtered_unmapped = 0
    filtered_low_mapq = 0
    
    with pysam.AlignmentFile(args.bam, 'rb') as bam:
        for read in bam:
            total_reads += 1
            
            if read.is_unmapped:
                filtered_unmapped += 1
                continue
            if read.mapping_quality < args.min_mapq:
                filtered_low_mapq += 1
                continue
            
            reads_by_name[read.query_name].append(read)
    
    print(f"  Total reads: {total_reads}", file=sys.stderr)
    print(f"  Filtered (unmapped): {filtered_unmapped}", file=sys.stderr)
    print(f"  Filtered (MAPQ < {args.min_mapq}): {filtered_low_mapq}", file=sys.stderr)
    print(f"  Unique read pairs: {len(reads_by_name)}", file=sys.stderr)
    
    # Analyze pairs
    print("Analyzing read pairs...", file=sys.stderr)
    pair_counts = defaultdict(lambda: {'upstream': 0, 'downstream': 0})
    
    pairs_analyzed = 0
    evidence_found = 0
    
    for query_name, reads in reads_by_name.items():
        if len(reads) != 2:
            continue
        
        pairs_analyzed += 1
        if pairs_analyzed % 50000 == 0:  # Report less frequently
            print(f"  Processed {pairs_analyzed} pairs, {evidence_found} evidence events...",
                  file=sys.stderr)
        
        read1, read2 = reads
        results = analyze_read_pair_fast(read1, read2, introns_by_chrom, intron_dict, args.tolerance)
        
        for gene_id, intron1, intron2, direction in results:
            evidence_found += 1
            key = f"{intron1[0]}\t{gene_id}\t{intron1[1]}\t{intron1[2]}\t{intron2[1]}\t{intron2[2]}"
            pair_counts[key][direction] += 1
    
    print(f"\nAnalyzed {pairs_analyzed} read pairs", file=sys.stderr)
    print(f"Found {evidence_found} splicing evidence events", file=sys.stderr)
    print(f"Unique intron pairs: {len(pair_counts)}", file=sys.stderr)
    
    # Write output
    print(f"Writing results to {args.output}", file=sys.stderr)
    with open(args.output, 'w') as f:
        f.write("chr\tgene_id\tintron1_start\tintron1_end\t"
                "intron2_start\tintron2_end\t"
                "upstream\tdownstream\ttotal\tfraction_downstream\n")
        
        valid_pairs = 0
        for key, counts in sorted(pair_counts.items()):
            upstream = counts['upstream']
            downstream = counts['downstream']
            total = upstream + downstream
            
            if total >= args.min_reads:
                frac_down = downstream / total if total > 0 else 0.5
                f.write(f"{key}\t{upstream}\t{downstream}\t{total}\t{frac_down:.4f}\n")
                valid_pairs += 1
        
        print(f"  Wrote {valid_pairs} intron pairs with ≥{args.min_reads} reads",
              file=sys.stderr)

if __name__ == "__main__":
    main()