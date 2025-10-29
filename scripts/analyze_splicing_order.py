#!/usr/bin/env python3
"""
Analyze splicing order from informative read pairs.

Input: BAM file where pairs are already filtered to be informative
       (one read has junction, one overlaps intron)

Output: For each intron pair in each gene, count:
        - upstream: intron1 spliced, intron2 retained
        - downstream: intron2 spliced, intron1 retained
"""

import pysam
import argparse
from collections import defaultdict
import sys
import gzip

def load_intron_annotations(intron_bed_file):
    """
    Load intron annotations from BED file.
    
    Returns:
        dict: {(chrom, start, end): gene_id}
    """
    introns = {}
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
                
                # Parse gene ID (adjust based on your BED format)
                if '_' in gene_info:
                    gene_id = gene_info.split('_')[0]
                elif '|' in gene_info:
                    gene_id = gene_info.split('|')[0]
                else:
                    gene_id = gene_info
                
                introns[(chrom, start, end)] = gene_id
    
    return introns

def get_spliced_introns(read, tolerance=10):
    """
    Extract intron coordinates that are SPLICED in this read (N in CIGAR).
    
    Returns:
        list of (chrom, start, end) tuples
    """
    if read.is_unmapped or not read.cigartuples:
        return []
    
    spliced = []
    ref_pos = read.reference_start
    
    for op, length in read.cigartuples:
        if op == 3:  # N operation = splice junction
            intron_start = ref_pos
            intron_end = ref_pos + length
            spliced.append((read.reference_name, intron_start, intron_end))
        
        if op in [0, 2, 3, 7, 8]:  # Operations consuming reference
            ref_pos += length
    
    return spliced

def get_overlapping_introns(read, introns, tolerance=10):
    """
    Find which annotated introns this read OVERLAPS (intron retained).
    
    Returns:
        list of (chrom, start, end, gene_id) tuples
    """
    if read.is_unmapped:
        return []
    
    overlapping = []
    read_start = read.reference_start
    read_end = read.reference_end
    
    for (chrom, intron_start, intron_end), gene_id in introns.items():
        if (chrom == read.reference_name and
            read_start < intron_end and read_end > intron_start):
            overlapping.append((chrom, intron_start, intron_end, gene_id))
    
    return overlapping

def coords_match(coord1, coord2, tolerance=10):
    """Check if two (chrom, start, end) tuples match within tolerance."""
    return (coord1[0] == coord2[0] and
            abs(coord1[1] - coord2[1]) <= tolerance and
            abs(coord1[2] - coord2[2]) <= tolerance)

def analyze_read_pair(read1, read2, introns, tolerance=10):
    """
    Analyze one informative read pair.
    
    Returns:
        list of (gene_id, intron1, intron2, direction) tuples
        where direction is 'upstream' or 'downstream'
    """
    results = []
    
    # Get spliced introns from both reads
    spliced1 = get_spliced_introns(read1, tolerance)
    spliced2 = get_spliced_introns(read2, tolerance)
    all_spliced = spliced1 + spliced2
    
    # Get retained introns (overlaps) from both reads
    retained1 = get_overlapping_introns(read1, introns, tolerance)
    retained2 = get_overlapping_introns(read2, introns, tolerance)
    all_retained = retained1 + retained2
    
    if not all_spliced or not all_retained:
        return results
    
    # Match spliced coordinates to annotated introns
    spliced_annotated = []
    for spliced_coord in all_spliced:
        for (chrom, start, end), gene_id in introns.items():
            if coords_match(spliced_coord, (chrom, start, end), tolerance):
                spliced_annotated.append((chrom, start, end, gene_id))
                break
    
    if not spliced_annotated:
        return results
    
    # For each (spliced, retained) pair, record evidence
    for spliced in spliced_annotated:
        for retained in all_retained:
            # Only consider pairs from same gene
            if spliced[3] != retained[3]:
                continue
            
            gene_id = spliced[3]
            
            # Determine which intron is upstream based on position
            if spliced[1] < retained[1]:  # spliced is upstream
                intron1 = (spliced[0], spliced[1], spliced[2])
                intron2 = (retained[0], retained[1], retained[2])
                direction = 'upstream'  # intron1 spliced first
            else:  # spliced is downstream
                intron1 = (retained[0], retained[1], retained[2])
                intron2 = (spliced[0], spliced[1], spliced[2])
                direction = 'downstream'  # intron2 spliced first
            
            results.append((gene_id, intron1, intron2, direction))
    
    return results

def main():
    parser = argparse.ArgumentParser(
        description='Analyze splicing order from informative read pairs'
    )
    parser.add_argument('--bam', required=True,
                       help='BAM file with informative read pairs (already filtered)')
    parser.add_argument('--intron-bed', required=True,
                       help='BED file with intron annotations')
    parser.add_argument('--output', required=True,
                       help='Output TSV file')
    parser.add_argument('--min-reads', type=int, default=5,
                       help='Minimum reads per intron pair (default: 5)')
    parser.add_argument('--min-mapq', type=int, default=10,
                       help='Minimum mapping quality (default: 10)')
    parser.add_argument('--tolerance', type=int, default=10,
                       help='Coordinate matching tolerance (bp, default: 10)')
    
    args = parser.parse_args()
    
    print(f"Loading intron annotations from {args.intron_bed}", file=sys.stderr)
    introns = load_intron_annotations(args.intron_bed)
    print(f"  Loaded {len(introns)} introns", file=sys.stderr)
    
    print(f"Loading informative read pairs from {args.bam}", file=sys.stderr)
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
    
    print(f"  Total reads in BAM: {total_reads}", file=sys.stderr)
    print(f"  Filtered (unmapped): {filtered_unmapped}", file=sys.stderr)
    print(f"  Filtered (MAPQ < {args.min_mapq}): {filtered_low_mapq}", file=sys.stderr)
    print(f"  Loaded {len(reads_by_name)} unique read pairs", file=sys.stderr)
    
    # Analyze each pair
    print("Analyzing read pairs...", file=sys.stderr)
    pair_counts = defaultdict(lambda: {'upstream': 0, 'downstream': 0})
    
    pairs_analyzed = 0
    evidence_found = 0
    
    for query_name, reads in reads_by_name.items():
        if len(reads) != 2:
            continue
        
        pairs_analyzed += 1
        if pairs_analyzed % 10000 == 0:
            print(f"  Processed {pairs_analyzed} pairs, {evidence_found} evidence events...",
                  file=sys.stderr)
        
        read1, read2 = reads
        results = analyze_read_pair(read1, read2, introns, args.tolerance)
        
        for gene_id, intron1, intron2, direction in results:
            evidence_found += 1
            # Create unique key for this intron pair
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
        
        if valid_pairs == 0:
            print(f"\nWARNING: No intron pairs met minimum threshold!",
                  file=sys.stderr)
            print(f"  Try reducing --min-reads (currently {args.min_reads})",
                  file=sys.stderr)

if __name__ == "__main__":
    main()