#!/usr/bin/env python3
"""
Minimal test: Extract one sequence and test RNAfold
"""

import pysam
import subprocess
import tempfile
import os
import sys

def extract_sequence(fasta_file, chrom, start, end):
    """Extract genomic sequence from fasta file."""
    with pysam.FastaFile(fasta_file) as fasta:
        seq = fasta.fetch(chrom, start, end)
        seq = seq.upper().replace('T', 'U')
    return seq

def test_rnafold(sequence, label="Test"):
    """Test RNAfold on a sequence"""
    print(f"\n{'='*60}")
    print(f"{label}")
    print(f"{'='*60}")
    print(f"Sequence length: {len(sequence)}")
    print(f"First 80bp: {sequence[:80]}")
    print(f"Last 80bp: {sequence[-80:]}")
    
    # Check for invalid characters
    valid_chars = set('AUCGN')
    seq_chars = set(sequence)
    invalid = seq_chars - valid_chars
    if invalid:
        print(f"⚠️  Invalid characters found: {invalid}")
    
    # Create temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False) as f:
        f.write(f">seq\n{sequence}\n")
        temp_file = f.name
    
    print(f"\nRunning RNAfold...")
    
    try:
        with open(temp_file, 'r') as f_in:
            result = subprocess.run(
                ['RNAfold', '--noPS', '-p'],
                stdin=f_in,
                capture_output=True,
                text=True,
                timeout=30
            )
        
        print(f"Return code: {result.returncode}")
        print(f"\nSTDOUT ({len(result.stdout)} chars):")
        print(result.stdout)
        print(f"\nSTDERR ({len(result.stderr)} chars):")
        print(result.stderr)
        
        # Try to parse
        lines = result.stdout.strip().split('\n')
        print(f"\nNumber of output lines: {len(lines)}")
        
        for i, line in enumerate(lines):
            print(f"Line {i}: {repr(line[:100])}{'...' if len(line) > 100 else ''}")
        
        if len(lines) >= 2:
            print(f"\n✅ RNAfold produced expected output")
            structure_line = lines[1]
            
            if '(' in structure_line:
                parts = structure_line.rsplit('(', 1)
                if len(parts) == 2:
                    structure = parts[0].strip()
                    mfe_str = parts[1].strip().rstrip(')')
                    try:
                        mfe = float(mfe_str)
                        print(f"✅ Successfully parsed: MFE = {mfe}")
                        print(f"   Structure length: {len(structure)}")
                        return True
                    except ValueError as e:
                        print(f"❌ Cannot parse MFE: {e}")
                        print(f"   MFE string was: '{mfe_str}'")
                        return False
        else:
            print(f"❌ Not enough output lines")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"❌ RNAfold TIMEOUT")
        return False
    except Exception as e:
        print(f"❌ Exception: {e}")
        import traceback
        traceback.print_exc()
        return False
    finally:
        os.unlink(temp_file)

def main():
    if len(sys.argv) < 3:
        print("Usage: python3 test_real_sequence.py <bam_file> <fasta_file>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    fasta_file = sys.argv[2]
    
    print("="*60)
    print("Testing RNAfold on REAL sequences from your BAM")
    print("="*60)
    
    # Get first read pair from BAM
    print("\nReading BAM file...")
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        reads_by_name = {}
        
        for read in bam:
            if read.is_unmapped:
                continue
            
            name = read.query_name
            if name not in reads_by_name:
                reads_by_name[name] = []
            reads_by_name[name].append(read)
            
            # Stop after finding first pair
            if len(reads_by_name) > 0:
                first_pair = list(reads_by_name.values())[0]
                if len(first_pair) == 2:
                    break
    
    if not reads_by_name or len(list(reads_by_name.values())[0]) < 2:
        print("ERROR: Could not find a read pair in BAM")
        sys.exit(1)
    
    read1, read2 = list(reads_by_name.values())[0]
    
    print(f"\nFound read pair: {read1.query_name}")
    print(f"  Read1: {read1.reference_name}:{read1.reference_start}-{read1.reference_end}")
    print(f"  Read2: {read2.reference_name}:{read2.reference_start}-{read2.reference_end}")
    
    # Test with different window sizes
    flank_size = 100
    
    for read_idx, read in enumerate([read1, read2], 1):
        print(f"\n{'='*60}")
        print(f"TESTING READ {read_idx}")
        print(f"{'='*60}")
        
        start = max(0, read.reference_start - flank_size)
        end = read.reference_end + flank_size
        
        print(f"Extracting: {read.reference_name}:{start}-{end}")
        
        try:
            sequence = extract_sequence(fasta_file, read.reference_name, start, end)
            print(f"✅ Extracted {len(sequence)}bp")
        except Exception as e:
            print(f"❌ Failed to extract sequence: {e}")
            continue
        
        # Test full sequence
        success = test_rnafold(sequence, f"Read{read_idx} - Full sequence ({len(sequence)}bp)")
        
        if not success:
            print(f"\n⚠️  Full sequence RNAfold FAILED")
            
            # Try shorter sequences
            print(f"\nTrying shorter sequences...")
            for test_len in [100, 50, 30]:
                if len(sequence) >= test_len:
                    test_seq = sequence[:test_len]
                    print(f"\n--- Testing first {test_len}bp ---")
                    test_rnafold(test_seq, f"First {test_len}bp")
        
        # Test window sizes from your script
        for window_size in [50, 100, 200]:
            if len(sequence) >= window_size:
                window_seq = sequence[:window_size]
                test_rnafold(window_seq, f"Read{read_idx} - Window {window_size}bp")
    
    print("\n" + "="*60)
    print("TEST COMPLETE")
    print("="*60)

if __name__ == "__main__":
    main()

python3 /users/dhan30/scratch/splicing_order/scripts/analyze_splicing_order.py \
    --bam /users/dhan30/scratch/splicing_order/results/GSM6754877/GSM6754877_informative_pairs.bam \
    --intron-bed /users/dhan30/reference/hg38.gencode.basic.v43.introns.bed.gz \
    --output pairwise_splicing_order.tsv \
    --min-reads 5 \
    --min-mapq 10 \
    --tolerance 10