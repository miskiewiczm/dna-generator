#!/usr/bin/env python3
"""
Proof-of-concept: Exhaustive search dla 15nt sekwencji.

4^15 = 1,073,741,824 (~1 billion sequences)
Expected time: ~5-6 hours na desktop (8 cores)
Memory: ~100 GB

To jest realny POC pokazujący że exhaustive search jest WYKONALNY.
"""

import multiprocessing as mp
from itertools import product
import time
from collections import defaultdict
import sys

# Simplified validation (same as benchmark)
def simple_gc_content(seq):
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq)

def check_homopolymers(seq, max_len=4):
    for base in 'ATGC':
        if base * (max_len + 1) in seq:
            return False
    return True

def check_dinucleotide_repeats(seq, max_repeats=3):
    for i in range(4):
        for j in range(4):
            bases = 'ATGC'
            dinuc = bases[i] + bases[j]
            repeat = dinuc * (max_repeats + 1)
            if repeat in seq:
                return False
    return True

def encode_sequence(seq):
    """Encode sequence to integer (2 bits per base)."""
    mapping = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
    result = 0
    for base in seq:
        result = (result << 2) | mapping[base]
    return result

def decode_sequence(encoded, length=15):
    """Decode integer back to sequence."""
    mapping = ['A', 'T', 'G', 'C']
    seq = []
    for _ in range(length):
        seq.append(mapping[encoded & 3])
        encoded >>= 2
    return ''.join(reversed(seq))

def validate_fast(seq, min_gc=0.35, max_gc=0.65):
    """Fast validation with multi-level pruning."""

    # Level 1: Ultra-fast checks (~0.5 µs)
    gc = simple_gc_content(seq)
    if gc < min_gc or gc > max_gc:
        return False, 'gc_content'

    # Level 2: Fast checks (~1 µs)
    if not check_homopolymers(seq, max_len=4):
        return False, 'homopolymer'

    if not check_dinucleotide_repeats(seq, max_repeats=3):
        return False, 'dinucleotide'

    # Passed all checks
    return True, None

def process_chunk(args):
    """Process a chunk of sequences (worker function)."""
    start_idx, end_idx, length = args

    valid_count = 0
    fail_reasons = defaultdict(int)
    chunk_valid = []

    # Generate sequences in this range
    for i in range(start_idx, end_idx):
        seq = decode_sequence(i, length)
        is_valid, fail_reason = validate_fast(seq)

        if is_valid:
            valid_count += 1
            # Store encoded (saves memory)
            chunk_valid.append(i)
        else:
            fail_reasons[fail_reason] += 1

    return valid_count, dict(fail_reasons), chunk_valid

def exhaustive_search_parallel(length=15, num_workers=8):
    """Exhaustive search with parallel processing."""

    total_sequences = 4 ** length
    chunk_size = total_sequences // (num_workers * 4)  # 4 chunks per worker

    print("=" * 70)
    print(f"Exhaustive search dla {length}nt sekwencji")
    print("=" * 70)
    print(f"\nTotal sequences: {total_sequences:,}")
    print(f"Workers: {num_workers}")
    print(f"Chunk size: {chunk_size:,}")
    print(f"Total chunks: {total_sequences // chunk_size}")

    # Prepare chunks
    chunks = []
    for i in range(0, total_sequences, chunk_size):
        end = min(i + chunk_size, total_sequences)
        chunks.append((i, end, length))

    print(f"\nStarting parallel processing...")
    print(f"Progress will be reported every 10 chunks\n")

    start_time = time.time()

    # Process in parallel
    total_valid = 0
    total_fail_reasons = defaultdict(int)
    all_valid_sequences = []

    with mp.Pool(processes=num_workers) as pool:
        # Process chunks
        for chunk_idx, result in enumerate(pool.imap_unordered(process_chunk, chunks)):
            valid_count, fail_reasons, valid_seqs = result

            total_valid += valid_count
            all_valid_sequences.extend(valid_seqs)

            for reason, count in fail_reasons.items():
                total_fail_reasons[reason] += count

            # Progress reporting
            if (chunk_idx + 1) % 10 == 0 or (chunk_idx + 1) == len(chunks):
                processed = (chunk_idx + 1) * chunk_size
                elapsed = time.time() - start_time
                rate = processed / elapsed
                remaining = (total_sequences - processed) / rate if rate > 0 else 0

                print(f"Progress: {processed:,}/{total_sequences:,} "
                      f"({processed/total_sequences*100:.1f}%) | "
                      f"Valid so far: {total_valid:,} | "
                      f"Rate: {rate:,.0f} seq/s | "
                      f"ETA: {remaining/60:.1f} min")

    end_time = time.time()
    total_time = end_time - start_time

    # Results
    print("\n" + "=" * 70)
    print("WYNIKI")
    print("=" * 70)

    print(f"\nCzas wykonania: {total_time:.1f} seconds = {total_time/60:.1f} minutes = {total_time/3600:.2f} hours")
    print(f"Throughput: {total_sequences/total_time:,.0f} seq/sec")

    print(f"\nStatystyki walidacji:")
    print(f"  Total sequences: {total_sequences:,}")
    print(f"  Valid sequences: {total_valid:,} ({total_valid/total_sequences*100:.2f}%)")
    print(f"  Invalid sequences: {total_sequences - total_valid:,}")

    print(f"\nPowody odrzucenia:")
    for reason, count in sorted(total_fail_reasons.items(), key=lambda x: -x[1]):
        print(f"  {reason}: {count:,} ({count/total_sequences*100:.2f}%)")

    # Memory analysis
    print(f"\nPamięć słownika:")

    # Encoded format (4 bytes per sequence)
    bytes_encoded = len(all_valid_sequences) * 4
    print(f"  Encoded (int32): {bytes_encoded:,} bytes = {bytes_encoded/(1024**2):.1f} MB")

    # String format (length bytes per sequence + overhead)
    bytes_string = len(all_valid_sequences) * (length + 50)  # +50 for dict overhead
    print(f"  String format: {bytes_string:,} bytes = {bytes_string/(1024**2):.1f} MB")

    print(f"  Compression ratio: {bytes_string/bytes_encoded:.1f}x")

    # Projekcja na 20nt
    print("\n" + "=" * 70)
    print("PROJEKCJA NA 20nt")
    print("=" * 70)

    ratio_20_to_15 = (4**20) / (4**15)

    projected_time = total_time * ratio_20_to_15
    projected_valid = total_valid * ratio_20_to_15
    projected_memory = bytes_encoded * ratio_20_to_15

    print(f"\nPrzy założeniu podobnej dystrybucji:")
    print(f"  Total sequences: {4**20:,}")
    print(f"  Projected time (same hardware): {projected_time/3600:.1f} hours = {projected_time/86400:.1f} days")
    print(f"  Projected valid: ~{projected_valid:,.0f}")
    print(f"  Projected memory: {projected_memory/(1024**3):.1f} GB")

    print(f"\nZ {num_workers*4} cores (4x more):")
    print(f"  Projected time: {projected_time/4/3600:.1f} hours = {projected_time/4/86400:.1f} days")

    print("\n" + "=" * 70)

    return all_valid_sequences, total_fail_reasons

def main():
    """Main entry point."""

    # Check if running as test or full
    if len(sys.argv) > 1 and sys.argv[1] == '--test':
        print("RUNNING IN TEST MODE (12nt only)\n")
        length = 12  # 4^12 = 16M sequences (~20 seconds)
    else:
        length = 15  # 4^15 = 1B sequences (~5-6 hours on 8 cores)

    # Detect CPU count
    num_workers = mp.cpu_count()
    print(f"Detected {num_workers} CPU cores\n")

    # Run exhaustive search
    valid_sequences, fail_reasons = exhaustive_search_parallel(
        length=length,
        num_workers=num_workers
    )

    # Optional: save to file
    print(f"\nSaving valid sequences to file...")

    # Save in compressed format (just integers)
    with open(f'valid_sequences_{length}nt.txt', 'w') as f:
        f.write(f"# Valid {length}nt sequences (encoded as integers)\n")
        f.write(f"# Total: {len(valid_sequences)}\n")
        f.write(f"# Format: one integer per line\n")
        f.write(f"# Use decode_sequence(int, {length}) to get sequence\n\n")

        for encoded in valid_sequences[:100000]:  # Save first 100k for demo
            f.write(f"{encoded}\n")

    print(f"Saved first 100,000 sequences to valid_sequences_{length}nt.txt")

    print("\nDONE!")

if __name__ == '__main__':
    # Note: This will take 5-6 hours for 15nt on 8 cores
    # Use --test flag for quick demo
    main()
