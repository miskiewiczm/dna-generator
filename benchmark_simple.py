#!/usr/bin/env python3
"""
Uproszczony benchmark - oszacowanie kosztów obliczeniowych.
Symuluje walidację bez pełnego stacku.
"""

import time
import random

def simple_gc_content(seq):
    """Oblicz zawartość GC."""
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq)

def simple_tm_calculation(seq):
    """Uproszczone obliczenie Tm (nearest neighbor approximation)."""
    # Symulacja - w rzeczywistości to jest kosztowne
    gc = simple_gc_content(seq)
    # Wallace rule (bardzo uproszczone)
    tm = 4 * (seq.count('G') + seq.count('C')) + 2 * (seq.count('A') + seq.count('T'))
    # Dodatkowe operacje stringowe (symulacja primer3)
    for i in range(len(seq)-1):
        dinuc = seq[i:i+2]  # Nearest neighbor
    return tm

def check_homopolymers(seq, max_len=4):
    """Sprawdź homopolimery."""
    for base in 'ATGC':
        if base * (max_len + 1) in seq:
            return False
    return True

def check_dinucleotide_repeats(seq, max_repeats=3):
    """Sprawdź powtórzenia dinukleotydowe."""
    for i in range(4):
        for j in range(4):
            bases = 'ATGC'
            dinuc = bases[i] + bases[j]
            repeat = dinuc * (max_repeats + 1)
            if repeat in seq:
                return False
    return True

def check_three_prime_stability(seq):
    """Sprawdź stabilność końca 3'."""
    last_5 = seq[-5:]
    gc_count = last_5.count('G') + last_5.count('C')
    return gc_count <= 3

def simple_hairpin_tm(seq):
    """Symulowana kalkulacja hairpin Tm (kosztowna)."""
    # W rzeczywistości to wywołanie primer3 - bardzo kosztowne
    # Symulujemy iterację po możliwych strukturach
    n = len(seq)
    score = 0
    for i in range(n-3):
        for j in range(i+3, n):
            # Sprawdź komplementarność
            if seq[i] == 'A' and seq[j] == 'T':
                score += 2
            elif seq[i] == 'T' and seq[j] == 'A':
                score += 2
            elif seq[i] == 'G' and seq[j] == 'C':
                score += 4
            elif seq[i] == 'C' and seq[j] == 'G':
                score += 4
    return score * 0.5  # Simplified Tm

def simple_homodimer_tm(seq):
    """Symulowana kalkulacja homodimer Tm."""
    # Podobnie kosztowne jak hairpin
    n = len(seq)
    max_score = 0
    for offset in range(1, n):
        score = 0
        for i in range(min(n-offset, n)):
            if i+offset >= n:
                break
            if seq[i] == 'A' and seq[i+offset] == 'T':
                score += 2
            elif seq[i] == 'T' and seq[i+offset] == 'A':
                score += 2
            elif seq[i] == 'G' and seq[i+offset] == 'C':
                score += 4
            elif seq[i] == 'C' and seq[i+offset] == 'G':
                score += 4
        max_score = max(max_score, score)
    return max_score * 0.3

def validate_sequence_simple(seq, min_gc=0.4, max_gc=0.6, min_tm=54, max_tm=66):
    """Uproszczona walidacja - symuluje pełny pipeline."""

    # 1. GC content (fast)
    gc = simple_gc_content(seq)
    if gc < min_gc or gc > max_gc:
        return False

    # 2. Tm calculation (medium cost)
    tm = simple_tm_calculation(seq)
    if tm < min_tm or tm > max_tm:
        return False

    # 3. Homopolymers (fast)
    if not check_homopolymers(seq):
        return False

    # 4. Dinucleotide repeats (medium)
    if not check_dinucleotide_repeats(seq):
        return False

    # 5. 3' stability (fast)
    if not check_three_prime_stability(seq):
        return False

    # 6. Hairpin Tm (EXPENSIVE)
    hairpin = simple_hairpin_tm(seq)
    if hairpin > 32:
        return False

    # 7. Homodimer Tm (EXPENSIVE)
    homodimer = simple_homodimer_tm(seq)
    if homodimer > 32:
        return False

    return True

def generate_random_sequence(length=20):
    """Generuj losową sekwencję."""
    bases = 'ATGC'
    return ''.join(random.choice(bases) for _ in range(length))

def benchmark():
    """Benchmark walidacji."""

    print("=" * 70)
    print("Uproszczony benchmark walidacji 20nt sekwencji")
    print("=" * 70)

    # Generate test sequences
    test_sequences = [generate_random_sequence(20) for _ in range(100)]

    # Warm-up
    for seq in test_sequences[:10]:
        validate_sequence_simple(seq)

    # Benchmark
    print("\nBenchmarking 1000 walidacji...")
    start = time.perf_counter()

    valid_count = 0
    for _ in range(1000):
        seq = random.choice(test_sequences)
        if validate_sequence_simple(seq):
            valid_count += 1

    end = time.perf_counter()

    total_time = end - start
    avg_time = total_time / 1000

    print(f"\nWyniki:")
    print(f"  Total time: {total_time:.3f} seconds")
    print(f"  Avg time per validation: {avg_time*1e6:.2f} µs")
    print(f"  Throughput: {1000/total_time:,.0f} validations/sec")
    print(f"  Valid sequences: {valid_count}/1000 ({valid_count/10:.1f}%)")

    # Projekcje
    print("\n" + "=" * 70)
    print("PROJEKCJE dla 4^20 = 1,099,511,627,776 sekwencji")
    print("=" * 70)

    total_sequences = 4**20

    print(f"\nPodstawowe liczby:")
    print(f"  Total sequences: {total_sequences:,}")
    print(f"  Avg time per seq: {avg_time*1e6:.2f} µs = {avg_time*1e3:.2f} ms")

    # Single thread
    single_thread_seconds = total_sequences * avg_time
    single_thread_hours = single_thread_seconds / 3600
    single_thread_days = single_thread_seconds / 86400
    single_thread_years = single_thread_days / 365

    print(f"\n[1] SINGLE THREAD:")
    print(f"  Total time: {single_thread_seconds:,.0f} seconds")
    print(f"            = {single_thread_hours:,.0f} hours")
    print(f"            = {single_thread_days:,.0f} days")
    print(f"            = {single_thread_years:,.1f} years")
    print(f"  Throughput: {1/avg_time:,.0f} seq/sec")

    # Parallel scenarios
    print(f"\n[2] PARALLEL PROCESSING:")

    scenarios = [
        ("Desktop (8 cores)", 8),
        ("Workstation (16 cores)", 16),
        ("Server (32 cores)", 32),
        ("Server (64 cores)", 64),
        ("Small cluster (128 cores)", 128),
        ("Medium cluster (256 cores)", 256),
        ("Large cluster (512 cores)", 512),
        ("Large cluster (1024 cores)", 1024),
        ("Supercomputer (10,000 cores)", 10000),
    ]

    for desc, cores in scenarios:
        parallel_seconds = single_thread_seconds / cores
        parallel_days = parallel_seconds / 86400
        parallel_years = parallel_days / 365
        throughput = cores / avg_time

        print(f"\n  {desc}:")
        print(f"    Cores: {cores}")
        print(f"    Time: {parallel_days:,.0f} days = {parallel_years:.2f} years")
        print(f"    Throughput: {throughput:,.0f} seq/sec")

    # Memory requirements
    print(f"\n[3] PAMIĘĆ:")

    # Minimal: sequence string (20 bytes) + is_valid (1 byte) + dict overhead (~40 bytes)
    bytes_per_minimal = 70
    total_bytes_minimal = total_sequences * bytes_per_minimal
    total_gb_minimal = total_bytes_minimal / (1024**3)
    total_tb_minimal = total_bytes_minimal / (1024**4)
    total_pb_minimal = total_bytes_minimal / (1024**5)

    print(f"\n  Minimalna struktura (sequence + is_valid):")
    print(f"    Per entry: ~{bytes_per_minimal} bytes")
    print(f"    Total: {total_gb_minimal:,.0f} GB")
    print(f"         = {total_tb_minimal:,.0f} TB")
    print(f"         = {total_pb_minimal:.2f} PB (petabytes)")

    # With full metrics
    bytes_per_full = 150
    total_bytes_full = total_sequences * bytes_per_full
    total_pb_full = total_bytes_full / (1024**5)

    print(f"\n  Pełna struktura (+ metryki jakości):")
    print(f"    Per entry: ~{bytes_per_full} bytes")
    print(f"    Total: {total_pb_full:.2f} PB (petabytes)")

    # Pruning analysis
    print(f"\n[4] ANALIZA Z PRUNINGIEM:")

    # Assume 80% fail fast checks (GC content, homopolymers, etc)
    fast_fail_rate = 0.80
    remaining_after_fast = total_sequences * (1 - fast_fail_rate)

    # Of remaining, 70% fail expensive checks (hairpin, homodimer)
    expensive_fail_rate = 0.70
    final_valid = remaining_after_fast * (1 - expensive_fail_rate)

    # Avg time with pruning (80% fail fast at 1µs, 20% full check at avg_time)
    fast_check_time = 1e-6  # 1 microsecond
    avg_time_with_pruning = (fast_fail_rate * fast_check_time +
                             (1 - fast_fail_rate) * avg_time)

    pruned_total_seconds = total_sequences * avg_time_with_pruning
    pruned_years = pruned_total_seconds / 86400 / 365

    print(f"\n  Założenia:")
    print(f"    - {fast_fail_rate*100:.0f}% fail fast checks (~1µs)")
    print(f"    - {(1-fast_fail_rate)*100:.0f}% reach expensive checks")
    print(f"    - Final valid rate: ~{final_valid/total_sequences*100:.2f}%")

    print(f"\n  Z pruningiem:")
    print(f"    Avg time per seq: {avg_time_with_pruning*1e6:.2f} µs")
    print(f"    Single thread: {pruned_years:.1f} years")
    print(f"    1024 cores: {pruned_years/1024*365:.0f} days")

    print(f"\n  Final dictionary size:")
    print(f"    Valid sequences: ~{final_valid:,.0f}")
    print(f"    Memory: ~{final_valid * bytes_per_full / (1024**3):.0f} GB")

    print("\n" + "=" * 70)

if __name__ == '__main__':
    benchmark()
