#!/usr/bin/env python3
"""
Benchmark walidacji sekwencji - oszacowanie czasu dla single sequence.
"""

import time
import sys
sys.path.insert(0, '/home/user/dna-generator')

from dna_generator import DNAValidator, GeneratorConfig, DeterministicRandom
from dna_commons import ValidationRules, ThermodynamicParams

def benchmark_validation():
    """Benchmark pojedynczej walidacji."""

    # Config jak w projekcie
    config = GeneratorConfig(validation_profile='pcr_friendly')

    # Create validator
    rules = ValidationRules(
        gc_content=config.validation_rules.get('gc_content', True),
        melting_temperature=config.validation_rules.get('melting_temperature', True),
        homopolymer_runs=config.validation_rules.get('homopolymer_runs', True),
        dinucleotide_repeats=config.validation_rules.get('dinucleotide_repeats', True),
        three_prime_stability=config.validation_rules.get('three_prime_stability', True),
        hairpin_structures=config.validation_rules.get('hairpin_structures', True),
        homodimer_structures=config.validation_rules.get('homodimer_structures', True),
        min_gc=config.min_gc,
        max_gc=config.max_gc,
        min_tm=config.min_tm,
        max_tm=config.max_tm,
        max_hairpin_tm=config.max_hairpin_tm,
        max_homodimer_tm=config.max_homodimer_tm,
        max_homopolymer_length=config.max_homopolymer_length,
        max_dinucleotide_repeats=config.max_dinucleotide_repeats
    )

    thermoparams = ThermodynamicParams.load_default()
    validator = DNAValidator(rules, thermoparams)

    # Testowe sekwencje 20nt
    test_sequences = [
        "ATGCATGCATGCATGCATGC",
        "GCTAGCTAGCTAGCTAGCTA",
        "AATTCCGGAATTCCGGAATT",
        "CGCGATATCGCGATCGCGAT",
        "TACGTACGTACGTACGTACG",
    ]

    print("=" * 70)
    print("Benchmark walidacji sekwencji 20nt")
    print("=" * 70)

    # Warm-up
    for seq in test_sequences[:2]:
        validator.validate_window(seq)

    # Benchmark validate_window (używane w backtracking)
    print("\n[1] validate_window() - używane w backtracking loop")
    times_window = []
    for seq in test_sequences:
        start = time.perf_counter()
        for _ in range(1000):
            result = validator.validate_window(seq)
        end = time.perf_counter()
        avg_time = (end - start) / 1000
        times_window.append(avg_time)
        print(f"  {seq}: {avg_time*1e6:.2f} µs per validation (valid={result})")

    avg_window = sum(times_window) / len(times_window)
    print(f"\n  Average: {avg_window*1e6:.2f} µs per validation")

    # Benchmark validate_sequence (pełna walidacja)
    print("\n[2] validate_sequence() - pełna walidacja z metrykami")
    times_full = []
    for seq in test_sequences:
        start = time.perf_counter()
        for _ in range(1000):
            metrics = validator.validate_sequence(seq)
        end = time.perf_counter()
        avg_time = (end - start) / 1000
        times_full.append(avg_time)
        print(f"  {seq}: {avg_time*1e6:.2f} µs per validation")

    avg_full = sum(times_full) / len(times_full)
    print(f"\n  Average: {avg_full*1e6:.2f} µs per validation")

    # Projekcje
    print("\n" + "=" * 70)
    print("Projekcje czasowe dla 4^20 = 1,099,511,627,776 sekwencji")
    print("=" * 70)

    total_sequences = 4**20

    print(f"\nUsługa validate_window() (szybsza, używana w backtracking):")
    print(f"  Avg time per seq: {avg_window*1e6:.2f} µs")
    single_thread_seconds_window = total_sequences * avg_window
    print(f"  Single thread: {single_thread_seconds_window:,.0f} seconds")
    print(f"               = {single_thread_seconds_window/3600:,.0f} hours")
    print(f"               = {single_thread_seconds_window/86400:,.0f} days")
    print(f"               = {single_thread_seconds_window/86400/365:.1f} years")

    print(f"\nUsługa validate_sequence() (pełna walidacja):")
    print(f"  Avg time per seq: {avg_full*1e6:.2f} µs")
    single_thread_seconds_full = total_sequences * avg_full
    print(f"  Single thread: {single_thread_seconds_full:,.0f} seconds")
    print(f"               = {single_thread_seconds_full/86400:,.0f} days")
    print(f"               = {single_thread_seconds_full/86400/365:.1f} years")

    # Parallel scenarios
    print("\n" + "=" * 70)
    print("Scenariusze równoległe (validate_window)")
    print("=" * 70)

    for cores in [8, 16, 32, 64, 128, 256, 1024]:
        parallel_seconds = single_thread_seconds_window / cores
        print(f"\n  {cores:4d} cores:")
        print(f"    Time: {parallel_seconds/86400:,.0f} days = {parallel_seconds/86400/365:.1f} years")
        print(f"    Throughput: {cores/avg_window:,.0f} seq/sec")

    # Pamięć
    print("\n" + "=" * 70)
    print("Wymagania pamięciowe")
    print("=" * 70)

    # Minimalna struktura: (sequence_str, is_valid)
    # sequence: 20 bytes (string) + overhead ~50 bytes
    # is_valid: 1 byte (bool)
    # Dict overhead: ~40 bytes per entry
    # Total per entry: ~110 bytes (conservative)

    bytes_per_entry = 110
    total_bytes = total_sequences * bytes_per_entry
    total_gb = total_bytes / (1024**3)
    total_tb = total_bytes / (1024**4)
    total_pb = total_bytes / (1024**5)

    print(f"\n  Minimalna struktura (sequence + is_valid):")
    print(f"    ~{bytes_per_entry} bytes per entry")
    print(f"    Total: {total_gb:,.0f} GB = {total_tb:,.0f} TB = {total_pb:.1f} PB")

    # Z pełnymi metrykami
    bytes_per_full = 200  # sequence + QualityMetrics
    total_bytes_full = total_sequences * bytes_per_full
    total_pb_full = total_bytes_full / (1024**5)

    print(f"\n  Pełna struktura (sequence + QualityMetrics):")
    print(f"    ~{bytes_per_full} bytes per entry")
    print(f"    Total: {total_pb_full:.1f} PB")

    print("\n" + "=" * 70)

    return avg_window, avg_full

if __name__ == '__main__':
    benchmark_validation()
