"""
DNA sequence library generation example.

This example demonstrates how to generate a library of multiple sequences:
- Generate multiple diverse sequences
- Use random mode for sequence variety
- Analyze diversity across the library
- Export sequences for downstream use
"""

from dna_generator import DNAGenerator, GeneratorConfig, GenerationMode

def analyze_library_diversity(sequences):
    """Calculate diversity metrics for a sequence library."""
    unique_sequences = set(sequences)

    # Calculate nucleotide diversity at each position
    if not sequences or not sequences[0]:
        return None

    length = len(sequences[0])
    position_diversity = []

    for pos in range(length):
        nucleotides = [seq[pos] for seq in sequences if len(seq) > pos]
        unique_nt = len(set(nucleotides))
        diversity = unique_nt / 4.0  # Normalized by max possible (A,T,G,C)
        position_diversity.append(diversity)

    avg_diversity = sum(position_diversity) / len(position_diversity) if position_diversity else 0

    return {
        'total_sequences': len(sequences),
        'unique_sequences': len(unique_sequences),
        'uniqueness_ratio': len(unique_sequences) / len(sequences),
        'average_positional_diversity': avg_diversity,
        'min_diversity': min(position_diversity) if position_diversity else 0,
        'max_diversity': max(position_diversity) if position_diversity else 0,
    }

def main():
    print("=" * 70)
    print("DNA Sequence Library Generation")
    print("=" * 70)

    # Library parameters
    initial_sequence = "CCTGTCATCACGCTAGTAAC"
    target_length = 150
    library_size = 10

    print(f"\nLibrary Parameters:")
    print(f"  Initial sequence: {initial_sequence}")
    print(f"  Target length: {target_length} bp")
    print(f"  Library size: {library_size} sequences")

    # Use random mode for diversity
    config = GeneratorConfig(
        validation_profile='sequence_only',
        generation_mode=GenerationMode.RANDOM,  # Random for diversity
        window_size=10,
        max_backtrack_attempts=50000,
    )

    generator = DNAGenerator(config)

    # Generate library
    print(f"\nGenerating {library_size} sequences...")
    results = generator.generate_multiple(
        initial_sequence,
        target_length,
        count=library_size
    )

    # Collect successful sequences
    sequences = [r.sequence for r in results if r.success]

    print(f"\n{'=' * 70}")
    print(f"Generation Results")
    print(f"{'=' * 70}")
    print(f"Successful: {len(sequences)}/{library_size}")

    if sequences:
        # Display first few sequences
        print(f"\nFirst 3 sequences:")
        for i, seq in enumerate(sequences[:3], 1):
            print(f"  {i}. {seq}")
        if len(sequences) > 3:
            print(f"  ... and {len(sequences) - 3} more")

        # Analyze diversity
        diversity = analyze_library_diversity(sequences)

        print(f"\n{'=' * 70}")
        print(f"Library Diversity Analysis")
        print(f"{'=' * 70}")
        print(f"Total sequences: {diversity['total_sequences']}")
        print(f"Unique sequences: {diversity['unique_sequences']}")
        print(f"Uniqueness ratio: {diversity['uniqueness_ratio']:.2%}")
        print(f"Average positional diversity: {diversity['average_positional_diversity']:.2%}")
        print(f"Diversity range: {diversity['min_diversity']:.2%} - {diversity['max_diversity']:.2%}")

        # Quality summary
        print(f"\n{'=' * 70}")
        print(f"Quality Summary")
        print(f"{'=' * 70}")

        valid_count = sum(1 for r in results if r.success and r.quality_metrics.is_valid)
        avg_gc = sum(r.quality_metrics.gc_content for r in results if r.success) / len(sequences)
        avg_tm = sum(r.quality_metrics.melting_temperature for r in results if r.success) / len(sequences)

        print(f"Valid sequences: {valid_count}/{len(sequences)}")
        print(f"Average GC content: {avg_gc:.2%}")
        print(f"Average Tm: {avg_tm:.1f}°C")

        # Export option
        print(f"\n{'=' * 70}")
        print(f"Export")
        print(f"{'=' * 70}")
        print("To export sequences to a file, use:")
        print("  with open('library.fasta', 'w') as f:")
        print("      for i, seq in enumerate(sequences, 1):")
        print("          f.write(f'>sequence_{i}\\n{seq}\\n')")

    else:
        print("\n✗ No sequences were successfully generated")

    print("=" * 70)

if __name__ == '__main__':
    main()
