"""
Basic DNA sequence generation example.

This example demonstrates the most basic usage of the DNA generator:
- Generate a single DNA sequence
- Use default (relaxed) validation profile
- Display sequence and quality metrics
"""

from dna_generator import DNAGenerator, GeneratorConfig, GenerationMode

def main():
    print("=" * 70)
    print("DNA Sequence Generator - Basic Example")
    print("=" * 70)

    # Initial sequence (e.g., start of a primer or known sequence)
    initial_sequence = "GAGACGAGCATGAGACATAC"  # Real validated primer
    target_length = 200  # Generate a 200bp sequence

    print(f"\nInitial sequence: {initial_sequence}")
    print(f"Target length: {target_length} bp")

    # Create generator with default configuration
    # By default, uses 'relaxed' profile and deterministic mode
    config = GeneratorConfig(
        validation_profile='sequence_only',  # Start with basic validation
        generation_mode=GenerationMode.DETERMINISTIC,
        seed=42,  # For reproducibility
        max_backtrack_attempts=50000  # Increase attempts for reliability
    )

    generator = DNAGenerator(config)

    # Generate sequence
    print("\nGenerating sequence...")
    result = generator.generate(initial_sequence, target_length)

    # Display results
    print("\n" + "=" * 70)
    if result.success:
        print("✓ Generation successful!")
        print(f"\nGenerated sequence ({len(result.sequence)} bp):")
        print(result.sequence)

        # Display quality metrics
        print("\nQuality Metrics:")
        metrics = result.quality_metrics
        print(f"  GC Content: {metrics.gc_content:.2%}")
        print(f"  Melting Temperature: {metrics.melting_temperature:.1f}°C")
        print(f"  Valid: {metrics.is_valid}")
        print(f"  Has Homopolymers: {metrics.has_homopolymers}")
        if metrics.longest_homopolymer:
            base, length = metrics.longest_homopolymer
            print(f"  Longest Homopolymer: {base} × {length}")

        # Display generation statistics
        print("\nGeneration Statistics:")
        stats = result.generation_stats
        print(f"  Generation time: {result.generation_time:.3f}s")
        print(f"  Backtrack count: {stats.get('backtrack_count', 0)}")

    else:
        print("✗ Generation failed!")
        print(f"Error: {result.error_message}")

    print("=" * 70)

if __name__ == '__main__':
    main()
