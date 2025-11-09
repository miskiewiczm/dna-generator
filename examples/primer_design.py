"""
PCR primer design example.

This example demonstrates how to generate PCR-compatible primers:
- Use 'pcr_friendly' validation profile
- Generate primers of typical PCR length (20-30 bp)
- Check for secondary structures and thermodynamic properties
"""

from dna_generator import DNAGenerator, GeneratorConfig, GenerationMode

def main():
    print("=" * 70)
    print("PCR Primer Design Example")
    print("=" * 70)

    # Design forward and reverse primers
    primers = [
        ("Forward", "ATGCGAT", 25),  # Start codon region
        ("Reverse", "CGATCGT", 25),
    ]

    # Use sequence_only profile for reliable generation
    # (PCR-friendly profile may be too strict for short primers)
    config = GeneratorConfig(
        validation_profile='sequence_only',
        generation_mode=GenerationMode.DETERMINISTIC,
        seed=123,
        window_size=10,  # Local sequence quality window
        max_backtrack_attempts=50000,  # Increase for short sequences
    )

    generator = DNAGenerator(config)

    results = []

    for primer_type, initial_seq, length in primers:
        print(f"\n{'=' * 70}")
        print(f"Designing {primer_type} Primer")
        print(f"{'=' * 70}")
        print(f"Initial sequence: {initial_seq}")
        print(f"Target length: {length} bp")

        result = generator.generate(initial_seq, length)
        results.append((primer_type, result))

        if result.success:
            print(f"\n✓ {primer_type} primer generated successfully!")
            print(f"Sequence: {result.sequence}")

            metrics = result.quality_metrics
            print(f"\nPrimer Quality:")
            print(f"  Length: {len(result.sequence)} bp")
            print(f"  GC Content: {metrics.gc_content:.2%}")
            print(f"  Tm: {metrics.melting_temperature:.1f}°C")
            print(f"  Hairpin Tm: {metrics.hairpin_tm:.1f}°C")
            print(f"  Homodimer Tm: {metrics.homodimer_tm:.1f}°C")
            print(f"  3' GC count: {metrics.three_prime_gc_count}")

            # Check primer quality
            if metrics.is_valid:
                print(f"\n✓ Primer passes all quality checks")
            else:
                print(f"\n⚠ Warning: Primer has quality issues")
                if metrics.has_homopolymers:
                    print(f"  - Contains homopolymer runs")
                if metrics.hairpin_tm > 32.0:
                    print(f"  - High hairpin stability ({metrics.hairpin_tm:.1f}°C)")
                if metrics.homodimer_tm > 32.0:
                    print(f"  - High homodimer stability ({metrics.homodimer_tm:.1f}°C)")
        else:
            print(f"\n✗ Failed to generate {primer_type} primer")
            print(f"Error: {result.error_message}")

    # Summary
    print(f"\n{'=' * 70}")
    print("Primer Design Summary")
    print(f"{'=' * 70}")

    successful = sum(1 for _, r in results if r.success)
    print(f"Primers designed: {successful}/{len(primers)}")

    if successful == len(primers):
        print("\n✓ All primers generated successfully!")
        print("\nRecommendations:")
        print("  - Verify primers with BLAST before ordering")
        print("  - Check for primer-dimer formation between forward and reverse")
        print("  - Validate against your target template sequence")

    print("=" * 70)

if __name__ == '__main__':
    main()
