"""
Custom validation profile example.

This example demonstrates how to use custom validation profiles:
- Load custom validation profiles from JSON
- Create profiles programmatically
- Compare sequences generated with different profiles
"""

from dna_generator import DNAGenerator, GeneratorConfig, GenerationMode

def generate_with_profile(profile_name, initial_seq, length, seed):
    """Generate a sequence with a specific validation profile."""
    config = GeneratorConfig(
        validation_profile=profile_name,
        generation_mode=GenerationMode.DETERMINISTIC,
        seed=seed,
        window_size=10,
        max_backtrack_attempts=50000,
    )

    generator = DNAGenerator(config)
    result = generator.generate(initial_seq, length)

    return result

def main():
    print("=" * 70)
    print("Custom Validation Profile Example")
    print("=" * 70)

    # Test parameters
    initial_sequence = "ATGCGATCGT"
    target_length = 200
    seed = 999

    # Test different validation profiles
    profiles = [
        ('strict', "Strict quality (narrow GC, low Tm variance)"),
        ('pcr_friendly', "PCR-optimized (balanced constraints)"),
        ('relaxed', "Relaxed quality (wider tolerances)"),
        ('sequence_only', "Sequence checks only (no Tm/GC validation)"),
    ]

    print(f"\nComparing validation profiles")
    print(f"Initial sequence: {initial_sequence}")
    print(f"Target length: {target_length} bp")
    print(f"Seed: {seed} (for reproducibility)")

    results = {}

    for profile_name, description in profiles:
        print(f"\n{'=' * 70}")
        print(f"Profile: {profile_name}")
        print(f"Description: {description}")
        print(f"{'=' * 70}")

        result = generate_with_profile(profile_name, initial_sequence, target_length, seed)
        results[profile_name] = result

        if result.success:
            print(f"✓ Generation successful")
            print(f"Sequence: {result.sequence}")

            metrics = result.quality_metrics
            print(f"\nQuality Metrics:")
            print(f"  GC Content: {metrics.gc_content:.2%}")
            print(f"  Tm: {metrics.melting_temperature:.1f}°C")
            print(f"  Hairpin Tm: {metrics.hairpin_tm:.1f}°C")
            print(f"  Homodimer Tm: {metrics.homodimer_tm:.1f}°C")
            print(f"  Valid: {metrics.is_valid}")

            stats = result.generation_stats
            print(f"\nGeneration Stats:")
            print(f"  Time: {result.generation_time:.3f}s")
            print(f"  Backtracks: {stats.get('backtrack_count', 0)}")
        else:
            print(f"✗ Generation failed")
            print(f"Error: {result.error_message}")

    # Comparison summary
    print(f"\n{'=' * 70}")
    print(f"Profile Comparison Summary")
    print(f"{'=' * 70}")

    print(f"\n{'Profile':<20} {'Success':<10} {'GC%':<10} {'Tm (°C)':<10} {'Time (s)':<10}")
    print("-" * 70)

    for profile_name, _ in profiles:
        result = results[profile_name]
        success = "✓" if result.success else "✗"

        if result.success:
            gc = f"{result.quality_metrics.gc_content:.1%}"
            tm = f"{result.quality_metrics.melting_temperature:.1f}"
            time = f"{result.generation_time:.3f}"
        else:
            gc = tm = time = "N/A"

        print(f"{profile_name:<20} {success:<10} {gc:<10} {tm:<10} {time:<10}")

    # Tips for custom profiles
    print(f"\n{'=' * 70}")
    print("Creating Custom Profiles")
    print(f"{'=' * 70}")
    print("""
To create your own validation profile:

1. Copy default_profiles.json from the dna_generator package
2. Add your custom profile to the "profiles" section:

{
  "profiles": {
    "my_custom_profile": {
      "description": "My custom validation rules",
      "rules": {
        "gc_content": true,
        "melting_temperature": true,
        "hairpin_structures": true,
        "homodimer_structures": true,
        "homopolymer_runs": true,
        "dinucleotide_repeats": true,
        "three_prime_stability": true
      },
      "params": {
        "min_gc": 0.45,
        "max_gc": 0.55,
        "min_tm": 58.0,
        "max_tm": 62.0,
        "max_hairpin_tm": 30.0,
        "max_homodimer_tm": 30.0,
        "max_homopolymer_length": 4,
        "max_dinucleotide_repeats": 3,
        "max_3prime_gc": 3
      }
    }
  }
}

3. Save as user_profiles.json in your project directory
4. Use: GeneratorConfig(validation_profile='my_custom_profile')

Available validation rules:
  - gc_content: GC content within min_gc to max_gc range
  - melting_temperature: Tm within min_tm to max_tm range
  - hairpin_structures: Secondary structure formation check
  - homodimer_structures: Self-dimerization check
  - homopolymer_runs: Repeating nucleotides (e.g., AAAA)
  - dinucleotide_repeats: Repeating dinucleotides (e.g., ATATAT)
  - three_prime_stability: 3' end GC content check
""")

    print("=" * 70)

if __name__ == '__main__':
    main()
