#!/usr/bin/env python3
"""
Command-line entry point for the DNA sequence generator.

Examples:
    python -m dna_generator --initial ATGCATGC --length 100 --mode deterministic
    python -m dna_generator --initial ATGCATGC --length 100 --mode random --count 3
"""

import argparse
import logging
import sys
from pathlib import Path

from .generator import DNAGenerator
from .config import GeneratorConfig, GenerationMode
from .exceptions import GeneratorError
from dna_commons import PRIMER3_AVAILABLE


def create_parser() -> argparse.ArgumentParser:
    """Create command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="DNA sequence generator with biochemical quality control",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:

  # Deterministic mode (default)
  python -m dna_generator --initial ATGCATGC --length 100

  # Random mode with multiple sequences
  python -m dna_generator --initial ATGCATGC --length 100 --mode random --count 3
  
  # Custom quality parameters
  python -m dna_generator --initial ATGCATGC --length 100 --min-gc 0.4 --max-gc 0.6
  
  # With a specific seed
  python -m dna_generator --initial ATGCATGC --length 100 --seed 12345
        """
    )
    
    # Podstawowe parametry
    parser.add_argument(
        "--initial", "-i",
        required=True,
        help="Initial DNA sequence"
    )
    
    parser.add_argument(
        "--length", "-l",
        type=int,
        required=True,
        help="Target sequence length"
    )
    
    parser.add_argument(
        "--mode", "-m",
        choices=["deterministic", "random"],
        default="deterministic",
        help="Generation mode (default: deterministic)"
    )
    
    parser.add_argument(
        "--count", "-c",
        type=int,
        default=1,
        help="Number of sequences to generate (default: 1)"
    )
    
    parser.add_argument(
        "--seed", "-s",
        type=int,
        help="Seed for the generator (deterministic mode only)"
    )
    
    # Parametry jakości
    quality_group = parser.add_argument_group("Sequence quality parameters")
    
    quality_group.add_argument(
        "--min-gc",
        type=float,
        default=0.45,
        help="Minimum GC content (0.0-1.0, default: 0.45)"
    )
    
    quality_group.add_argument(
        "--max-gc",
        type=float,
        default=0.55,
        help="Maximum GC content (0.0-1.0, default: 0.55)"
    )
    
    quality_group.add_argument(
        "--min-tm",
        type=float,
        default=55.0,
        help="Minimum melting temperature (°C, default: 55.0)"
    )
    
    quality_group.add_argument(
        "--max-tm",
        type=float,
        default=65.0,
        help="Maximum melting temperature (°C, default: 65.0)"
    )
    
    quality_group.add_argument(
        "--max-hairpin",
        type=float,
        default=30.0,
        help="Maximum hairpin Tm (°C, default: 30.0)"
    )
    
    quality_group.add_argument(
        "--max-homodimer",
        type=float,
        default=30.0,
        help="Maximum homodimer Tm (°C, default: 30.0)"
    )
    
    quality_group.add_argument(
        "--window-size",
        type=int,
        default=20,
        help="Analysis window size (default: 20)"
    )
    
    # Parametry walidacji
    validation_group = parser.add_argument_group("Validation control")
    validation_group.add_argument(
        "--profile",
        choices=['pcr_friendly', 'strict', 'relaxed', 'sequence_only', 'none', 'user'],
        help="Validation profile (pcr_friendly/strict/relaxed/sequence_only/none/user). With 'user', load JSON via --profile-file."
    )
    validation_group.add_argument(
        "--profile-file",
        type=str,
        help="Path to JSON profile when using --profile user (keys: rules, params)"
    )
    
    validation_group.add_argument(
        "--no-gc-check",
        action="store_true",
        help="Disable GC content check"
    )
    
    validation_group.add_argument(
        "--no-tm-check",
        action="store_true",
        help="Disable melting temperature check"
    )
    
    validation_group.add_argument(
        "--no-hairpin-check",
        action="store_true",
        help="Disable hairpin structures check"
    )
    
    validation_group.add_argument(
        "--no-homodimer-check",
        action="store_true",
        help="Disable homodimer structures check"
    )
    
    validation_group.add_argument(
        "--no-homopolymer-check",
        action="store_true",
        help="Disable homopolymer check"
    )
    
    validation_group.add_argument(
        "--no-dinucleotide-check",
        action="store_true",
        help="Disable dinucleotide repeats check"
    )
    
    validation_group.add_argument(
        "--no-3prime-check",
        action="store_true",
        help="Disable 3' end stability check"
    )
    
    # Parametry algorytmu
    algo_group = parser.add_argument_group("Algorithm parameters")
    
    algo_group.add_argument(
        "--max-attempts",
        type=int,
        default=10000,
        help="Maximum number of backtracking attempts (default: 10000)"
    )
    
    # Heuristics control (mutually exclusive)
    heuristics_group = algo_group.add_mutually_exclusive_group()
    heuristics_group.add_argument(
        "--heuristics",
        dest="heuristics_override",
        action="store_const",
        const=True,
        help="Force enable backtracking heuristics (overrides profile auto-disable)"
    )
    heuristics_group.add_argument(
        "--no-heuristics",
        dest="heuristics_override",
        action="store_const",
        const=False,
        help="Force disable backtracking heuristics (overrides profile defaults)"
    )
    
    # Parametry wyjścia
    output_group = parser.add_argument_group("Output parameters")
    
    output_group.add_argument(
        "--output", "-o",
        type=str,
        help="Output file (default: stdout)"
    )
    
    output_group.add_argument(
        "--format",
        choices=["text", "fasta", "json"],
        default="text",
        help="Output format (default: text)"
    )
    
    output_group.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    
    output_group.add_argument(
        "--quiet", "-q",
        action="store_true",
        help="Disable progress logging"
    )

    output_group.add_argument(
        "--sequences-only",
        action="store_true",
        help="Output sequences only (no comments/statistics)"
    )

    output_group.add_argument(
        "--json-raw-metrics",
        action="store_true",
        help="In JSON output, return raw numeric metrics"
    )

    output_group.add_argument(
        "--show-primer3-status",
        action="store_true",
        help="Show primer3 availability status"
    )

    output_group.add_argument(
        "--csv-file",
        type=str,
        help="Export sliding window analysis to CSV file (20nt windows, step=1)"
    )

    return parser


def create_config_from_args(args) -> GeneratorConfig:
    """
    Create configuration from CLI arguments with clear priority hierarchy.

    Priority (lowest to highest):
    1. Defaults (from dataclass)
    2. Profile (strict/relaxed/none/etc.)
    3. Explicit CLI arguments (user's explicit choice)

    This ensures profiles set smart defaults, but users can always override.
    """
    mode = GenerationMode.DETERMINISTIC if args.mode == "deterministic" else GenerationMode.RANDOM

    # Prepare validation rules only when no profile is chosen
    validation_rules = None
    if not args.profile:  # indywidualne flagi wyłączeń
        validation_rules = {
            'gc_content': not getattr(args, 'no_gc_check', False),
            'melting_temperature': not getattr(args, 'no_tm_check', False),
            'hairpin_structures': not getattr(args, 'no_hairpin_check', False),
            'homodimer_structures': not getattr(args, 'no_homodimer_check', False),
            'homopolymer_runs': not getattr(args, 'no_homopolymer_check', False),
            'dinucleotide_repeats': not getattr(args, 'no_dinucleotide_check', False),
            'three_prime_stability': not getattr(args, 'no_3prime_check', False)
        }

    base_kwargs = dict(
        generation_mode=mode,
        seed=args.seed,
        window_size=args.window_size,
        max_backtrack_attempts=args.max_attempts,
        enable_progress_logging=not args.quiet,
        log_level="DEBUG" if args.verbose else "INFO",
    )

    # If a built-in profile is provided, use it and ignore CLI thresholds
    if args.profile and args.profile != 'user':
        config = GeneratorConfig(validation_profile=args.profile, **base_kwargs)
    else:
        # No profile or custom 'user' profile handled here
        config = GeneratorConfig(
            min_gc=args.min_gc,
            max_gc=args.max_gc,
            min_tm=args.min_tm,
            max_tm=args.max_tm,
            max_hairpin_tm=args.max_hairpin,
            max_homodimer_tm=args.max_homodimer,
            **base_kwargs
        )

        # Apply individual flags if no profile at all
        if validation_rules and not args.profile:
            config.validation_rules = validation_rules

        # Apply user profile from file: overrides thresholds/rules if provided
        if args.profile == 'user':
            if not args.profile_file:
                raise ValueError("--profile user requires --profile-file <path to JSON>.")
            import json
            with open(args.profile_file, 'r') as f:
                data = json.load(f)
            rules = data.get('rules')
            if isinstance(rules, dict):
                config.validation_rules = {**config.validation_rules, **rules}
            params = data.get('params', {})
            if isinstance(params, dict):
                for k, v in params.items():
                    if hasattr(config, k):
                        setattr(config, k, v)

    # CLI OVERRIDES: Only apply if user explicitly specified the flag
    # This respects profile auto-settings but allows user override

    # Heuristics: Apply explicit override if user specified --heuristics or --no-heuristics
    # If neither flag is specified (heuristics_override=None), respect profile defaults
    if args.heuristics_override is not None:
        config.enable_backtrack_heuristics = args.heuristics_override

    # Log active rules and thresholds when verbose
    if args.verbose:
        import json as _json
        print("Active rules:", _json.dumps(config.validation_rules, indent=2, ensure_ascii=False))
        params_view = {
            'min_gc': config.min_gc,
            'max_gc': config.max_gc,
            'min_tm': config.min_tm,
            'max_tm': config.max_tm,
            'max_hairpin_tm': config.max_hairpin_tm,
            'max_homodimer_tm': config.max_homodimer_tm,
            'max_homopolymer_length': config.max_homopolymer_length,
            'max_dinucleotide_repeats': config.max_dinucleotide_repeats,
            'max_3prime_gc': config.max_3prime_gc,
            'window_size': config.window_size,
        }
        print("Thresholds:", _json.dumps(params_view, indent=2, ensure_ascii=False))
    
    return config


def export_windows_to_csv(sequence: str, window_size: int, validator, csv_path: str):
    """
    Export sliding window analysis to CSV file.

    Args:
        sequence: Full DNA sequence to analyze
        window_size: Size of sliding window (typically 20)
        validator: DNAValidator instance for quality checks
        csv_path: Path to output CSV file
    """
    import csv

    rows = []

    # Analyze each window (sliding with step=1)
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]

        # Get comprehensive metrics for this window
        metrics = validator.validate_sequence(window)

        row = {
            'window_start': i,
            'window_end': i + window_size,
            'sequence': window,
            'gc_content': f"{metrics.gc_content:.4f}",
            'melting_temperature': f"{metrics.melting_temperature:.2f}",
            'hairpin_tm': f"{metrics.hairpin_tm:.2f}",
            'homodimer_tm': f"{metrics.homodimer_tm:.2f}",
            'has_homopolymers': metrics.has_homopolymers,
            'has_dinucleotide_repeats': metrics.has_dinucleotide_repeats,
            'three_prime_gc_count': metrics.three_prime_gc_count,
            'is_valid': metrics.is_valid,
            'longest_homopolymer': str(metrics.longest_homopolymer) if metrics.longest_homopolymer else '',
            'max_dinucleotide_repeat': str(metrics.max_dinucleotide_repeat) if metrics.max_dinucleotide_repeat else '',
        }
        rows.append(row)

    # Write to CSV
    fieldnames = [
        'window_start', 'window_end', 'sequence',
        'gc_content', 'melting_temperature', 'hairpin_tm', 'homodimer_tm',
        'has_homopolymers', 'has_dinucleotide_repeats', 'three_prime_gc_count',
        'is_valid', 'longest_homopolymer', 'max_dinucleotide_repeat'
    ]

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"CSV analysis exported to: {csv_path}")
    print(f"  Total windows analyzed: {len(rows)}")
    print(f"  Window size: {window_size}bp")
    print(f"  Step size: 1bp (sliding window)")


def format_output(results, format_type, initial_sequence, sequences_only=False, json_raw_metrics=False):
    """Format results for requested output type.
    
    Args:
        results: List of generation results
        format_type: Output format type
        initial_sequence: Initial sequence
        sequences_only: If True, outputs sequences only (no comments/statistics)
    """
    if format_type == "fasta":
        output = []
        for i, result in enumerate(results, 1):
            if result.success:
                if sequences_only:
                    header = f">sequence_{i}"
                else:
                    header = f">sequence_{i}|length={result.actual_length}|gc={result.quality_metrics.gc_content:.2%}"
                output.append(header)
                output.append(result.sequence)
        return "\n".join(output)
    
    elif format_type == "json":
        import json
        if sequences_only:
            sequences = [result.sequence for result in results if result.success]
            return json.dumps(sequences, indent=2)
        else:
            return json.dumps([result.to_dict(raw=json_raw_metrics) for result in results], indent=2)
    
    else:  # text format
        if sequences_only:
            # Sequences only, one per line
            sequences = [result.sequence for result in results if result.success]
            return "\n".join(sequences)
        else:
            # Full format with comments and statistics
            output = []
            for i, result in enumerate(results, 1):
                output.append(f"=== SEQUENCE {i} ===")
                if result.success:
                    output.append(f"Success: YES")
                    output.append(f"Sequence: {result.sequence}")
                    output.append(f"Length: {result.actual_length}")
                    output.append(f"Generation time: {result.generation_time:.2f}s")
                    output.append(f"GC content: {result.quality_metrics.gc_content:.2%}")
                    output.append(f"Melting temperature: {result.quality_metrics.melting_temperature:.1f}°C")
                    output.append(f"Hairpin Tm: {result.quality_metrics.hairpin_tm:.1f}°C")
                    output.append(f"Homodimer Tm: {result.quality_metrics.homodimer_tm:.1f}°C")
                    # Window compliance / global compliance info
                    output.append(f"Meets window criteria: YES")
                    output.append(f"Global compliance (informational): {'YES' if result.quality_metrics.is_valid else 'NO'}")
                    
                    if result.generation_stats:
                        stats = result.generation_stats
                        output.append(f"Backtracking attempts: {stats.get('backtrack_count', 0)}")
                        output.append(f"Total attempts: {stats.get('total_attempts', 0)}")
                        roll = stats.get('window_rollup', {})
                        gc_min = roll.get('gc_min')
                        gc_max = roll.get('gc_max')
                        tm_min = roll.get('tm_min')
                        tm_max = roll.get('tm_max')
                        hp_max = roll.get('hairpin_tm_max')
                        hd_max = roll.get('homodimer_tm_max')
                        # Display if available
                        roll_lines = []
                        if gc_min is not None and gc_max is not None:
                            roll_lines.append(f"Window GC: {gc_min:.2%} - {gc_max:.2%}")
                        if tm_min is not None and tm_max is not None:
                            roll_lines.append(f"Window Tm: {tm_min:.1f} - {tm_max:.1f}°C")
                        if hp_max is not None:
                            roll_lines.append(f"Window max hairpin Tm: {hp_max:.1f}°C")
                        if hd_max is not None:
                            roll_lines.append(f"Window max homodimer Tm: {hd_max:.1f}°C")
                        if roll_lines:
                            output.append("; ".join(roll_lines))
                else:
                    output.append(f"Success: NO")
                    output.append(f"Error: {result.error_message}")
                    if result.generation_time:
                        output.append(f"Generation time: {result.generation_time:.2f}s")
                
                output.append("")  # Empty line
            
            return "\n".join(output)


def main():
    """Main entry function."""
    parser = create_parser()
    args = parser.parse_args()

    # Show primer3 status if requested
    if getattr(args, 'show_primer3_status', False):
        print("\nPrimer3 Status:")
        print(f"  Available: {'YES' if PRIMER3_AVAILABLE else 'NO'}")
        if not PRIMER3_AVAILABLE:
            print("  Mode: Fallback calculations (simplified)")
            print("  Note: Install primer3-py for full thermodynamic calculations")
        else:
            print("  Mode: Full thermodynamic calculations")
        print()

    try:
        # Configure logging on CLI side (the library does not configure global logging)
        log_level = logging.DEBUG if args.verbose else (logging.ERROR if args.quiet else logging.INFO)
        logging.basicConfig(level=log_level, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # Ensure package loggers respect the chosen level
        logging.getLogger('dna_generator').setLevel(log_level)

        # Tworzenie konfiguracji
        config = create_config_from_args(args)
        
        # Show warning if thermodynamic checks are enabled but primer3 is not available
        if not PRIMER3_AVAILABLE and args.verbose:
            thermo_status = config.get_thermodynamic_status()
            if any([thermo_status['melting_temperature_enabled'],
                   thermo_status['hairpin_structures_enabled'],
                   thermo_status['homodimer_structures_enabled']]):
                print("\nNote: Using fallback thermodynamic calculations (primer3 not available)\n")

        # Inicjalizacja generatora
        generator = DNAGenerator(config)

        # Generowanie sekwencji
        if args.count == 1:
            results = [generator.generate(args.initial, args.length, args.seed)]
        else:
            results = generator.generate_multiple(args.initial, args.length, args.count, args.seed)

        # Export CSV if requested (only for successful sequences)
        if args.csv_file:
            successful_results = [r for r in results if r.success]
            if successful_results:
                # Export first successful sequence
                seq = successful_results[0].sequence
                export_windows_to_csv(
                    sequence=seq,
                    window_size=config.window_size,
                    validator=generator.validator,
                    csv_path=args.csv_file
                )
                if len(successful_results) > 1 and not args.quiet:
                    print(f"Note: CSV export uses first successful sequence only ({len(successful_results)} total)\n")
            else:
                print(f"Warning: No successful sequences to export to CSV", file=sys.stderr)

        # Formatowanie wyjścia
        output = format_output(
            results,
            args.format,
            args.initial,
            sequences_only=args.sequences_only,
            json_raw_metrics=args.json_raw_metrics,
        )
        
        # Write results
        if args.output:
            with open(args.output, 'w') as f:
                f.write(output)
            print(f"Results written to file: {args.output}")
        else:
            print(output)
        
        # Summary
        successful = sum(1 for r in results if r.success)
        if not args.quiet:
            print(f"\nSummary: {successful}/{len(results)} sequences generated successfully")
            
            if config.generation_mode == GenerationMode.DETERMINISTIC and len(results) > 1:
                sequences = [r.sequence for r in results if r.success]
                if sequences:
                    all_identical = all(seq == sequences[0] for seq in sequences)
                    print(f"All sequences identical: {all_identical}")
        
        # Exit code
        return 0 if successful > 0 else 1
        
    except GeneratorError as e:
        print(f"Generator error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
