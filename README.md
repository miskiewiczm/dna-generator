# DNA Generator

[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Deterministic DNA sequence generator with biochemical quality control. Generates high-quality DNA sequences using a backtracking algorithm with configurable validation profiles for PCR, primer design, and general molecular biology applications.

## Features

- **Deterministic Mode**: Reproducible sequences with seed-based generation
- **Random Mode**: Generate diverse sequence libraries
- **Quality Control**: Validation of GC content, melting temperature, secondary structures, homopolymer runs
- **Backtracking Algorithm**: Intelligent generation that meets biochemical constraints
- **Validation Profiles**: Pre-configured profiles (strict, relaxed, PCR-friendly) or custom JSON profiles
- **CLI Interface**: Command-line tool for batch generation
- **Built on dna-commons**: Uses shared validation and thermodynamic calculation library

## Installation

```bash
pip install dna-generator
```

For development:

```bash
git clone https://github.com/miskiewiczm/dna-generator.git
cd dna-generator
pip install -e .[dev]
```

## Quick Start

### Basic Usage

```python
from dna_generator import DNAGenerator, GeneratorConfig, GenerationMode

# Create generator with default configuration
config = GeneratorConfig(
    validation_profile='sequence_only',
    generation_mode=GenerationMode.DETERMINISTIC,
    seed=42
)

generator = DNAGenerator(config)

# Generate a 200bp sequence
result = generator.generate("ATGCGATCGTAGC", 200)

if result.success:
    print(f"Generated sequence: {result.sequence}")
    print(f"GC Content: {result.quality_metrics.gc_content:.2%}")
    print(f"Melting Temp: {result.quality_metrics.melting_temperature:.1f}°C")
else:
    print(f"Generation failed: {result.error_message}")
```

### Generate Multiple Sequences (Library)

```python
# Random mode for diverse sequences
config = GeneratorConfig(
    validation_profile='sequence_only',
    generation_mode=GenerationMode.RANDOM
)

generator = DNAGenerator(config)

# Generate 10 diverse sequences
results = generator.generate_multiple("ATGCGATCGTAGC", 150, count=10)

sequences = [r.sequence for r in results if r.success]
print(f"Generated {len(sequences)} sequences")
```

### Command Line Usage

```bash
# Generate a single sequence
dna-generator --initial ATGCGATCGT --length 200 --profile sequence_only --seed 42

# Generate multiple sequences
dna-generator --initial ATGCGATCGT --length 200 --count 10 --mode random

# Use custom validation profile
dna-generator --initial ATGCGATCGT --length 200 --profile pcr_friendly

# Save to FASTA file
dna-generator --initial ATGCGATCGT --length 200 --count 5 --format fasta --output sequences.fasta
```

## Validation Profiles

DNA Generator includes several pre-configured validation profiles:

- **`sequence_only`**: Basic sequence checks (homopolymers, dinucleotide repeats) without thermodynamic validation
- **`relaxed`**: Wider GC/Tm tolerances for general-purpose use
- **`pcr_friendly`**: Balanced constraints optimized for PCR applications
- **`strict`**: Narrow GC/Tm ranges for high-quality sequences
- **`none`**: No validation (random generation)

### Custom Profiles

Create custom validation profiles in JSON format:

```json
{
  "profiles": {
    "my_custom_profile": {
      "description": "Custom validation rules",
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
```

Save as `user_profiles.json` and use with `--profile-file user_profiles.json`.

## Examples

The `examples/` directory contains comprehensive usage examples:

- `basic_generation.py`: Simple sequence generation
- `primer_design.py`: PCR primer design workflow
- `sequence_library.py`: Generate diverse sequence libraries
- `custom_validation.py`: Using custom validation profiles

Run examples:

```bash
python examples/basic_generation.py
python examples/primer_design.py
```

## Documentation

For detailed documentation, see [DOCUMENTATION.md](DOCUMENTATION.md).

Topics covered:
- Algorithm details
- Configuration options
- Validation rules and parameters
- API reference
- Troubleshooting

## Testing

```bash
# Run tests with coverage
pytest tests/ --cov=dna_generator --cov-report=term-missing

# Run specific test file
pytest tests/test_generator.py -v
```

## Dependencies

- **dna-commons** (>=0.1.0): Shared validation and thermodynamic calculations
- **primer3-py** (>=2.0.0): Thermodynamic calculations (via dna-commons)

## Citation

If you use DNA Generator in your research, please cite:

```bibtex
@software{dna_generator,
  title = {DNA Generator: Deterministic DNA Sequence Generator},
  author = {[Author Name]},
  year = {2025},
  url = {https://github.com/miskiewiczm/dna-generator},
  version = {1.0.0}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Related Projects

- [dna-commons](https://github.com/miskiewiczm/dna-commons): Shared DNA analysis and validation library

## Contact

For questions, issues, or suggestions, please open an issue on GitHub.
