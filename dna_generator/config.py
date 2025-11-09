"""
Configuration for the DNA sequence generator.

Holds all configuration and parameters used by the generator, including
primer3 settings, sequence quality limits, and generation modes.

Validation profiles are loaded from JSON files:
- default_profiles.json (shipped with package)
- user_profiles.json (optional, in current directory)
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Dict, Any
import logging

logger = logging.getLogger(__name__)


class GenerationMode(Enum):
    """DNA sequence generation mode."""
    DETERMINISTIC = "deterministic"  # Deterministyczny z seedem
    RANDOM = "random"  # Czysto losowy


@dataclass
class GeneratorConfig:
    """
    Configuration for the DNA sequence generator.
    
    Attributes:
        generation_mode: Generation mode (deterministic/random)
        seed: Seed for the generator (None = auto in deterministic mode)
        
        primer3_mv_conc: Monovalent ion concentration [mM]
        primer3_dv_conc: Divalent ion concentration [mM]
        primer3_dntp_conc: dNTP concentration [mM]
        primer3_dna_conc: DNA concentration [nM]
        
        min_gc: Minimum GC content (0.0–1.0)
        max_gc: Maximum GC content (0.0–1.0)
        min_tm: Minimum melting temperature [°C]
        max_tm: Maximum melting temperature [°C]
        max_hairpin_tm: Maximum hairpin Tm [°C]
        max_homodimer_tm: Maximum homodimer Tm [°C]
        
        window_size: Analysis window size for quality checks
        max_homopolymer_length: Maximum length of identical nucleotide runs
        max_dinucleotide_repeats: Maximum count of dinucleotide repeats
        max_3prime_gc: Maximum G/C count in the last 5 nucleotides
        
        max_backtrack_attempts: Maximum number of backtracking attempts
        enable_progress_logging: Enable progress logging
        log_level: Logging level
    """
    
    # Tryb generowania
    generation_mode: GenerationMode = GenerationMode.DETERMINISTIC
    seed: Optional[int] = None
    
    # Parametry primer3 dla obliczeń termodynamicznych
    primer3_mv_conc: float = 50.0
    primer3_dv_conc: float = 4.0
    primer3_dntp_conc: float = 0.5
    primer3_dna_conc: float = 50.0
    
    # Ograniczenia zawartości GC
    min_gc: float = 0.45
    max_gc: float = 0.55
    
    # Ograniczenia temperatury topnienia
    min_tm: float = 55.0
    max_tm: float = 65.0
    
    # Ograniczenia struktur wtórnych
    max_hairpin_tm: float = 30.0
    max_homodimer_tm: float = 30.0
    
    # Parametry jakości sekwencji
    window_size: int = 20
    max_homopolymer_length: int = 4
    max_dinucleotide_repeats: int = 2
    max_3prime_gc: int = 3
    
    # Parametry algorytmu
    max_backtrack_attempts: int = 10000
    enable_progress_logging: bool = True

    # Logowanie (konfiguracja logowania nie jest wykonywana przez bibliotekę)
    log_level: str = "INFO"
    configure_logging: bool = False

    # Heurystyki backtrackingu
    enable_backtrack_heuristics: bool = True
    
    # Flagi włączania/wyłączania poszczególnych sprawdzeń
    validation_profile: Optional[str] = None  # Nazwa profilu lub None dla niestandardowego
    validation_rules: Dict[str, bool] = field(default_factory=lambda: {
        'gc_content': True,
        'melting_temperature': True,
        'hairpin_structures': True,
        'homodimer_structures': True,
        'homopolymer_runs': True,
        'dinucleotide_repeats': True,
        'three_prime_stability': True
    })
    
    def __post_init__(self):
        """Validate and initialize after instance creation."""
        # Jeśli podano profil, zastosuj jego ustawienia (reguły + zalecane progi)
        if self.validation_profile:
            # Load profiles from JSON files
            from .profile_loader import get_profile_loader

            loader = get_profile_loader()
            profile = loader.get_profile(self.validation_profile)

            if profile:
                # Reguły
                rules = profile.get('rules')
                if rules:
                    self.validation_rules = rules.copy()
                # Parametry
                params = profile.get('params', {})
                for k, v in params.items():
                    if hasattr(self, k):
                        setattr(self, k, v)
                logger.debug(f"Applied validation profile: {self.validation_profile}")

                # Auto-disable heuristics for 'none' profile to ensure truly random generation
                if self.validation_profile == 'none':
                    self.enable_backtrack_heuristics = False
                    logger.debug("Auto-disabled heuristics for 'none' profile")
            else:
                available = ', '.join(loader.list_profiles().keys())
                raise ValueError(f"Unknown validation profile: {self.validation_profile}. Available: {available}")

        # Auto-disable thermodynamic checks if primer3 is not available
        # Import here to avoid circular import
        from dna_commons import PRIMER3_AVAILABLE

        if not PRIMER3_AVAILABLE:
            thermo_rules = ['melting_temperature', 'hairpin_structures', 'homodimer_structures']
            disabled_rules = []
            for rule in thermo_rules:
                if self.validation_rules.get(rule, False):
                    self.validation_rules[rule] = False
                    disabled_rules.append(rule)

            if disabled_rules:
                logger.warning(
                    f"primer3 not available - automatically disabled thermodynamic checks: {', '.join(disabled_rules)}. "
                    f"Install primer3-py to enable these checks or use 'sequence_only' profile."
                )

        self._validate_config()
        # Nie konfigurujemy globalnego logowania w bibliotece.
        # Jeśli naprawdę wymagane, można ustawić configure_logging=True,
        # jednak zalecane jest konfigurowanie logowania w aplikacji/CLI.
        if getattr(self, 'configure_logging', False):
            self._setup_logging()
    
    def _validate_config(self) -> None:
        """
        Validate configuration parameters.
        
        Raises:
            ValueError: If parameters have invalid values
        """
        # GC content validation
        if not (0.0 <= self.min_gc <= 1.0):
            raise ValueError(f"min_gc must be in [0.0, 1.0], got {self.min_gc}")
        if not (0.0 <= self.max_gc <= 1.0):
            raise ValueError(f"max_gc must be in [0.0, 1.0], got {self.max_gc}")
        if self.min_gc > self.max_gc:
            raise ValueError(f"min_gc ({self.min_gc}) cannot be greater than max_gc ({self.max_gc})")
        
        # Melting temperature validation
        if self.min_tm < 0:
            raise ValueError(f"min_tm cannot be negative, got {self.min_tm}")
        if self.min_tm > self.max_tm:
            raise ValueError(f"min_tm ({self.min_tm}) cannot be greater than max_tm ({self.max_tm})")
        
        # primer3 concentrations validation
        if self.primer3_mv_conc <= 0:
            raise ValueError(f"primer3_mv_conc must be positive, got {self.primer3_mv_conc}")
        if self.primer3_dv_conc < 0:
            raise ValueError(f"primer3_dv_conc cannot be negative, got {self.primer3_dv_conc}")
        if self.primer3_dntp_conc < 0:
            raise ValueError(f"primer3_dntp_conc cannot be negative, got {self.primer3_dntp_conc}")
        if self.primer3_dna_conc <= 0:
            raise ValueError(f"primer3_dna_conc must be positive, got {self.primer3_dna_conc}")
        
        # Quality parameters validation
        if self.window_size < 5:
            raise ValueError(f"window_size must be >= 5, got {self.window_size}")
        if self.max_homopolymer_length < 2:
            raise ValueError(f"max_homopolymer_length must be >= 2, got {self.max_homopolymer_length}")
        if self.max_3prime_gc < 0 or self.max_3prime_gc > 5:
            raise ValueError(f"max_3prime_gc must be in [0, 5], got {self.max_3prime_gc}")
        
        # Algorithm parameters validation
        if self.max_backtrack_attempts < 1:
            raise ValueError(f"max_backtrack_attempts must be >= 1, got {self.max_backtrack_attempts}")
    
    def _setup_logging(self) -> None:
        """Konfiguruje logowanie zgodnie z ustawieniami."""
        numeric_level = getattr(logging, self.log_level.upper(), None)
        if not isinstance(numeric_level, int):
            raise ValueError(f"Invalid log level: {self.log_level}")
        
        logging.basicConfig(
            level=numeric_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
    
    def get_primer3_params(self) -> Dict[str, float]:
        """
        Return primer3 parameters as a dictionary.
        
        Returns:
            Dict of primer3 parameters
        """
        return {
            "mv_conc": self.primer3_mv_conc,
            "dv_conc": self.primer3_dv_conc,
            "dntp_conc": self.primer3_dntp_conc,
            "dna_conc": self.primer3_dna_conc
        }
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert configuration to a dictionary.
        
        Returns:
            Dict with configuration parameters
        """
        return {
            "generation_mode": self.generation_mode.value,
            "seed": self.seed,
            "primer3_mv_conc": self.primer3_mv_conc,
            "primer3_dv_conc": self.primer3_dv_conc,
            "primer3_dntp_conc": self.primer3_dntp_conc,
            "primer3_dna_conc": self.primer3_dna_conc,
            "min_gc": self.min_gc,
            "max_gc": self.max_gc,
            "min_tm": self.min_tm,
            "max_tm": self.max_tm,
            "max_hairpin_tm": self.max_hairpin_tm,
            "max_homodimer_tm": self.max_homodimer_tm,
            "window_size": self.window_size,
            "max_homopolymer_length": self.max_homopolymer_length,
            "max_dinucleotide_repeats": self.max_dinucleotide_repeats,
            "max_3prime_gc": self.max_3prime_gc,
            "max_backtrack_attempts": self.max_backtrack_attempts,
            "enable_progress_logging": self.enable_progress_logging,
            "log_level": self.log_level
        }
    
    def is_primer3_available(self) -> bool:
        """Check if primer3 is available for thermodynamic calculations."""
        from dna_commons import PRIMER3_AVAILABLE
        return PRIMER3_AVAILABLE

    def get_thermodynamic_status(self) -> Dict[str, Any]:
        """Get status of thermodynamic calculation capabilities."""
        from dna_commons import PRIMER3_AVAILABLE
        return {
            'primer3_available': PRIMER3_AVAILABLE,
            'melting_temperature_enabled': self.validation_rules.get('melting_temperature', False),
            'hairpin_structures_enabled': self.validation_rules.get('hairpin_structures', False),
            'homodimer_structures_enabled': self.validation_rules.get('homodimer_structures', False),
            'fallback_mode': not PRIMER3_AVAILABLE
        }

    def __str__(self) -> str:
        """String representation of the configuration."""
        from dna_commons import PRIMER3_AVAILABLE
        primer3_status = "primer3" if PRIMER3_AVAILABLE else "fallback"
        return (
            f"GeneratorConfig(mode={self.generation_mode.value}, "
            f"GC={self.min_gc:.0%}-{self.max_gc:.0%}, "
            f"Tm={self.min_tm:.0f}-{self.max_tm:.0f}°C, "
            f"window={self.window_size}bp, "
            f"thermo={primer3_status})"
        )


# Domyślna konfiguracja
DEFAULT_CONFIG = GeneratorConfig()
